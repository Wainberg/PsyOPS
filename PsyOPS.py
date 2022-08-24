import numpy as np, os, pandas as pd
from argparse import ArgumentParser
from scipy.stats import sem
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from subprocess import run
from scipy.stats import norm

# Dependencies: Python >= 3.6, numpy, scipy, pandas, scikit-learn

# Utility functions

def auPRC(Y_test, predictions):
    precision, recall, thresholds = precision_recall_curve(Y_test, predictions)
    auPRC = auc(recall, precision)
    return auPRC

def logistic_regression(X, y, return_model=False, report_intercept=False,
                        **kwargs):
    # Based on https://pingouin-stats.org/_modules/pingouin/regression.html
    # and https://pingouin-stats.org/generated/pingouin.logistic_regression.html
    assert isinstance(X, pd.DataFrame), type(X)
    assert isinstance(y, pd.Series), type(y)
    assert np.issubdtype(X.values.dtype, np.number), X.dtypes
    assert y.dtype == bool, y.dtype
    model = LogisticRegression(penalty='none', solver='newton-cg',
                               random_state=0, **kwargs).fit(X, y)
    predictions = model.decision_function(X)
    denom = np.tile(2 * (1 + np.cosh(predictions)), (X.shape[1] + 1, 1)).T
    X_design = np.column_stack((np.ones(len(X)), X))
    variance_covariance_matrix = np.linalg.pinv((X_design / denom).T @ X_design)
    index = X.columns.insert(0, 'Intercept') if report_intercept else X.columns
    betas = pd.Series(np.concatenate([model.intercept_, model.coef_[0]])
                      if report_intercept else model.coef_[0], index=index)
    ORs = np.exp(betas)
    sigmas = pd.Series(np.sqrt(np.diag(variance_covariance_matrix)[
                               0 if report_intercept else 1:]), index=index)
    lower_CIs = np.exp(betas - norm.isf(0.025) * sigmas)
    upper_CIs = np.exp(betas + norm.isf(0.025) * sigmas)
    pvalues = pd.Series(2 * norm.sf(np.abs(betas / sigmas)), index=index)
    results = pd.DataFrame({'OR': ORs, 'lower_CI': lower_CIs,
                            'upper_CI': upper_CIs, 'p': pvalues})
    return (results, model) if return_model else results

# Two arguments:
# 1) --GWAS-hit-file: TSV file of GWAS hits with columns rs, chrom, bp_hg19
# 2) --output-file: TSV file where results will be output

parser = ArgumentParser()
parser.add_argument('--GWAS-hit-file', required=True)
parser.add_argument('--output-file', required=True)
args = parser.parse_args()
GWAS_hits = pd.read_table(args.GWAS_hit_file, index_col='rs', 
                          dtype={'chrom': str})\
    .assign(chrom=lambda df: df.chrom.where(df.chrom.str.startswith('chr'),
                                            'chr' + df.chrom))\
    .query('chrom != "chrX" and chrom != "chrY" and chrom != "chrM"')\
    .assign(chrom_int=lambda df: df.chrom.str[3:].astype(int))\
    .sort_values(['chrom_int', 'bp_hg19'])\
    .drop('chrom_int', axis=1)
print(f'Loaded {len(GWAS_hits)} autosomal GWAS hits from {args.GWAS_hit_file}')

# Get universe of genes - autosomal Gencode coding genes with at least one
# non-readthrough transcript - and their hg19 coordinates

if not os.path.exists('non_readthrough_genes.txt'):
    run('curl -s ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/'
        'GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz | zcat | '
        'tr -d ";\\"" | awk \'$3 == "transcript" && $0 !~ '
        '/readthrough_transcript/ {print $15 == "gene_name" ? $16 : $17 == '
        '"gene_name" ? $18 : "ERROR"}\' | sort -u > non_readthrough_genes.txt',
        check=True, shell=True)

non_readthrough_genes = pd.read_table('non_readthrough_genes.txt',
                                      header=None, index_col=0).index

if not os.path.exists('coding_genes.bed'):
    run('curl -s ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/'
        'GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz | zcat | '
        'tr -d ";\\"" | awk \'$3 == "gene" && $0 ~ /protein_coding/ {print $1, '
        '$4, $5, $13 == "gene_name" ? $14 : $15 == "gene_name" ? $16 : '
        '"ERROR"}\' OFS="\t" | sort -k1,1V -k2,3n > coding_genes.bed',
        check=True, shell=True)

genes = pd.read_table('coding_genes.bed', header=None,
                      names=['chrom', 'start', 'end', 'gene'])\
    .groupby('gene')\
    .agg(dict(chrom='first', start='min', end='max'))\
    .reset_index()
assert len(genes) == 20084, len(genes)
non_autosomal_and_readthrough_query = \
    'chrom != "chrX" and chrom != "chrY" and chrom != "chrM" and ' \
    'gene in @non_readthrough_genes'
non_autosomal_and_readthrough_genes = \
    genes.query(f'not ({non_autosomal_and_readthrough_query})').gene
assert len(non_autosomal_and_readthrough_genes) == 1500, \
    len(non_autosomal_and_readthrough_genes)
genes = genes.query(non_autosomal_and_readthrough_query)
assert len(genes) == 18584, len(genes)

# Load gene aliases, and manually add a few that are missing

gene_alias_file = 'gene_aliases.tsv'
if not os.path.exists(gene_alias_file):
    run(f'curl -s ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/'
        f'locus_groups/protein-coding_gene.txt | cut -f2,11 | tail -n +2 > '
        f'{gene_alias_file}', check=True, shell=True)

alias_to_gene = pd.read_table(gene_alias_file, header=None)\
    .dropna()\
    .apply('|'.join, axis=1)\
    .str.split('|')\
    .map(set)
alias_to_gene = pd.Series({gene: aliases.difference({gene})
                           for aliases in alias_to_gene.values
                           for gene in aliases}).map(tuple).explode()
alias_to_gene = alias_to_gene[alias_to_gene.isin(genes.gene)]
alias_to_gene['DKFZP761J1410'] = 'PRG4'  # pLI
alias_to_gene['RP11-1055B8.7'] = 'BAHCC1'  # pLI
alias_to_gene['AC012309.1'] = 'ZNF585B'  # HPA
alias_to_gene['AL136373.1'] = 'TMEM275'  # HPA
alias_to_gene['H3.Y'] = 'H3Y1'  # HPA
alias_to_gene['SGK494'] = 'RSKR'  # HPA
alias_to_gene['TMEM155'] = 'SMIM43'  # HPA

def unalias(gene_set, return_non_matches=False, target_gene_set=genes.gene):
    matches = gene_set.intersection(target_gene_set)
    not_yet_matches = gene_set.difference(target_gene_set)
    remapped = alias_to_gene[alias_to_gene.index.isin(not_yet_matches)]
    assert not remapped.index.duplicated().any()
    new_matches = matches.union(pd.Index(remapped)).sort_values()
    assert new_matches.isin(target_gene_set).all()
    if return_non_matches:
        non_matches = not_yet_matches[
            ~not_yet_matches.isin(alias_to_gene.index)]
        return new_matches, non_matches
    else:
        return new_matches

# Load extreme-pLI genes (pLI > 0.99)

if not os.path.exists('constraint_scores.tsv'):
    run('curl -s storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/'
        'constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat | '
        'awk -F"\t" -v OFS="\t" \'$21 != "NA" && $68 == "protein_coding" '
        '{print $1, $21, $24, $30, $66, $75, $76, $77}\' | sort -u > '
        'constraint_scores.tsv', check=True, shell=True)

extreme_pLI_genes = pd.read_table('constraint_scores.tsv', header=None, names=[
    'gene', 'pLI', 'oe', 'LOEUF', 'cds_length', 'chrom', 'start', 'end'])\
    [['gene', 'pLI']]\
    .groupby('gene')\
    .agg(dict(pLI='max'))\
    .query('pLI >= 0.99')\
    .index
assert len(extreme_pLI_genes) == 1983, len(extreme_pLI_genes)
extreme_pLI_genes = extreme_pLI_genes[
    ~extreme_pLI_genes.isin(non_autosomal_and_readthrough_genes)]
assert len(extreme_pLI_genes) == 1841, len(extreme_pLI_genes)
extreme_pLI_genes = unalias(extreme_pLI_genes)
assert len(extreme_pLI_genes) == 1833, len(extreme_pLI_genes)

# Load brain-enriched expression genes from Human Protein Atlas
# https://www.proteinatlas.org/humanproteome/brain, red part of pie chart
# https://www.proteinatlas.org/humanproteome/tissue/tissue+specific for methods

brain_enriched_file = 'brain_enriched.txt'
if not os.path.exists(brain_enriched_file):
    run('curl -s https://www.proteinatlas.org/search/tissue_category_rna%3A'
        'brain%3BTissue+enriched%2CGroup+enriched%2CTissue+enhanced?'
        'format=tsv | cut -f1 | tail -n +2 > ' + brain_enriched_file,
        check=True, shell=True)

brain_enriched_genes = pd.read_table(
    brain_enriched_file, header=None, index_col=0).index
assert len(brain_enriched_genes) == 2587, len(brain_enriched_genes)
brain_enriched_genes = brain_enriched_genes[
    ~brain_enriched_genes.isin(non_autosomal_and_readthrough_genes)]
assert len(brain_enriched_genes) == 2441, len(brain_enriched_genes)
brain_enriched_genes = unalias(brain_enriched_genes)
assert len(brain_enriched_genes) == 2386, len(brain_enriched_genes)

# Load Genomics England NDD genes (monoallelic or unspecified inheritance)
# from https://panelapp.genomicsengland.co.uk/panels:
# Autism (v0.22, panelapp.genomicsengland.co.uk/panels/657)
# Intellectual disability (v3.2, panelapp.genomicsengland.co.uk/panels/285)
# Genetic epilepsy syndromes (v2.434, panelapp.genomicsengland.co.uk/panels/402)

get_GE_gene_set = lambda panel_file: pd.read_table(
    panel_file, usecols=['Entity Name', 'Entity type', 'Model_Of_Inheritance'],
    index_col='Entity Name')\
    .fillna('NA')\
    .query('(Model_Of_Inheritance.str.lower().str.contains("monoallelic") or '
           'Model_Of_Inheritance == "NA") and `Entity type` == "gene"',
           engine='python')\
    .index
GE_autism_genes = get_GE_gene_set('GE_autism.tsv')
GE_ID_genes = get_GE_gene_set('GE_ID.tsv')
GE_epilepsy_genes = get_GE_gene_set('GE_epilepsy.tsv')
assert len(GE_autism_genes) == 732, len(GE_autism_genes)
assert len(GE_ID_genes) == 977, len(GE_ID_genes)
assert len(GE_epilepsy_genes) == 279, len(GE_epilepsy_genes)
NDD_genes = GE_autism_genes.union(GE_ID_genes).union(GE_epilepsy_genes)
assert len(NDD_genes) == 1556, len(NDD_genes)
NDD_genes = NDD_genes[~NDD_genes.isin(non_autosomal_and_readthrough_genes)]
assert len(NDD_genes) == 1379, len(NDD_genes)
NDD_genes = unalias(NDD_genes)
assert len(NDD_genes) == 1370, len(NDD_genes)

# Get gene annotations

gene_annotations = {
    'Extreme pLI': extreme_pLI_genes,
    'Brain-enriched expression': brain_enriched_genes,
    'Neurodevelopmental disorder': NDD_genes}
assert len(gene_annotations) == 3, len(gene_annotations)
annotation_matrix = pd.DataFrame({
    annotation_name: pd.Index(genes.gene).to_series().isin(annotation)
    for annotation_name, annotation in gene_annotations.items()})
num_features_annotations = ({
    f'{n} feature{"s" if n > 1 else ""}':
        annotation_matrix.index[annotation_matrix.sum(axis=1) == n]
    for n in range(1, len(gene_annotations) + 1)})
assert len(num_features_annotations) == 3, len(num_features_annotations)

# Define features for all genes within 500 kb

nearest_genes, non_nearest_genes = pd.Index([]), pd.Index([])
locus_genes = pd.Series(index=GWAS_hits.index, dtype=float)
for rs, chrom, bp in GWAS_hits.itertuples():
    genes_within_500_kb = \
        genes[['gene', 'chrom', 'start', 'end']]\
        .query('chrom == @chrom')\
        .assign(distance=lambda df: (df[['start', 'end']] - bp)
                                    .abs()
                                    .min(axis=1)
                                    .where((bp < df.start) |
                                           (df.end < bp), 0),
                gene_length=lambda df: df.end - df.start)\
        .query('distance < 500_000')\
        .sort_values(['distance', 'gene_length'])\
        [['gene', 'distance']]\
        .set_index('gene')
    locus_genes[rs] = genes_within_500_kb.index
    nearest_genes = nearest_genes.union(genes_within_500_kb.iloc[:1].index)
    non_nearest_genes = non_nearest_genes.union(
        genes_within_500_kb.iloc[1:].index)

non_nearest_genes = non_nearest_genes.difference(nearest_genes)
y = pd.concat([pd.Series(True, index=nearest_genes),
               pd.Series(False, index=non_nearest_genes)]).sort_index()
X = pd.DataFrame({
    gene_set_name: pd.Series({gene: gene in gene_set for gene in y.index})
    for gene_set_name, gene_set in gene_annotations.items()}).astype(float)
print(f'\n{y.sum()} of {len(y)} genes ({100 * y.mean():.1f}%) are positives '
      f'(i.e. nearest genes)')

# Just for info, print multivariate odds ratios, without cross-validation

print('\nModel coefficients (when trained without cross-validation):')
print(logistic_regression(X, y, report_intercept=True))

# Using leave-one-chromosome-out cross-validation, get P(nearest gene),
# i.e. the PsyOPS score; also print model performance

PsyOPS_scores = pd.Series(index=X.index, dtype=float)
chroms = genes.set_index('gene').chrom[y.index]
metrics = {'AUC': roc_auc_score, 'AUPRC': auPRC}
performance = {metric: [] for metric in metrics}
for fold_index, (train_fold, test_fold) in enumerate(
        LeaveOneGroupOut().split(X, y, chroms), start=1):
    X_train = X.iloc[train_fold]
    Y_train = y.iloc[train_fold]
    X_test = X.iloc[test_fold]
    y_test = y.iloc[test_fold]
    scaler = StandardScaler()
    X_train = pd.DataFrame(scaler.fit_transform(X_train),
                           index=X_train.index, columns=X_train.columns)
    X_test = pd.DataFrame(scaler.transform(X_test),
                          index=X_test.index, columns=X_test.columns)
    results, model = logistic_regression(X_train, Y_train,
                                         return_model=True)
    PsyOPS_scores.iloc[test_fold] = y_pred = model.predict_proba(X_test)[:, 1]
    if y_test.nunique() > 1:
        for metric_name, metric in metrics.items():
            performance[metric_name].append(metric(y_test, y_pred))
    else:
        print(f'\nWARNING: only one class represented for chromosome '
              f'{chroms[test_fold].iloc[0]}; excluding from accuracy '
              f'calculation')

print(f'\nAverage performance across chromosomes: ' + ', '.join(
      f'{metric_name} {np.mean(performance[metric_name]):.2f} +/- '
      f'{sem(performance[metric_name]):.2f}' for metric_name in metrics))

# Print genes prioritized by PsyOPS (i.e. top-scoring at each locus, breaking
# ties by taking the first/nearest gene, which idxmax() does by default)

PsyOPS_genes = pd.Index({
    PsyOPS_scores.reindex(locus_genes[rs]).idxmax()
    for rs in GWAS_hits.index if len(locus_genes[rs]) > 0}).sort_values()
print(f'\nPrioritized genes:\n{", ".join(PsyOPS_genes)}')

# Export results as TSV

results = pd.concat([
    genes[['gene', 'chrom', 'start', 'end']]
        .query('chrom == @chrom')
        .assign(distance=lambda df: (df[['start', 'end']] - bp)
                                    .abs()
                                    .min(axis=1)
                                    .where((bp < df.start) | (df.end < bp), 0),
                gene_length=lambda df: df.end - df.start)
        .query('distance < 500_000')
        .sort_values(['distance', 'gene_length'])
        .assign(lead_variant=rs,
                extreme_pLI=lambda df: df.gene.isin(
                    gene_annotations['Extreme pLI']),
                brain_enriched_expression=lambda df: df.gene.isin(
                    gene_annotations['Brain-enriched expression']),
                neurodevelopmental_disorder=lambda df: df.gene.isin(
                    gene_annotations['Neurodevelopmental disorder']),
                PsyOPS_score=lambda df: df.gene.apply(
                    PsyOPS_scores.__getitem__))
        [['lead_variant', 'distance', 'gene', 'chrom', 'start', 
          'end', 'extreme_pLI', 'brain_enriched_expression', 
          'neurodevelopmental_disorder', 'PsyOPS_score']]
    for rs, chrom, bp in GWAS_hits.itertuples()])
results.to_csv(args.output_file, sep='\t', index=False)
print(f'\nExported results to {args.output_file}')
