# PsyOPS 

PsyOPS (short for **Psy**chiatric **O**mnilocus **P**rioritization **S**core) 
is a software package for prioritizing causal genes from psychiatric genome-wide 
association studies (GWAS). 

For a detailed description of PsyOPS, check out [our paper](https://www.nature.com/articles/s41380-022-01542-6)!

## Getting started

Download this repository, which contains the PsyOPS code (`PsyOPS.py`) and all 
necessary data files:

```  
git clone https://github.com/Wainberg/PsyOPS
cd PsyOPS
```

PsyOPS requires Python >= 3.6 with the `numpy`, `scipy`, `pandas`, and 
`scikit-learn` packages. You can install Python and the required packages with 
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) as follows:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -bp ~/.anaconda3
rm Miniconda3-latest-Linux-x86_64.sh
export PATH="~/.anaconda3/bin:$PATH"
conda install -y numpy scipy pandas scikit-learn
```

## Running PsyOPS

PsyOPS.py takes two arguments, both of which are mandatory:

1) `--GWAS-hit-file`: a tab-separated file of lead variants. Must contain 3 
   columns: "rs" (the variant's rs number), "chrom" (the variant's chromosome; 
   only autosomal variants will be analyzed), and "bp_hg19" (the variant's hg19/
   GRCh37 base pair coordinates). 
2) `--output-file`: a tab-separated file where results will be output. 
   Each row lists a gene within 500 kilobases of one of the lead variants. 
   Columns denote the lead variant, the distance to the gene, the gene, the 
   gene's chromosome, transcription start site and transcription end site (in 
   hg19/GRCh37 coordinates), whether the gene has each of the three PsyOPS 
   features (extreme pLI, brain-enriched expression, neurodevelopmental 
   disorder), and the gene's PsyOPS score.

As an example, let's run PsyOPS on the 102 lead variants from 
[a recent major depressive disorder (MDD) GWAS
](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6522363), provided in the file 
`Howard_102_MDD.tsv` in this repository:

```
python PsyOPS.py --GWAS-hit-file Howard_102_MDD.tsv --output-file PsyOPS_MDD.tsv
```

This takes about 5 seconds to run and prints the following output:

```
Loaded 102 autosomal GWAS hits from Howard_102_MDD.tsv

88 of 708 genes (12.4%) are positives (i.e. nearest genes)

Model coefficients (when trained without cross-validation):
                                   OR  lower_CI   upper_CI             p
Intercept                    0.071024  0.051358   0.098221  1.531120e-57
Extreme pLI                  2.710843  1.491015   4.928634  1.076882e-03
Brain-enriched expression    2.927626  1.661044   5.160004  2.033198e-04
Neurodevelopmental disorder  5.458068  2.909110  10.240418  1.249739e-07

WARNING: only one class represented for chromosome chr15; excluding from accuracy calculation

Average performance across chromosomes: AUC 0.71 +/- 0.05, AUPRC 0.53 +/- 0.07

Prioritized genes:
ADARB2, AGBL4, ASTN2, ATP2A2, BRINP2, CABP1, CACNA1D, CACNA1E, CDC42BPB, CELF4, CHD6, CNTLN, CRB1, CRYBA1, CTNND1, CTTNBP2, CYP7B1, DAAM1, DAGLA, DCC, DDX6, DNAH3, DOCK9, DPF2, DRD2, EIF3B, ELAVL2, EP300, ESR2, EYS, FHIT, GRIK2, GRIK3, GRM5, H4C12, HIVEP2, KCNG2, KIF2A, KIRREL3, KLF7, LHX2, LIN28B, LRFN5, MEF2C, MEIS2, MPPED2, MRTFB, MSANTD1, NCOA2, NEGR1, NR2F1, NR4A2, PAX5, PCDH8, PCDH9, PCLO, PHOX2B, PPFIA1, PROX2, QRICH1, RAB3B, RBFOX1, REEP1, RERE, RFX3, RSRC1, SEMA6D, SGIP1, SHISA9, SLC12A5, SORCS3, SOX5, SPRY2, STXBP5, TCF4, TEX26, THSD7A, TNR, TRPC3, TUSC1, VRK2, WDR17, ZNF445, ZNF536

Exported results to PsyOPS_MDD.tsv
```

We also provide lead variant files for the other two psychiatric GWAS analyzed 
in the PsyOPS paper, for [bipolar disorder
](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8192451) (`Mullins_64_BIP.tsv`) 
and [schizophrenia](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692) 
(`Pardinas_145_SCZ.tsv`).
