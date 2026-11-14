# Molecular Subtypes of Breast Cancer: Discovery, Outcomes and Prediction from TCGA Gene Expression

An end-to-end Python analysis of ~1,100 breast tumours from The Cancer Genome Atlas (TCGA): discovering molecular subtypes without labels, testing whether they predict patient outcomes, and asking whether gene expression improves risk prediction beyond standard clinical variables.

## Motivation

Patients with the same clinical diagnosis can have very different diseases at the molecular level, and so respond differently to treatment and have different outcomes. Breast cancer is the classic example: the PAM50 gene-expression subtypes (Luminal A, Luminal B, HER2-enriched, Basal-like, Normal-like) are now used alongside clinical stage to guide care.

This project uses breast cancer as a test case for three general questions in predictive medicine:

1. **Discovery.** Can we recover meaningful disease subtypes from molecular data alone, without being told the answer?
2. **Clinical relevance.** Do these molecular subtypes differ in patient outcomes, independently of age and stage?
3. **Prediction.** Does molecular information actually improve outcome prediction over clinical variables, when evaluated honestly with cross-validation?

## Data

All data are public, from the [UCSC Xena](https://xenabrowser.net) TCGA Hub (cohort: TCGA Breast Cancer, BRCA):

| File | Contents |
|------|----------|
| `HiSeqV2` | RNA-seq gene expression, log2(normalised count + 1), ~20,000 genes |
| `BRCA_clinicalMatrix` | Clinical data including PAM50 subtype, age and pathologic stage |
| `BRCA_survival.txt` | Curated survival endpoints from the TCGA Pan-Cancer Clinical Data Resource ([Liu et al. 2018](https://doi.org/10.1016/j.cell.2018.02.052)) |

One primary tumour per patient is kept. The main endpoint is the **progression-free interval (PFI)**, which Liu et al. recommend for breast cancer because overall survival has relatively few events in TCGA's follow-up period; overall survival (OS) is reported as a secondary endpoint.

## Methods

| Step | Script | What it does |
|------|--------|--------------|
| 1 | `01_download_data.py` | Download the three files from UCSC Xena |
| 2 | `02_preprocess.py` | Keep primary tumours, filter lowly expressed genes, harmonise clinical variables, merge survival data |
| 3 | `03_unsupervised_subtypes.py` | Top 2,000 variable genes → PCA → **consensus clustering** (k-means on 50 random 80% subsamples, k = 2–8). k chosen by the proportion of ambiguous clustering (PAC). Clusters compared with PAM50 using ARI and NMI |
| 4 | `04_survival_analysis.py` | Kaplan–Meier curves and log-rank tests by cluster and by PAM50; **multivariable Cox model** adjusting for age and stage |
| 5 | `05_predict_subtype.py` | **Nested cross-validation** (5 outer × 3 inner folds) comparing elastic-net logistic regression, random forest and a baseline for predicting PAM50 subtype; inspection of which genes the sparse model selects and whether known PAM50 genes are enriched (hypergeometric test) |
| 6 | `06_survival_prediction.py` | Penalised Cox models with five feature sets (clinical; clinical + PAM50; clinical + discovered clusters; clinical + expression PCs; expression only), compared by **C-index over 5 × 5-fold repeated cross-validation** |

**Avoiding data leakage.** In steps 5 and 6, gene selection, scaling and PCA are fitted on the training fold only, so the test fold never influences model building. Without this, cross-validated performance would be optimistically biased.


## How to reproduce

```bash
git clone https://github.com/abdullahikrambutt/brca-molecular-subtypes.git
cd brca-molecular-subtypes

python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python run_all.py --test           # optional: 3-minute check on synthetic data
python run_all.py                  # full analysis (~10-20 minutes on a laptop)
```

Figures go to `figures/` and tables to `results/`. Settings such as the number of genes, the range of k, or a fixed choice of k are in `scripts/config.py`.

## Project structure

```
├── scripts/
│   ├── config.py                     # paths and settings
│   ├── 01_download_data.py
│   ├── 02_preprocess.py
│   ├── 03_unsupervised_subtypes.py
│   ├── 04_survival_analysis.py
│   ├── 05_predict_subtype.py
│   └── 06_survival_prediction.py
├── tests/
│   └── make_synthetic_data.py        # synthetic data for the --test run
├── figures/
├── results/
├── run_all.py
└── requirements.txt
```

## Limitations

- **Circularity in subtype prediction.** PAM50 labels were derived from gene expression, so predicting them from expression is expected to work well. The interest lies in model sparsity and interpretability, not raw accuracy.
- **Single omics layer.** Only gene expression is used. Real molecular endotyping would integrate genomics, methylation, proteomics and other layers.
- **Cohort and follow-up.** TCGA is not a population-based cohort, follow-up is relatively short for breast cancer, and treatment information is incomplete, so outcome associations may be confounded.
- **No external validation.** All performance estimates come from cross-validation within TCGA; a proper evaluation would test models on an independent cohort such as METABRIC.

## Possible extensions

- Add copy-number or methylation data and compare single-omic and multi-omic clustering
- Validate the discovered subtypes and risk models in METABRIC
- Replace PCA with an autoencoder or other representation-learning approach

## Resources

- Cancer Genome Atlas Network (2012). [Comprehensive molecular portraits of human breast tumours](https://doi.org/10.1038/nature11412). *Nature*
- Parker et al. (2009). [Supervised risk predictor of breast cancer based on intrinsic subtypes](https://doi.org/10.1200/JCO.2008.18.1370). *J Clin Oncol* (PAM50)
- Monti et al. (2003). Consensus clustering. *Machine Learning*
- Şenbabaoğlu et al. (2014). [Critical limitations of consensus clustering in class discovery](https://doi.org/10.1038/srep06207). *Sci Rep* (PAC)
- Liu et al. (2018). [An integrated TCGA pan-cancer clinical data resource](https://doi.org/10.1016/j.cell.2018.02.052). *Cell*

Data: The results shown here are based on data generated by the TCGA Research Network: https://www.cancer.gov/tcga

## License

MIT
