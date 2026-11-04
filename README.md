# Evaluating EVEE on diabetes-related genes

## Question

Can EVEE distinguish harmful (pathogenic) from harmless (benign) variants across genes involved in insulin production and blood-sugar control, and does its performance differ between **missense** and **splice-site** variants?

We focus on five genes:

| Gene   | Role                                               | Associated conditions                      |
| ------ | -------------------------------------------------- | ------------------------------------------ |
| GCK    | Glucose sensor in pancreatic beta cells            | GCK-MODY (a mild, stable form of diabetes) |
| HNF1A  | Transcription factor regulating beta-cell function | HNF1A-MODY                                 |
| HNF4A  | Transcription factor regulating beta-cell function | HNF4A-MODY                                 |
| KCNJ11 | Subunit of the beta-cell K-ATP channel (Kir6.2)    | Neonatal diabetes, hyperinsulinism         |
| ABCC8  | Subunit of the beta-cell K-ATP channel (SUR1)      | Neonatal diabetes, hyperinsulinism         |

EVEE is _[one-line description of EVEE and what it predicts]_.

## Why this matters

- **Variant interpretation affects treatment.** In monogenic diabetes, the specific gene and variant often change management. For example, GCK-MODY usually needs no medication, HNF1A- and HNF4A-MODY typically respond well to sulfonylureas, and many people with K-ATP channel neonatal diabetes (KCNJ11, ABCC8) can switch from insulin to sulfonylureas. A correct pathogenic/benign call can therefore avoid misdiagnosis and unnecessary treatment.
- **Many variants remain unclassified.** A large share of variants in these genes in ClinVar are variants of uncertain significance (VUS) or have conflicting interpretations. Reliable computational predictors could help prioritize them for follow-up.
- **Missense and splice-site variants act through different mechanisms.** Missense variants change a single amino acid in the protein. Canonical splice-site variants disrupt how the transcript is processed, and their effect is not always visible at the protein level. Predictors may therefore perform differently on the two classes, and we want to know whether EVEE does.
- **These genes cover different disease mechanisms.** HNF1A, HNF4A, and GCK cause diabetes through loss of function, while KCNJ11 and ABCC8 can cause either diabetes (gain of function) or hyperinsulinism (loss of function). This lets us check whether performance holds across different genes and mechanisms.

## Approach

1. Build a labeled benchmark from ClinVar: pathogenic/likely pathogenic = 1, benign/likely benign = 0, restricted to missense and canonical splice-site variants (see below).
2. Score each variant with EVEE.
3. Compare EVEE's scores with the ClinVar labels overall, per gene, and separately for missense and splice-site variants.

## Dataset preparation

Variants for GCK, HNF1A, HNF4A, KCNJ11, and ABCC8 were downloaded from ClinVar (download date: `YYYY-MM-DD`) and filtered with `prepare_evee_data.py`:

```bash
python3 prepare_evee_data.py raw/
```

### Filtering criteria

- **Genes:** GCK, HNF1A, HNF4A, KCNJ11, ABCC8. Variants listed under more than one gene are excluded.
- **Labels:** `Target = 1` for Pathogenic / Likely pathogenic, `Target = 0` for Benign / Likely benign. Variants of uncertain significance and conflicting classifications are excluded.
- **Review quality:** at least 1 ClinVar review star (single submitter or better).
- **Variant type:** `missense` or canonical `splice-site` (splice donor/acceptor). Splice-region variants are not treated as splice-site, and records matching both types are excluded as ambiguous.
- **Duplicates:** removed by ClinVar `VariationID`.

### Results

7,374 ClinVar records were loaded, and 21 duplicates were removed. **1,080 variants** passed all filters.

| Variant type | Pathogenic (1) | Benign (0) | Total     |
| ------------ | -------------- | ---------- | --------- |
| Missense     | 791            | 81         | 872       |
| Splice-site  | 206            | 2          | 208       |
| **All**      | **997**        | **83**     | **1,080** |

Breakdown by gene (pathogenic / benign):

| Gene   | Missense | Splice-site | Total |
| ------ | -------- | ----------- | ----- |
| GCK    | 396 / 12 | 64 / 0      | 472   |
| ABCC8  | 138 / 11 | 75 / 2      | 226   |
| HNF1A  | 141 / 33 | 43 / 0      | 217   |
| HNF4A  | 72 / 12  | 24 / 0      | 108   |
| KCNJ11 | 44 / 13  | none        | 57    |

Excluded records: 3,055 had no clear benign/pathogenic classification, 2,817 were neither missense nor canonical splice-site, 323 were multi-gene variants, 72 had insufficient review stars, and 6 were ambiguous in variant type. Each record is counted under the first reason that applied. Reasons for every exclusion are in `evee_excluded.csv`.

### Limitations

- **Class imbalance.** About 92% of included variants are pathogenic (91% of missense, 99% of splice-site). Accuracy alone is misleading, so report AUPRC, balanced accuracy, or sensitivity/specificity alongside AUROC.
- **Splice-site variants have almost no benign examples.** Only 2 of 208 are benign (both ABCC8), and GCK, HNF1A, and HNF4A have none. AUROC and specificity cannot be estimated for this group. The missense vs. splice-site comparison is therefore limited to sensitivity (fraction of pathogenic variants correctly detected) unless benign splice-site controls are added from another source.
- **Few benign variants per gene.** Benign missense counts range from 11 to 33 per gene, so per-gene metrics will have wide confidence intervals. Pool genes and use bootstrap confidence intervals.
- **KCNJ11 has no canonical splice-site variants** in this dataset.
- **Disease mechanism differs by gene.** GCK, HNF1A, and HNF4A cause MODY, while ABCC8 and KCNJ11 cause both neonatal diabetes (gain of function) and hyperinsulinism (loss of function). Pathogenic labels therefore mean different things across genes, so report per-gene results as well as pooled ones.
- **Possible label circularity.** ClinVar classifications can rely on computational predictors or literature that overlap with EVEE's training data. Check whether EVEE was trained or calibrated on ClinVar.
- **Single-star evidence.** Calls from a single submitter are included. Re-run with `--min-stars 2` as a sensitivity analysis.

### Output files

| File                   | Contents                                                                       |
| ---------------------- | ------------------------------------------------------------------------------ |
| `evee_dataset.csv`     | All included variants (`SNP`, `Gene`, `Variant_Type`, `Target`, plus metadata) |
| `evee_missense.csv`    | Missense variants only                                                         |
| `evee_splice_site.csv` | Splice-site variants only                                                      |
| `evee_excluded.csv`    | Excluded records with reasons                                                  |
| `evee_summary.csv`     | Counts by gene, variant type, and target                                       |
