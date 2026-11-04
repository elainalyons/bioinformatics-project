# GCK ClinVar Variant Analysis

## Overview

This project analyzes genetic variants associated with the **GCK gene** using data from [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).

The goal is to identify and organize GCK variants based on their clinical classifications and reported associated conditions, with a particular focus on diabetes-related variants.

## Methods

A Python script (`filter_gck.py`) processes the downloaded ClinVar dataset and:

1. Identifies records involving the GCK gene.
2. Removes duplicate variants using their Variation IDs.
3. Groups variants by clinical classification (pathogenic, benign, uncertain significance, conflicting, etc.).
4. Evaluates ClinVar review status to account for the level of supporting evidence.
5. Examines associated conditions to identify diabetes- and MODY-related variants.
6. Separates relevant variants from records involving multiple genes, unrelated conditions, or insufficient review evidence.

The script exports the processed records into CSV files for further analysis. Excluded records and variants of uncertain significance are also saved separately.

## Results

The original dataset contained **1,454 ClinVar records involving GCK**.

Initial classification breakdown:

| Classification              | Number of records |
| --------------------------- | ----------------: |
| Pathogenic                  |               369 |
| Likely pathogenic           |               302 |
| Uncertain significance      |               357 |
| Conflicting classifications |               144 |
| Likely benign               |               107 |
| Benign                      |                40 |
| Other                       |                27 |

After filtering, the primary dataset contained **855 records**:

| Category                       | Number of records |
| ------------------------------ | ----------------: |
| Pathogenic / likely pathogenic |               666 |
| Benign / likely benign         |               189 |
| **Total**                      |           **855** |

An additional 357 records classified as variants of uncertain significance (VUS) were saved separately.

Records excluded from the primary dataset included variants affecting multiple genes, variants associated with conditions unrelated to the project's focus, and variants with insufficient review evidence.

## Data Source

- [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
