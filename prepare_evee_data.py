#!/usr/bin/env python3
"""
Prepare a ClinVar dataset for evaluating EVEE across diabetes-related genes
(GCK, HNF1A, HNF4A, KCNJ11, ABCC8).

Usage:
    python3 prepare_evee_data.py raw/

    python prepare_evee_data.py raw/GCK.tsv raw/HNF1A.tsv raw/HNF4A.tsv \
        raw/KCNJ11.tsv raw/ABCC8.tsv

    # Sensitivity analysis with stricter review requirements:
    python prepare_evee_data.py raw/ --min-stars 2 --outdir data_2star

Outputs (written to --outdir, default: data/):
    evee_dataset.csv          Combined missense + splice-site dataset
    evee_missense.csv         Missense variants only
    evee_splice_site.csv      Splice-site variants only
    evee_excluded.csv         Excluded records and reasons
    evee_summary.csv          Counts by gene, variant type, and target

Target:
    1 = Pathogenic / Likely pathogenic
    0 = Benign / Likely benign

The original TSV files are never modified.
"""

import argparse
import re
from pathlib import Path

import pandas as pd


# ---------------- CONFIG ----------------

TARGET_GENES = {"GCK", "HNF1A", "HNF4A", "KCNJ11", "ABCC8"}

PATHOGENIC = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
}

BENIGN = {
    "benign",
    "likely benign",
    "benign/likely benign",
}

# Minimum ClinVar review stars (override with --min-stars).
# Set to 0 to include all review levels.
MIN_STARS = 1

# Only canonical splice donor/acceptor variants are treated as splice-site.
# Splice-region variants are intentionally NOT counted as splice-site because
# they do not necessarily disrupt normal splicing. Missense variants that also
# carry a splice-region consequence are kept as missense but flagged in the
# `Has_Splice_Region` column so they can be dropped in a sensitivity analysis.

# ----------------------------------------


def normalize(value):
    """Normalize text for consistent matching."""
    return re.sub(r"\s+", " ", str(value).strip().lower())


def parse_genes(value):
    """Extract gene symbols from ClinVar's Gene(s) column."""
    parts = re.split(r"[;,|]", str(value))

    return {
        p.split(":")[0].strip().upper()
        for p in parts
        if p.strip()
    }


def assign_target_gene(genes):
    """
    Deterministic gene assignment.
    Single-gene records get that gene; multi-gene records are labelled MULTI
    (they are excluded later); records with no target gene get NA.
    """
    hits = sorted(genes & TARGET_GENES)

    if not hits:
        return pd.NA

    if len(genes) > 1:
        return "MULTI"

    return hits[0]


def label_from_class(value):
    """Convert ClinVar classification into a binary target."""
    classification = normalize(value)

    if classification in PATHOGENIC:
        return 1

    if classification in BENIGN:
        return 0

    return pd.NA


def review_stars(status):
    """Convert ClinVar review status into a 0-4 star score."""
    s = normalize(status)

    if "practice guideline" in s:
        return 4
    if "expert panel" in s:
        return 3
    if "multiple submitters, no conflicts" in s:
        return 2
    if "single submitter" in s or "conflicting" in s:
        return 1

    return 0


def split_consequences(consequence):
    """Split ClinVar's Molecular consequence column into a set of terms."""
    c = normalize(consequence)
    return {x.strip() for x in re.split(r"[|;,]", c) if x.strip()}


def classify_variant_type(consequence):
    """
    Identify missense and canonical splice-site variants
    using ClinVar's Molecular consequence column.
    """
    consequences = split_consequences(consequence)

    is_missense = any(
        "missense_variant" in x or "missense variant" in x
        for x in consequences
    )

    is_splice = any(
        any(term in x for term in [
            "splice_donor_variant",
            "splice_acceptor_variant",
            "splice donor variant",
            "splice acceptor variant",
            "splice-site variant",
            "splice site variant",
        ])
        for x in consequences
    )

    # Exclude records with ambiguous multiple classifications.
    if is_missense and is_splice:
        return "ambiguous"

    if is_missense:
        return "missense"

    if is_splice:
        return "splice-site"

    return pd.NA


def has_splice_region(consequence):
    """True if the consequence list includes a splice-region term."""
    return any(
        "splice_region" in x or "splice region" in x
        for x in split_consequences(consequence)
    )


def read_clinvar_tsv(filepath):
    """
    Read a ClinVar TSV robustly.

    - Quote characters are NOT interpreted (ClinVar free text can contain
      stray quotes that would otherwise merge columns).
    - Trailing empty fields caused by extra tabs at the end of the header
      and/or data rows are stripped, no matter how many there are.
    - Data rows with extra NON-empty fields raise an error (real misalignment).
    - Rows that are short (missing trailing values) are padded with empty
      strings, and a warning is printed with the count.
    - Blank lines are skipped.
    """
    with open(filepath, "r", encoding="utf-8-sig", newline="") as f:
        lines = f.read().splitlines()

    if not lines:
        raise ValueError(f"{filepath}: file is empty")

    header = lines[0].split("\t")

    # Drop trailing empty header cells (extra trailing tab(s)).
    while header and header[-1].strip() == "":
        header.pop()

    header = [h.strip() for h in header]
    n = len(header)

    rows = []
    padded = 0
    trimmed = 0

    for line_num, line in enumerate(lines[1:], start=2):

        if not line.strip():
            continue

        row = line.split("\t")

        if len(row) > n:
            extra = row[n:]

            if any(cell.strip() for cell in extra):
                raise ValueError(
                    f"{filepath}, line {line_num}: expected {n} columns, "
                    f"got {len(row)} with non-empty extra fields: {extra}"
                )

            row = row[:n]
            trimmed += 1

        elif len(row) < n:
            row = row + [""] * (n - len(row))
            padded += 1

        rows.append(row)

    if trimmed:
        print(f"  {Path(filepath).name}: trimmed trailing empty field(s) "
              f"from {trimmed} rows")

    if padded:
        print(f"  WARNING {Path(filepath).name}: padded {padded} short rows")

    return pd.DataFrame(rows, columns=header, dtype=str)


def find_input_files(inputs):
    """Accept individual TSV files or directories containing TSVs."""
    files = []

    for item in inputs:
        path = Path(item)

        if path.is_dir():
            files.extend(sorted(path.glob("*.tsv")))
        elif path.is_file():
            files.append(path)
        else:
            raise FileNotFoundError(path)

    files = list(dict.fromkeys(files))

    if not files:
        raise FileNotFoundError("No TSV files found.")

    return files


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "inputs",
        nargs="+",
        help="ClinVar TSV files or directories containing TSV files"
    )

    parser.add_argument(
        "--outdir",
        default="data",
        help="Output directory (default: data)"
    )

    parser.add_argument(
        "--min-stars",
        type=int,
        default=MIN_STARS,
        help=f"Minimum ClinVar review stars (default: {MIN_STARS})"
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    min_stars = args.min_stars

    files = find_input_files(args.inputs)

    print(f"Found {len(files)} input files.")

    # --- Load all files ---

    frames = []

    for filepath in files:
        df = read_clinvar_tsv(filepath)
        df["source_file"] = filepath.name
        frames.append(df)

        print(f"Loaded {filepath.name}: {len(df)} rows")

    df = pd.concat(frames, ignore_index=True)

    print(f"\nTotal rows loaded: {len(df)}")

    # --- Validate required columns ---

    required = [
        "Gene(s)",
        "Name",
        "Molecular consequence",
        "Germline classification",
        "Germline review status",
        "VariationID",
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # --- Gene filtering ---

    df["parsed_genes"] = df["Gene(s)"].apply(parse_genes)

    df["target_gene"] = df["parsed_genes"].apply(assign_target_gene)

    df = df[df["target_gene"].notna()].copy()

    df["is_multigene"] = df["parsed_genes"].apply(
        lambda genes: len(genes) > 1
    )

    df["genes_listed"] = df["parsed_genes"].apply(
        lambda genes: ";".join(sorted(genes))
    )

    print(f"Rows involving target genes: {len(df)}")

    # --- Remove duplicate variants ---

    before = len(df)

    df["_variant_key"] = df.apply(
        lambda r: (
            f"VariationID:{r['VariationID']}"
            if str(r["VariationID"]).strip()
            else f"{r['target_gene']}:{r['Name']}"
        ),
        axis=1
    )

    df = df.drop_duplicates(subset="_variant_key").copy()

    print(f"Duplicate variants removed: {before - len(df)}")

    # --- Assign labels ---

    df["Target"] = df["Germline classification"].apply(
        label_from_class
    )

    df["review_stars"] = df["Germline review status"].apply(
        review_stars
    )

    df["Variant_Type"] = df["Molecular consequence"].apply(
        classify_variant_type
    )

    df["Has_Splice_Region"] = df["Molecular consequence"].apply(
        has_splice_region
    )

    # --- Determine exclusions ---

    def exclusion_reason(row):

        if row["is_multigene"]:
            return "Multi-gene variant"

        if pd.isna(row["Target"]):
            return "Not a clear benign/pathogenic classification"

        if row["review_stars"] < min_stars:
            return "Insufficient ClinVar review stars"

        if pd.isna(row["Variant_Type"]):
            return "Not missense or canonical splice-site"

        if row["Variant_Type"] == "ambiguous":
            return "Ambiguous variant type"

        return ""

    df["exclusion_reason"] = df.apply(
        exclusion_reason,
        axis=1
    )

    primary = df[df["exclusion_reason"] == ""].copy()
    excluded = df[df["exclusion_reason"] != ""].copy()

    # Convert target to integer.
    primary["Target"] = primary["Target"].astype(int)

    # --- Create final dataset ---

    # Name is retained as the ClinVar variant identifier.
    # Protein change and genomic coordinates are included for EVEE input.
    output_columns = [
        "Name",
        "target_gene",
        "Variant_Type",
        "Target",
        "Has_Splice_Region",
        "Protein change",
        "VariationID",
        "GRCh38Chromosome",
        "GRCh38Location",
        "Canonical SPDI",
        "dbSNP ID",
        "Condition(s)",
        "Molecular consequence",
        "Germline classification",
        "Germline review status",
        "review_stars",
        "Germline date last evaluated",
        "source_file",
    ]

    # Only retain columns that exist in the downloaded files.
    output_columns = [
        c for c in output_columns if c in primary.columns
    ]

    missing_optional = [
        c for c in ["GRCh38Chromosome", "GRCh38Location",
                    "Canonical SPDI", "dbSNP ID", "Protein change"]
        if c not in df.columns
    ]

    if missing_optional:
        print(f"\nNote: optional columns not in downloads: {missing_optional}")

    primary = primary[output_columns].copy()

    primary = primary.rename(columns={
        "Name": "SNP",
        "target_gene": "Gene",
        "Protein change": "Protein_Change",
        "Condition(s)": "Conditions",
    })

    # Excluded records: keep useful columns for auditing.
    excluded_columns = [
        c for c in [
            "Name", "target_gene", "genes_listed", "Variant_Type",
            "Target", "Has_Splice_Region", "VariationID",
            "Molecular consequence", "Germline classification",
            "Germline review status", "review_stars",
            "exclusion_reason", "source_file",
        ] if c in excluded.columns
    ]

    excluded = excluded[excluded_columns].copy()

    # --- Split by variant type ---

    missense = primary[
        primary["Variant_Type"] == "missense"
    ].copy()

    splice = primary[
        primary["Variant_Type"] == "splice-site"
    ].copy()

    # --- Summary ---

    summary = (
        primary.groupby(
            ["Gene", "Variant_Type", "Target"]
        )
        .size()
        .reset_index(name="Count")
    )

    # --- Save outputs ---

    primary.to_csv(outdir / "evee_dataset.csv", index=False)
    missense.to_csv(outdir / "evee_missense.csv", index=False)
    splice.to_csv(outdir / "evee_splice_site.csv", index=False)
    excluded.to_csv(outdir / "evee_excluded.csv", index=False)
    summary.to_csv(outdir / "evee_summary.csv", index=False)

    # --- Print results ---

    print(f"\nFINAL DATASET (min review stars = {min_stars})")
    print(f"Total included: {len(primary)}")

    print("\nTarget counts:")
    print(
        primary["Target"]
        .value_counts()
        .rename({1: "Pathogenic", 0: "Benign"})
    )

    print("\nVariant type counts:")
    print(primary["Variant_Type"].value_counts())

    print("\nGene counts:")
    print(primary["Gene"].value_counts())

    print("\nGene / Variant Type / Target breakdown:")
    print(summary.to_string(index=False))

    n_sr = int(missense["Has_Splice_Region"].sum())
    print(f"\nMissense variants also flagged splice-region: {n_sr}")

    # Warn about sparse classes that make AUROC/AUPRC unreliable.
    by_type = primary.groupby(["Variant_Type", "Target"]).size().unstack(fill_value=0)

    for vtype, counts in by_type.iterrows():
        for target, label in [(0, "benign"), (1, "pathogenic")]:
            if counts.get(target, 0) < 10:
                print(f"WARNING: only {counts.get(target, 0)} {label} "
                      f"{vtype} variants - metrics will be unreliable.")

    print("\nExcluded records:")
    print(excluded["exclusion_reason"].value_counts())

    print(f"\nResults saved to: {outdir.resolve()}")


if __name__ == "__main__":
    main()