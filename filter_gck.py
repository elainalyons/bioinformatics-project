#!/usr/bin/env python3
"""Filter and annotate a ClinVar GCK download (TSV).

Usage:
    python filter_gck.py raw/clinvar_gck.tsv --outdir data/

Outputs (in --outdir):
    gck_annotated_all.csv  every GCK row, with labels and flags (nothing dropped)
    gck_primary.csv        rows kept for the main analysis
    gck_vus.csv            variants of uncertain significance (for Phase F)
    gck_excluded.csv       labeled rows excluded from the primary set, with reason

The raw file is never modified. Edit the CONFIG section to change rules.
"""
import argparse
import re
from pathlib import Path

import pandas as pd

# ---------------- CONFIG ----------------
MIN_STARS = 1  # minimum review stars for the primary set (0 = keep everything)

PATHOGENIC = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
BENIGN = {"benign", "likely benign", "benign/likely benign"}
VUS = {"uncertain significance"}

# Condition matching (lowercase substrings). Adjust after looking at your data.
MODY_TERMS = ["maturity-onset diabetes of the young", "mody"]
OTHER_DIABETES_TERMS = ["diabetes"]  # checked after MODY, e.g. neonatal diabetes
HYPERINSULINISM_TERMS = ["hyperinsulin", "hypoglycemia"]
NONSPECIFIC_TERMS = ["not provided", "not specified"]
# ----------------------------------------


def review_stars(status):
    """Map ClinVar review status text to the usual 0-4 star scale."""
    s = str(status).lower()
    if "practice guideline" in s:
        return 4
    if "expert panel" in s:
        return 3
    if "multiple submitters, no conflicts" in s:
        return 2
    if "single submitter" in s or "conflicting" in s:
        return 1
    return 0


def label_from_class(cls):
    c = str(cls).strip().lower()
    if c in PATHOGENIC:
        return "pathogenic"
    if c in BENIGN:
        return "benign"
    if c in VUS:
        return "vus"
    if "conflicting" in c:
        return "conflicting"
    return "other"


def condition_flags(cond):
    conds = [c.strip().lower() for c in re.split(r"\|", str(cond)) if c.strip()]
    text = " | ".join(conds)
    has_mody = any(t in text for t in MODY_TERMS)
    has_other_diabetes = (not has_mody) and any(t in text for t in OTHER_DIABETES_TERMS)
    has_hyper = any(t in text for t in HYPERINSULINISM_TERMS)
    nonspecific = bool(conds) and all(any(t in c for t in NONSPECIFIC_TERMS) for c in conds)
    return pd.Series(
        {
            "has_mody": has_mody,
            "has_other_diabetes": has_other_diabetes,
            "has_hyperinsulinism": has_hyper,
            "condition_nonspecific": nonspecific,
        }
    )


def protein_position(row):
    m = re.search(r"(\d+)", str(row.get("Protein change", "")))
    if m:
        return int(m.group(1))
    m = re.search(r"p\.[A-Za-z]{3}(\d+)", str(row.get("Name", "")))
    return int(m.group(1)) if m else pd.NA


def exclusion_reason(r):
    """Return a reason string if a labeled row should be excluded, else ''."""
    if r["is_multigene"]:
        return "multi-gene variant"
    if r["review_stars"] < MIN_STARS:
        return f"review stars < {MIN_STARS}"
    if r["label"] == "pathogenic":
        if r["has_mody"] or r["has_other_diabetes"] or r["condition_nonspecific"]:
            return ""
        if r["has_hyperinsulinism"]:
            return "pathogenic only for hyperinsulinism/hypoglycemia"
        return "pathogenic only for unrelated condition"
    return ""  # benign rows: condition is not required


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv")
    ap.add_argument("--outdir", default="data")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.tsv, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8-sig")

    # --- read TSV and correct extra trailing empty column ---
    import csv

    with open(args.tsv, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")

        header = next(reader)
        rows = []
        malformed = []

        for line_num, row in enumerate(reader, start=2):

            # Remove one extra trailing empty field
            if len(row) == len(header) + 1 and row[-1] == "":
                row = row[:-1]

            if len(row) != len(header):
                malformed.append((line_num, len(row)))
                continue

            rows.append(row)

        if malformed:
            raise ValueError(
                f"Found malformed rows: {malformed[:10]}"
            )

    df = pd.DataFrame(rows, columns=header, dtype=str)

    print(f"Successfully loaded {len(df)} rows")
    print(f"Number of columns: {len(df.columns)}")
    print("\nClassification counts:")
    print(df["Germline classification"].value_counts())

    print("\nReview status counts:")
    print(df["Germline review status"].value_counts())

    df.columns = [c.strip() for c in df.columns]

    # --- keep only rows that involve GCK ---
    # Handles values like "GCK", "GCK:2645", "YKT6:10652|GCK:2645", "GCK, YKT6"
    def parse_genes(g):
        parts = re.split(r"[;,|]", g)
        return {p.split(":")[0].strip().upper() for p in parts if p.strip()}

    genes = df["Gene(s)"].apply(parse_genes)
    df = df[genes.apply(lambda s: "GCK" in s)].copy()
    df["is_multigene"] = genes[df.index].apply(lambda s: len(s) > 1)
    print(f"Rows involving GCK: {len(df)} ({df['is_multigene'].sum()} multi-gene)")

    # --- remove exact duplicate variants ---
    before = len(df)
    df = df.drop_duplicates(subset="VariationID")
    print(f"Duplicate VariationIDs removed: {before - len(df)}")

    # --- labels and annotations ---
    df["label"] = df["Germline classification"].apply(label_from_class)
    df["label_binary"] = df["label"].map({"pathogenic": 1, "benign": 0})
    df["review_stars"] = df["Germline review status"].apply(review_stars)
    df["last_evaluated"] = pd.to_datetime(df["Germline date last evaluated"], errors="coerce")
    df["protein_position"] = df.apply(protein_position, axis=1)
    df = pd.concat([df, df["Condition(s)"].apply(condition_flags)], axis=1)

    print("\nLabel counts (all GCK rows):")
    print(df["label"].value_counts().to_string())

    # --- primary set ---
    labeled = df["label"].isin(["pathogenic", "benign"])
    df["exclusion_reason"] = ""
    df.loc[labeled, "exclusion_reason"] = df[labeled].apply(exclusion_reason, axis=1)
    df["in_primary"] = labeled & (df["exclusion_reason"] == "")

    primary = df[df["in_primary"]]
    excluded = df[labeled & ~df["in_primary"]]
    vus = df[df["label"] == "vus"]

    print("\nPrimary set:")
    print(primary["label"].value_counts().to_string())
    print("\nExclusion reasons (labeled rows):")
    print(excluded["exclusion_reason"].value_counts().to_string())
    print(f"\nVUS rows saved separately: {len(vus)}")

    df.to_csv(out / "gck_annotated_all.csv", index=False)
    primary.to_csv(out / "gck_primary.csv", index=False)
    excluded.to_csv(out / "gck_excluded.csv", index=False)
    vus.to_csv(out / "gck_vus.csv", index=False)
    print(f"\nWrote files to {out.resolve()}")

if __name__ == "__main__":
    main()