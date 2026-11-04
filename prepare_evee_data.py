#!/usr/bin/env python3
"""
Prepare a ClinVar dataset for evaluating EVEE across diabetes-related genes
(GCK, HNF1A, HNF4A, KCNJ11, ABCC8).

Usage:
    python3 prepare_evee_data.py raw/
    
    python prepare_evee_dataset.py raw/GCK.tsv raw/HNF1A.tsv raw/HNF4A.tsv \
        raw/KCNJ11.tsv raw/ABCC8.tsv

    # Sensitivity analysis with stricter review requirements:
    python prepare_evee_dataset.py raw/ --min-stars 2 --outdir data_2star

Outputs (written to --outdir, default: data/):

  Benchmark (labeled, single-gene, >= --min-stars):
    evee_dataset.csv          Missense + splice-site benchmark
    evee_missense.csv         Missense variants only
    evee_splice_site.csv      Splice-site variants only

  For the wider analyses (VUS, review-star / domain / conservation strata):
    evee_all_variants.csv     ALL single-gene missense + splice-site variants,
                              every ClinVar class and review level, with
                              Class_Group and In_Benchmark flags
    evee_vus.csv              Variants of uncertain significance only

  Audit and summary:
    evee_excluded.csv         Records left out of the benchmark, with reasons
    evee_summary.csv          Benchmark counts by gene, variant type, and target
    evee_overview.md          Readable results overview (also printed)
    evee_overview.png         Overview figure (only if matplotlib is installed)

Target:
    1 = Pathogenic / Likely pathogenic
    0 = Benign / Likely benign
    (blank in evee_all_variants.csv / evee_vus.csv for VUS, conflicting, etc.)

The original TSV files are never modified.
"""

import argparse
import re
from pathlib import Path

import pandas as pd


# ---------------- CONFIG ----------------

TARGET_GENES = {"GCK", "HNF1A", "HNF4A", "KCNJ11", "ABCC8"}

# Display order for tables and figures.
GENE_ORDER = ["GCK", "HNF1A", "HNF4A", "KCNJ11", "ABCC8"]
TYPE_ORDER = ["missense", "splice-site"]

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

VUS = {"uncertain significance"}

# Minimum ClinVar review stars for the benchmark (override with --min-stars).
# Set to 0 to include all review levels.
MIN_STARS = 1

# Only canonical splice donor/acceptor variants are treated as splice-site.
# Splice-region variants are intentionally NOT counted as splice-site because
# they do not necessarily disrupt normal splicing. Missense variants that also
# carry a splice-region consequence are kept as missense but flagged in the
# `Has_Splice_Region` column so they can be dropped in a sensitivity analysis.

# Sample-size thresholds (per class) used to judge whether AUROC/AUPRC are
# meaningful for a group.
MIN_CLASS_OK = 30
MIN_CLASS_LIMITED = 10

EVIDENCE_LABELS = {
    0: "0 stars - no criteria / assertion only",
    1: "1 star - single submitter or conflicting",
    2: "2 stars - multiple submitters, no conflicts",
    3: "3 stars - expert panel",
    4: "4 stars - practice guideline",
}

REASON_MULTI = "Multi-gene variant"
REASON_NOT_CLEAR = "Not a clear benign/pathogenic classification"
REASON_STARS = "Insufficient ClinVar review stars"
REASON_TYPE = "Not missense or canonical splice-site"
REASON_AMBIGUOUS = "Ambiguous variant type"
REASON_ORDER = [REASON_MULTI, REASON_NOT_CLEAR, REASON_STARS,
                REASON_TYPE, REASON_AMBIGUOUS]

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

TRANSCRIPT_RE = re.compile(r"^(N[MR]_\d+(?:\.\d+)?)")
HGVS_C_RE = re.compile(r"\):(c\.[^\s]+)")
HGVS_P_RE = re.compile(r"\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*|=|\?)\)")
PROTEIN_CHANGE_RE = re.compile(r"^([A-Z])(\d+)([A-Z*])")

# Columns written to the output CSVs (only those present are kept).
OUTPUT_COLUMNS = [
    "Name",
    "target_gene",
    "Variant_Type",
    "Class_Group",
    "Target",
    "Has_Splice_Region",
    "Transcript",
    "HGVS_c",
    "Protein change",
    "Ref_AA",
    "Protein_Position",
    "Alt_AA",
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

RENAMES = {
    "Name": "SNP",
    "target_gene": "Gene",
    "Protein change": "Protein_Change",
    "Condition(s)": "Conditions",
}

OPTIONAL_COLUMNS = ["GRCh38Chromosome", "GRCh38Location",
                    "Canonical SPDI", "dbSNP ID", "Protein change"]

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


def class_group(value):
    """Group ClinVar's germline classification into broad classes."""
    c = normalize(value)

    if c in PATHOGENIC:
        return "Pathogenic"

    if c in BENIGN:
        return "Benign"

    if c in VUS:
        return "VUS"

    if c.startswith("conflicting"):
        return "Conflicting"

    return "Other"


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


def parse_transcript(name):
    m = TRANSCRIPT_RE.match(str(name))
    return m.group(1) if m else pd.NA


def parse_hgvs_c(name):
    m = HGVS_C_RE.search(str(name))
    return m.group(1) if m else pd.NA


def parse_protein(name, protein_change):
    """
    Return (ref_aa, position, alt_aa) for a protein-level change.
    Uses the HGVS p. notation in Name first, then the Protein change column.
    Position is NA for variants without a protein-level description
    (for example canonical splice-site variants).
    """
    m = HGVS_P_RE.search(str(name))

    if m:
        ref = THREE_TO_ONE.get(m.group(1).capitalize(), pd.NA)
        pos = int(m.group(2))
        alt = m.group(3)
        if len(alt) == 3:
            alt = THREE_TO_ONE.get(alt.capitalize(), pd.NA)
        return ref, pos, alt

    first = str(protein_change).split(",")[0].strip()
    m = PROTEIN_CHANGE_RE.match(first)

    if m:
        return m.group(1), int(m.group(2)), m.group(3)

    return pd.NA, pd.NA, pd.NA


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


def select_output(frame, extra=()):
    """Pick and rename the output columns that exist in the frame."""
    cols = [c for c in OUTPUT_COLUMNS + list(extra) if c in frame.columns]
    return frame[cols].rename(columns=RENAMES)


# ---------------- OVERVIEW HELPERS ----------------


def fmt(n):
    return f"{int(n):,}"


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "n/a"


def md_table(frame):
    """Render a DataFrame as a markdown table (no extra dependencies)."""
    cols = [str(c) for c in frame.columns]
    lines = [
        "| " + " | ".join(cols) + " |",
        "|" + "|".join("---" for _ in cols) + "|",
    ]

    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(v) for v in row) + " |")

    return "\n".join(lines)


def feasibility(pathogenic, benign):
    """Is there enough of each class for AUROC / AUPRC to be meaningful?"""
    smallest = min(pathogenic, benign)

    if smallest >= MIN_CLASS_OK:
        return f"Yes (>= {MIN_CLASS_OK} per class)"

    if smallest >= MIN_CLASS_LIMITED:
        return "Limited (wide CIs)"

    return f"No (< {MIN_CLASS_LIMITED} in a class)"


def count_table(frame, value_col, classes):
    """
    Count variants by (variant type, gene) and value_col, with pooled rows
    ("All genes" per variant type, and one overall row).
    """
    idx = pd.MultiIndex.from_product(
        [TYPE_ORDER, GENE_ORDER], names=["Variant_Type", "Gene"]
    )

    counts = (
        frame.groupby(["Variant_Type", "Gene", value_col])
        .size()
        .unstack(value_col, fill_value=0)
        .reindex(idx, fill_value=0)
    )

    for c in classes:
        if c not in counts.columns:
            counts[c] = 0

    counts = counts[classes]

    rows = []

    for vt in TYPE_ORDER:
        for gene in GENE_ORDER:
            vals = [int(counts.loc[(vt, gene), c]) for c in classes]
            if sum(vals) == 0:
                continue
            rows.append([vt, gene] + vals)

        sub = counts.loc[vt]
        rows.append([vt, "All genes"] + [int(sub[c].sum()) for c in classes])

    rows.append(["All", "All genes"] + [int(counts[c].sum()) for c in classes])

    out = pd.DataFrame(rows, columns=["Variant type", "Gene"] + classes)
    out["Total"] = out[classes].sum(axis=1)

    return out


def benchmark_table(bench):
    """Benchmark counts with % pathogenic and a metric-feasibility column."""
    bench = bench.copy()
    bench["Label"] = bench["Target"].map({1: "Pathogenic", 0: "Benign"})

    out = count_table(bench, "Label", ["Pathogenic", "Benign"])
    out["% Pathogenic"] = [
        pct(p, t) for p, t in zip(out["Pathogenic"], out["Total"])
    ]
    out["AUROC / AUPRC meaningful?"] = [
        feasibility(p, b) for p, b in zip(out["Pathogenic"], out["Benign"])
    ]

    return out


def stars_table(all_variants, min_stars):
    """Labeled variants by ClinVar review stars (all levels)."""
    labeled = all_variants[
        all_variants["Class_Group"].isin(["Pathogenic", "Benign"])
    ]

    counts = (
        labeled.groupby(["review_stars", "Class_Group"])
        .size()
        .unstack("Class_Group", fill_value=0)
        .reindex(range(5), fill_value=0)
    )

    for c in ["Pathogenic", "Benign"]:
        if c not in counts.columns:
            counts[c] = 0

    out = pd.DataFrame({
        "Stars": list(range(5)),
        "Evidence tier": [EVIDENCE_LABELS[s] for s in range(5)],
        "Pathogenic": counts["Pathogenic"].astype(int).values,
        "Benign": counts["Benign"].astype(int).values,
    })

    out["Total"] = out["Pathogenic"] + out["Benign"]
    out["In benchmark?"] = [
        "Yes" if s >= min_stars else "No (below --min-stars)"
        for s in out["Stars"]
    ]

    return out


def make_figure(bench_tbl, stars_tbl, path):
    """Two-panel overview figure. Returns False if matplotlib is missing."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    path_color = "#D55E00"
    ben_color = "#0072B2"

    rows = bench_tbl[
        (bench_tbl["Gene"] != "All genes") & (bench_tbl["Total"] > 0)
    ]

    labels = [f"{g} {vt}" for g, vt in zip(rows["Gene"], rows["Variant type"])]
    path_n = rows["Pathogenic"].tolist()
    ben_n = rows["Benign"].tolist()
    ypos = list(range(len(rows)))[::-1]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.5, 4.8), gridspec_kw={"width_ratios": [1.5, 1]}
    )

    ax1.barh(ypos, path_n, color=path_color, label="Pathogenic")
    ax1.barh(ypos, ben_n, left=path_n, color=ben_color, label="Benign")

    top = max((p + b for p, b in zip(path_n, ben_n)), default=1)

    for y, p, b in zip(ypos, path_n, ben_n):
        ax1.text(p + b + top * 0.01, y, f"{p} / {b}", va="center", fontsize=8)

    ax1.set_yticks(ypos)
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlim(0, top * 1.15)
    ax1.set_xlabel("Variants (labels show pathogenic / benign)")
    ax1.set_title("Benchmark by gene and variant type", fontsize=11)
    ax1.legend(frameon=False, loc="lower right")

    x = stars_tbl["Stars"].tolist()
    in_bench = [s.startswith("Yes") for s in stars_tbl["In benchmark?"]]
    alphas = [1.0 if k else 0.35 for k in in_bench]
    sp = stars_tbl["Pathogenic"].tolist()
    sb = stars_tbl["Benign"].tolist()

    for xi, p, b, a in zip(x, sp, sb, alphas):
        ax2.bar(xi, p, color=path_color, alpha=a)
        ax2.bar(xi, b, bottom=p, color=ben_color, alpha=a)
        ax2.text(xi, p + b, f"{p} / {b}", ha="center", va="bottom", fontsize=8)

    ax2.set_xticks(x)
    ax2.set_xlabel("ClinVar review stars (faded = below --min-stars)")
    ax2.set_ylabel("Labeled variants")
    ax2.set_title("Labeled variants by evidence level", fontsize=11)
    ax2.margins(y=0.15)

    for ax in (ax1, ax2):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return True


def build_overview(funnel, bench_tbl, class_tbl, stars_tbl, notes,
                   min_stars, figure_name):
    """Assemble the markdown results overview."""
    pooled = bench_tbl[bench_tbl["Gene"] == "All genes"].set_index("Variant type")

    def headline(vt, label):
        if vt not in pooled.index:
            return f"- {label}: none"
        r = pooled.loc[vt]
        return (f"- **{label}:** {fmt(r['Total'])} variants "
                f"({fmt(r['Pathogenic'])} pathogenic, {fmt(r['Benign'])} benign; "
                f"{r['% Pathogenic']} pathogenic)")

    lines = [
        "# EVEE benchmark overview",
        "",
        f"Genes: {', '.join(GENE_ORDER)}  |  "
        f"minimum review stars: {min_stars}",
        "",
        "## Headline",
        "",
        headline("All", "Benchmark"),
        headline("missense", "Missense"),
        headline("splice-site", "Splice-site"),
        "",
        "## 1. Where the records went",
        "",
        md_table(funnel),
        "",
        "## 2. Benchmark composition (labeled variants used for evaluation)",
        "",
        md_table(bench_tbl),
        "",
        f"Thresholds: 'Yes' needs at least {MIN_CLASS_OK} variants in the "
        f"smaller class; 'Limited' needs at least {MIN_CLASS_LIMITED}.",
        "",
        "## 3. Everything available for analysis (all review levels)",
        "",
        "All single-gene missense and canonical splice-site variants, "
        "including uncertain (VUS) and conflicting variants "
        "(`evee_all_variants.csv`; VUS only in `evee_vus.csv`).",
        "",
        md_table(class_tbl),
        "",
        "## 4. Evidence level of labeled variants",
        "",
        md_table(stars_tbl.drop(columns="Stars")),
        "",
        "## 5. Notes",
        "",
    ]

    lines += [f"- {n}" for n in notes]

    if figure_name:
        lines += ["", f"![Benchmark overview]({figure_name})"]

    return "\n".join(lines) + "\n"


# ---------------- MAIN ----------------


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

    n_loaded = len(df)

    print(f"\nTotal rows loaded: {n_loaded}")

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

    n_target = len(df)

    if df.empty:
        raise SystemExit("No records involve the target genes.")

    df["is_multigene"] = df["parsed_genes"].apply(
        lambda genes: len(genes) > 1
    )

    df["genes_listed"] = df["parsed_genes"].apply(
        lambda genes: ";".join(sorted(genes))
    )

    print(f"Rows involving target genes: {n_target}")

    # --- Remove duplicate variants ---

    df["_variant_key"] = df.apply(
        lambda r: (
            f"VariationID:{r['VariationID']}"
            if str(r["VariationID"]).strip()
            else f"{r['target_gene']}:{r['Name']}"
        ),
        axis=1
    )

    df = df.drop_duplicates(subset="_variant_key").copy()

    n_dups = n_target - len(df)
    n_unique = len(df)

    print(f"Duplicate variants removed: {n_dups}")

    # --- Assign labels and annotations ---

    df["Class_Group"] = df["Germline classification"].apply(class_group)

    df["Target"] = (
        df["Class_Group"]
        .map({"Pathogenic": 1, "Benign": 0})
        .astype("Int64")
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

    protein_col = (
        df["Protein change"] if "Protein change" in df.columns
        else pd.Series([""] * len(df), index=df.index)
    )

    parsed = [
        parse_protein(n, p) for n, p in zip(df["Name"], protein_col)
    ]

    df["Ref_AA"] = [p[0] for p in parsed]
    df["Protein_Position"] = pd.array([p[1] for p in parsed], dtype="Int64")
    df["Alt_AA"] = [p[2] for p in parsed]
    df["Transcript"] = df["Name"].apply(parse_transcript)
    df["HGVS_c"] = df["Name"].apply(parse_hgvs_c)

    # --- Determine exclusions from the benchmark ---

    def exclusion_reason(row):

        if row["is_multigene"]:
            return REASON_MULTI

        if row["Class_Group"] not in ("Pathogenic", "Benign"):
            return REASON_NOT_CLEAR

        if row["review_stars"] < min_stars:
            return REASON_STARS

        if pd.isna(row["Variant_Type"]):
            return REASON_TYPE

        if row["Variant_Type"] == "ambiguous":
            return REASON_AMBIGUOUS

        return ""

    df["exclusion_reason"] = df.apply(
        exclusion_reason,
        axis=1
    )

    df["In_Benchmark"] = df["exclusion_reason"] == ""

    # --- Benchmark (labeled, single-gene, >= min stars) ---

    primary = df[df["In_Benchmark"]].copy()
    primary["Target"] = primary["Target"].astype(int)

    # --- All single-gene missense / splice-site variants, any class ---

    analysis = df[
        (~df["is_multigene"])
        & df["Variant_Type"].isin(TYPE_ORDER)
    ].copy()

    vus = analysis[analysis["Class_Group"] == "VUS"].copy()

    # --- Excluded records (audit trail) ---

    excluded = df[~df["In_Benchmark"]].copy()

    excluded_columns = [
        c for c in [
            "Name", "target_gene", "genes_listed", "Variant_Type",
            "Class_Group", "Target", "Has_Splice_Region", "VariationID",
            "Molecular consequence", "Germline classification",
            "Germline review status", "review_stars",
            "exclusion_reason", "source_file",
        ] if c in excluded.columns
    ]

    excluded = excluded[excluded_columns].copy()

    # --- Final output tables ---

    primary_out = select_output(primary)
    analysis_out = select_output(analysis, extra=["In_Benchmark"])
    vus_out = select_output(vus, extra=["In_Benchmark"])

    missense = primary_out[primary_out["Variant_Type"] == "missense"].copy()
    splice = primary_out[primary_out["Variant_Type"] == "splice-site"].copy()

    if primary_out.empty:
        raise SystemExit("No variants passed the benchmark filters.")

    summary = (
        primary_out.groupby(["Gene", "Variant_Type", "Target"])
        .size()
        .reset_index(name="Count")
    )

    # --- Save data outputs ---

    primary_out.to_csv(outdir / "evee_dataset.csv", index=False)
    missense.to_csv(outdir / "evee_missense.csv", index=False)
    splice.to_csv(outdir / "evee_splice_site.csv", index=False)
    analysis_out.to_csv(outdir / "evee_all_variants.csv", index=False)
    vus_out.to_csv(outdir / "evee_vus.csv", index=False)
    excluded.to_csv(outdir / "evee_excluded.csv", index=False)
    summary.to_csv(outdir / "evee_summary.csv", index=False)

    # --- Build the overview ---

    # Funnel: where every unique record ended up.
    reason_counts = excluded["exclusion_reason"].value_counts()

    not_clear = excluded[
        excluded["exclusion_reason"] == REASON_NOT_CLEAR
    ]["Class_Group"].value_counts()

    funnel_rows = [
        ["Records loaded (all files)", fmt(n_loaded)],
        ["Not in target genes", fmt(n_loaded - n_target)],
        ["Duplicate records removed", fmt(n_dups)],
        ["**Unique variants**", f"**{fmt(n_unique)}**"],
    ]

    for reason in REASON_ORDER:
        funnel_rows.append(
            [f"Excluded: {reason}", fmt(reason_counts.get(reason, 0))]
        )

        if reason == REASON_NOT_CLEAR:
            for grp in ["VUS", "Conflicting", "Other"]:
                funnel_rows.append(
                    [f"↳ of which {grp}", fmt(not_clear.get(grp, 0))]
                )

    funnel_rows.append(
        ["**Included in benchmark**", f"**{fmt(len(primary_out))}**"]
    )

    funnel = pd.DataFrame(funnel_rows, columns=["Step", "Records"])

    bench_tbl = benchmark_table(primary_out)

    class_tbl = count_table(
        analysis_out, "Class_Group",
        ["Pathogenic", "Benign", "VUS", "Conflicting", "Other"],
    )

    stars_tbl = stars_table(analysis_out, min_stars)

    # Notes: automatic data-quality observations.
    notes = []

    per_type = bench_tbl[bench_tbl["Gene"] == "All genes"].set_index("Variant type")

    for vt in TYPE_ORDER:
        if vt in per_type.index:
            r = per_type.loc[vt]
            if r["Benign"] < MIN_CLASS_LIMITED:
                notes.append(
                    f"Only {int(r['Benign'])} benign {vt} variants in the "
                    f"benchmark: discrimination metrics (AUROC, AUPRC, "
                    f"specificity) are not meaningful for {vt}. Report "
                    f"sensitivity on pathogenic variants only, or add benign "
                    f"controls from another source."
                )

    weak_missense = []

    for vt in TYPE_ORDER:
        for gene in GENE_ORDER:
            sub = bench_tbl[
                (bench_tbl["Variant type"] == vt) & (bench_tbl["Gene"] == gene)
            ]
            if sub.empty:
                notes.append(f"No {vt} variants for {gene} in the benchmark.")
            elif (vt == "missense" and
                  min(sub["Pathogenic"].iloc[0], sub["Benign"].iloc[0])
                  < MIN_CLASS_OK):
                weak_missense.append(gene)

    if weak_missense:
        notes.append(
            f"Fewer than {MIN_CLASS_OK} benign or pathogenic missense variants "
            f"for {', '.join(weak_missense)}: per-gene metrics will have wide "
            f"confidence intervals, so pool genes and bootstrap."
        )

    overall = per_type.loc["All"] if "All" in per_type.index else None

    if overall is not None and overall["Total"]:
        share = overall["Pathogenic"] / overall["Total"]
        if share > 0.8:
            notes.append(
                f"The benchmark is {100 * share:.0f}% pathogenic overall, so "
                f"accuracy alone is misleading (predicting 'pathogenic' for "
                f"everything would score {100 * share:.0f}%). Use AUPRC, MCC, "
                f"balanced accuracy, or sensitivity/specificity with "
                f"bootstrap confidence intervals."
            )

    n_vus = int((analysis_out["Class_Group"] == "VUS").sum())

    notes.append(
        f"{fmt(n_vus)} single-gene missense/splice-site variants are VUS "
        f"(`evee_vus.csv`). They have no ground truth, but their EVEE score "
        f"distribution can be compared with the pathogenic and benign groups."
    )

    n_miss = int((primary_out["Variant_Type"] == "missense").sum())
    n_pos = int(primary_out[
        primary_out["Variant_Type"] == "missense"
    ]["Protein_Position"].notna().sum())

    notes.append(
        f"Protein position was parsed for {fmt(n_pos)} of {fmt(n_miss)} "
        f"benchmark missense variants (needed to join protein-domain and "
        f"conservation annotations)."
    )

    n_sr = int(primary_out[
        primary_out["Variant_Type"] == "missense"
    ]["Has_Splice_Region"].sum())

    notes.append(
        f"{fmt(n_sr)} benchmark missense variants also carry a splice-region "
        f"consequence (`Has_Splice_Region`); drop them for a sensitivity "
        f"analysis if the number is non-trivial."
    )

    missing_optional = [c for c in OPTIONAL_COLUMNS if c not in df.columns]

    if missing_optional:
        notes.append(
            f"Optional columns not present in the downloads: "
            f"{', '.join(missing_optional)}."
        )

    notes.append(
        "Pathogenic labels mean different things across genes (GCK, HNF1A, "
        "HNF4A: MODY; KCNJ11, ABCC8: neonatal diabetes and hyperinsulinism), "
        "so report per-gene results alongside pooled ones."
    )

    figure_name = "evee_overview.png"
    figure_made = make_figure(bench_tbl, stars_tbl, outdir / figure_name)

    if not figure_made:
        print("\nNote: matplotlib not installed - skipped overview figure.")

    overview = build_overview(
        funnel, bench_tbl, class_tbl, stars_tbl, notes, min_stars,
        figure_name if figure_made else None,
    )

    (outdir / "evee_overview.md").write_text(overview, encoding="utf-8")

    # --- Print results ---

    print("\n" + "=" * 70 + "\n")
    print(overview)
    print("=" * 70)
    print(f"\nResults saved to: {outdir.resolve()}")


if __name__ == "__main__":
    main()