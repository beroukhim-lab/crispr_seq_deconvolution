#!/usr/bin/env python3
"""
Traverse CRISPResso output folders, aggregate deletion statistics from
Alleles_frequency_table*.txt files, and emit a long-format summary.
"""

import argparse
from pathlib import Path
from typing import List, Tuple

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize deletion counts from CRISPResso Alleles_frequency_table files."
    )
    parser.add_argument(
        "-i",
        "--input_root",
        required=True,
        help="Root directory containing CRISPResso output folders.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="allele_deletion_summary.tsv",
        help="Path to write the long-format summary TSV (default: allele_deletion_summary.tsv).",
    )
    parser.add_argument(
        "-p",
        "--plot",
        help="Path to write stacked bar plot (default: summary output stem + '_stacked_bar.pdf').",
    )
    parser.add_argument(
        "--pie",
        help="Path to write per-sample pie charts PDF (default: summary output stem + '_pies.pdf').",
    )
    return parser.parse_args()


def infer_sgrna_from_filename(fname: str) -> str:
    """Extract an sgRNA label from the table filename."""
    stem = Path(fname).stem  # e.g., Alleles_frequency_table_around_sgRNA
    prefix = "Alleles_frequency_table"
    if stem.startswith(prefix):
        label = stem[len(prefix):]
        label = label.lstrip("_")
        return label if label else "NA"
    return "NA"


def collect_tables(input_root: Path) -> List[Path]:
    return list(input_root.rglob("Alleles_frequency_table*.txt"))


def summarize_file(path: Path) -> Tuple[str, str, int, int, int]:
    """
    Return (sample, sgrna, no_del_reads, in_frame_del_reads, out_frame_del_reads)
    for one Alleles_frequency_table file.
    """
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as exc:  # pragma: no cover - defensive
        raise RuntimeError(f"Failed to read {path}: {exc}")

    required_cols = {"n_deleted", "#Reads"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {', '.join(sorted(missing))}")

    df["n_deleted"] = pd.to_numeric(df["n_deleted"], errors="coerce").fillna(0).astype(int)
    df["#Reads"] = pd.to_numeric(df["#Reads"], errors="coerce").fillna(0).astype(int)

    no_del_reads = int(df.loc[df["n_deleted"] == 0, "#Reads"].sum())
    in_frame_reads = int(df.loc[(df["n_deleted"] > 0) & (df["n_deleted"] % 3 == 0), "#Reads"].sum())
    out_frame_reads = int(df.loc[(df["n_deleted"] > 0) & (df["n_deleted"] % 3 != 0), "#Reads"].sum())

    sample_raw = path.parent.name
    sample = sample_raw
    for prefix in ["CRISPResso_on_crisprseq_", "CRISPResso_on_", "crisprseq_"]:
        if sample.startswith(prefix):
            sample = sample[len(prefix):]
            break
    sgrna = infer_sgrna_from_filename(path.name)
    return sample, sgrna, no_del_reads, in_frame_reads, out_frame_reads


def add_percentages(df: pd.DataFrame) -> pd.DataFrame:
    totals = (
        df["number of reads with no deletion"]
        + df["number of reads with in-frame deletion"]
        + df["number of reads with out-of-frame deletion"]
    )
    df["percent of reads with no deletion"] = (df["number of reads with no deletion"] / totals.replace(0, pd.NA) * 100).fillna(0)
    df["percent of reads with in-frame deletion"] = (df["number of reads with in-frame deletion"] / totals.replace(0, pd.NA) * 100).fillna(0)
    df["percent of reads with out-of-frame deletion"] = (df["number of reads with out-of-frame deletion"] / totals.replace(0, pd.NA) * 100).fillna(0)
    return df


def plot_stacked_bars(df: pd.DataFrame, plot_path: Path) -> None:
    plot_df = df.copy()
    plot_df["label"] = plot_df.apply(
        lambda row: f"{row['sample']} ({row['sgrna']})" if row["sgrna"] != "NA" else row["sample"], axis=1
    )
    melt = plot_df.melt(
        id_vars=["label"],
        value_vars=[
            "percent of reads with no deletion",
            "percent of reads with in-frame deletion",
            "percent of reads with out-of-frame deletion",
        ],
        var_name="class",
        value_name="percent",
    )
    melt["fraction"] = melt["percent"] / 100.0

    order = sorted(plot_df["label"].tolist())
    class_order = [
        "percent of reads with no deletion",
        "percent of reads with in-frame deletion",
        "percent of reads with out-of-frame deletion",
    ]
    colors = {
        "percent of reads with no deletion": "grey",
        "percent of reads with in-frame deletion": "orange",
        "percent of reads with out-of-frame deletion": "firebrick",
    }

    plt.figure(figsize=(10, max(6, len(order) * 0.5)))
    left = pd.Series(0, index=order, dtype=float)

    for cls in class_order:
        subset = melt[melt["class"] == cls].set_index("label").reindex(order)
        heights = subset["fraction"].fillna(0)
        plt.barh(order, heights, left=left, color=colors[cls], label=cls)
        left = left + heights

    plt.xlabel("Fraction of reads")
    plt.ylabel("Sample")
    plt.xlim(0, 1)
    plt.legend(title="Class", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()


def plot_pies(df: pd.DataFrame, pie_path: Path) -> None:
    colors = ["grey", "orange", "firebrick"]
    with PdfPages(pie_path) as pdf:
        for _, row in df.iterrows():
            vals = [
                row["number of reads with no deletion"],
                row["number of reads with in-frame deletion"],
                row["number of reads with out-of-frame deletion"],
            ]
            total = sum(vals)
            if total <= 0:
                vals = [1, 0, 0]  # avoid zero-size pie; interpret as all unedited
            fig, ax = plt.subplots(figsize=(5, 5))
            wedges, _ = ax.pie(
                vals,
                labels=None,  # use legend only
                colors=colors,
                autopct=None,
                startangle=90,
            )
            fig.subplots_adjust(left=0.05, right=0.7)
            legend_labels = []
            for label, count in zip(
                ["no deletion", "in-frame", "out-of-frame"], vals
            ):
                pct = 0 if total == 0 else (count / total * 100)
                legend_labels.append(f"{label} ({count}, {pct:.1f}%)")
            ax.legend(wedges, legend_labels, loc="center left", bbox_to_anchor=(1.02, 0.5))
            title = row["sample"] if row["sgrna"] == "NA" else f"{row['sample']} ({row['sgrna']})"
            ax.set_title(title)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def main() -> None:
    args = parse_args()
    input_root = Path(args.input_root).resolve()
    output_path = Path(args.output).resolve()
    plot_path = Path(args.plot).resolve() if args.plot else output_path.with_name(f"{output_path.stem}_stacked_bar.pdf")
    pie_path = Path(args.pie).resolve() if args.pie else output_path.with_name(f"{output_path.stem}_pies.pdf")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    pie_path.parent.mkdir(parents=True, exist_ok=True)

    if not input_root.exists():
        raise FileNotFoundError(f"Input root {input_root} does not exist.")

    tables = collect_tables(input_root)
    if not tables:
        raise FileNotFoundError(f"No Alleles_frequency_table*.txt files found under {input_root}")

    rows = []
    for table_path in sorted(tables):
        sample, sgrna, no_del, in_frame, out_frame = summarize_file(table_path)
        rows.append({
            "sample": sample,
            "sgrna": sgrna,
            "number of reads with no deletion": no_del,
            "number of reads with in-frame deletion": in_frame,
            "number of reads with out-of-frame deletion": out_frame,
        })

    df_out = pd.DataFrame(rows)
    df_out = add_percentages(df_out)

    written_path = output_path
    try:
        df_out.to_csv(output_path, sep="\t", index=False)
    except PermissionError:
        fallback = input_root / output_path.name
        fallback.parent.mkdir(parents=True, exist_ok=True)
        df_out.to_csv(fallback, sep="\t", index=False)
        written_path = fallback

    plot_written = plot_path
    try:
        plot_stacked_bars(df_out, plot_path)
    except PermissionError:
        fallback_plot = (written_path.parent / plot_path.name)
        fallback_plot.parent.mkdir(parents=True, exist_ok=True)
        plot_stacked_bars(df_out, fallback_plot)
        plot_written = fallback_plot
    pie_written = pie_path
    try:
        plot_pies(df_out, pie_path)
    except PermissionError:
        fallback_pie = (written_path.parent / pie_path.name)
        fallback_pie.parent.mkdir(parents=True, exist_ok=True)
        plot_pies(df_out, fallback_pie)
        pie_written = fallback_pie

    print(f"Wrote summary for {len(rows)} files to {written_path}")
    print(f"Wrote stacked bar plot to {plot_written}")
    print(f"Wrote pie charts to {pie_written}")


if __name__ == "__main__":
    main()
