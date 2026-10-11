"""ClassyFire count/composition plots and opposite-direction normalization.

Input statistics are the tables written by the repository's YMSTAT workflow.
No raw abundance matrix is required. All file locations are supplied in JSON;
the script contains no user-specific absolute paths.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator, PercentFormatter

P_THRESHOLD = 0.05
G_THRESHOLD = 0.7
FC_THRESHOLD = 1.5
LIPID_SUPERCLASS = "Lipids and lipid-like molecules"


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_excel(path) if path.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(path)


def find_stat_column(columns, suffix: str, first_tokens: list[str], second_tokens: list[str]) -> str:
    hits = [
        str(c) for c in columns
        if str(c).endswith(suffix)
        and any(x in str(c) for x in first_tokens)
        and any(x in str(c) for x in second_tokens)
    ]
    if len(hits) != 1:
        raise ValueError(f"Expected one {suffix} column, found: {hits}")
    return hits[0]


def clean_class(value) -> str:
    if pd.isna(value) or str(value).strip().lower() in {"", "0", "nan", "none", "na", "n/a"}:
        return "Unclassified"
    return str(value).strip()


def qualifies(p: pd.Series, g: pd.Series, fc: pd.Series) -> pd.Series:
    return ((p < P_THRESHOLD) | (g.abs() >= G_THRESHOLD)) & (
        (fc >= FC_THRESHOLD) | (fc <= 1 / FC_THRESHOLD)
    )


def load_results(config: dict, config_dir: Path):
    selected_rows, normalized_rows, normalization_summary = [], [], []
    comparisons = config["comparisons"]
    if len(comparisons) != 2:
        raise ValueError("Exactly two comparisons are required for normalization.")

    for organ, raw_path in config["statistics_files"].items():
        path = (config_dir / raw_path).resolve()
        frame = read_table(path)
        compound_col = config.get("compound_column", "Compound name")
        if compound_col not in frame:
            compound_col = str(frame.columns[0])
        tests = []
        for comp in comparisons:
            p_col = find_stat_column(frame.columns, "___u-test", comp["first_tokens"], comp["second_tokens"])
            g_col = find_stat_column(frame.columns, "___hedges_g", comp["first_tokens"], comp["second_tokens"])
            fc_col = find_stat_column(frame.columns, "___fold_change", comp["first_tokens"], comp["second_tokens"])
            p = pd.to_numeric(frame[p_col], errors="coerce")
            g = pd.to_numeric(frame[g_col], errors="coerce")
            fc = pd.to_numeric(frame[fc_col], errors="coerce")
            q = qualifies(p, g, fc)
            tests.append({"p": p, "g": g, "fc": fc, "q": q})
            for idx in frame.index[q.fillna(False)]:
                selected_rows.append({
                    "Organ": organ,
                    "Comparison": comp["label"],
                    "Compound": frame.at[idx, compound_col],
                    "Superclass": clean_class(frame.at[idx, "Superclass"]),
                    "Class": clean_class(frame.at[idx, "Class"]),
                    "Subclass": clean_class(frame.at[idx, "Subclass"]),
                    "Parent Level 1": clean_class(frame.at[idx, "Parent Level 1"]),
                    "p_value": p.at[idx], "Hedges_g": g.at[idx], "Fold_change": fc.at[idx],
                })

        opposite = (tests[0]["fc"] - 1) * (tests[1]["fc"] - 1) < 0
        normalized = tests[0]["q"] & tests[1]["q"] & opposite
        denominator = int(tests[0]["q"].sum())
        numerator = int(normalized.sum())
        normalization_summary.append({
            "Organ": organ,
            "AD_vs_WT_eligible": denominator,
            "Normalized_opposite": numerator,
            "Normalization_rate_pct": 100 * numerator / denominator if denominator else 0.0,
        })
        for idx in frame.index[normalized.fillna(False)]:
            normalized_rows.append({
                "Organ": organ,
                "Compound": frame.at[idx, compound_col],
                "AD_vs_WT_p": tests[0]["p"].at[idx],
                "AD_vs_WT_Hedges_g": tests[0]["g"].at[idx],
                "AD_vs_WT_fold_change": tests[0]["fc"].at[idx],
                "ADLP_vs_AD_p": tests[1]["p"].at[idx],
                "ADLP_vs_AD_Hedges_g": tests[1]["g"].at[idx],
                "ADLP_vs_AD_fold_change": tests[1]["fc"].at[idx],
            })
    return pd.DataFrame(selected_rows), pd.DataFrame(normalized_rows), pd.DataFrame(normalization_summary)


def save_figure(fig, output: Path, stem: str):
    fig.savefig(output / f"{stem}.svg", format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(output / f"{stem}.tiff", format="tiff", dpi=300, bbox_inches="tight",
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def superclass_plots(selected: pd.DataFrame, config: dict, output: Path):
    organs = list(config["statistics_files"])
    comparisons = [item["label"] for item in config["comparisons"]]
    counts = selected.groupby(["Comparison", "Superclass", "Organ"]).size().rename("Count").reset_index()
    keep = counts.groupby("Superclass")["Count"].max()
    retained = keep[keep >= int(config.get("minimum_superclass_count", 5))].index
    counts = counts[counts["Superclass"].isin(retained)].copy()
    order = counts.groupby("Superclass")["Count"].sum().sort_values(ascending=False).index
    colors = plt.get_cmap("tab10")(np.linspace(0, 0.8, len(organs)))
    for comparison in comparisons:
        matrix = counts[counts["Comparison"] == comparison].pivot_table(
            index="Superclass", columns="Organ", values="Count", fill_value=0
        ).reindex(index=order, columns=organs, fill_value=0)
        fig, ax = plt.subplots(figsize=(17, 10.5))
        x, total_width = np.arange(len(matrix)), 0.82
        for index, organ in enumerate(organs):
            width = total_width / len(organs)
            ax.bar(x - total_width / 2 + (index + 0.5) * width, matrix[organ],
                   width=width * 0.94, label=organ, color=colors[index])
        ax.set_ylabel("Qualifying metabolite count")
        ax.set_xticks(x, matrix.index, rotation=25, ha="right", fontweight="bold")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(axis="y", color="#D7DCE2"); ax.set_axisbelow(True)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=min(6, len(organs)), frameon=False)
        save_figure(fig, output, "superclass_" + comparison.lower().replace(" ", "_") )
    counts.to_csv(output / "superclass_counts.csv", index=False)


def lipid_composition(selected: pd.DataFrame, config: dict, output: Path):
    lipid = selected[selected["Superclass"].eq(LIPID_SUPERCLASS)].copy()
    n_top = int(config.get("top_lipid_subclasses", 5))
    top = lipid["Subclass"].value_counts().head(n_top).index.tolist()
    segments = top + ["Other"]
    rows = []
    for keys, subset in lipid.groupby(["Organ", "Comparison"]):
        counts, total = Counter(subset["Subclass"]), len(subset)
        for segment in segments:
            count = counts[segment] if segment != "Other" else total - sum(counts[x] for x in top)
            rows.append({"Organ": keys[0], "Comparison": keys[1], "Subclass": segment,
                         "Count": count, "Percent": 100 * count / total if total else 0})
    table = pd.DataFrame(rows)
    labels = [f"{o} | {c}" for o in config["statistics_files"] for c in [x["label"] for x in config["comparisons"]]]
    fig, ax = plt.subplots(figsize=(16, 12)); y = np.arange(len(labels)); left = np.zeros(len(labels))
    pastel = ["#F4A6B8", "#CDB4DB", "#A8DADC", "#FFD6A5", "#BDE0FE", "#D9DDE3"]
    for segment, color in zip(segments, pastel):
        values = []
        for label in labels:
            organ, comparison = label.split(" | ", 1)
            hit = table[(table.Organ == organ) & (table.Comparison == comparison) & (table.Subclass == segment)]
            values.append(float(hit.Percent.iloc[0]) if len(hit) else 0.0)
        ax.barh(y, values, left=left, label=segment, color=color); left += values
    ax.set_yticks(y, labels); ax.invert_yaxis(); ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(PercentFormatter()); ax.grid(axis="x", color="#D7DCE2")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False)
    save_figure(fig, output, "lipid_subclass_composition")
    table.to_csv(output / "lipid_subclass_composition.csv", index=False)


def normalization_plot(summary: pd.DataFrame, output: Path):
    table = summary.sort_values("Normalization_rate_pct")
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(table.Organ, table.Normalization_rate_pct, color="#84B6D7")
    for y, row in enumerate(table.itertuples(index=False)):
        ax.text(row.Normalization_rate_pct + 0.5, y,
                f"{row.Normalization_rate_pct:.1f}% ({row.Normalized_opposite}/{row.AD_vs_WT_eligible})",
                va="center", fontweight="bold")
    ax.set_xlabel("Normalization rate (%)")
    ax.grid(axis="x", color="#D7DCE2"); ax.set_axisbelow(True)
    save_figure(fig, output, "organ_normalization_opposite_direction")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    output = (config_path.parent / config["output_directory"]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "pdf.fonttype": 42})
    selected, normalized, summary = load_results(config, config_path.parent)
    selected.to_csv(output / "qualifying_metabolites.csv", index=False)
    normalized.to_csv(output / "normalized_metabolites.csv", index=False)
    summary.sort_values("Normalization_rate_pct", ascending=False).to_csv(
        output / "organ_normalization_summary.csv", index=False
    )
    superclass_plots(selected, config, output)
    lipid_composition(selected, config, output)
    normalization_plot(summary, output)


if __name__ == "__main__":
    main()
