"""One-way ANOVA followed by two-sided Dunnett comparisons with a control."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def stars(p: float) -> str:
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def analyze(data: pd.DataFrame, group_col: str, value_col: str, control: str, order: list[str]):
    arrays = {group: pd.to_numeric(data.loc[data[group_col] == group, value_col], errors="coerce").dropna().to_numpy() for group in order}
    if control not in arrays:
        raise ValueError(f"Control {control!r} is not in group order")
    if any(len(values) < 2 for values in arrays.values()):
        raise ValueError("Each group requires at least two finite observations")
    anova = stats.f_oneway(*(arrays[group] for group in order), equal_var=True)
    treatment_names = [group for group in order if group != control]
    dunnett = stats.dunnett(*(arrays[group] for group in treatment_names), control=arrays[control],
                            alternative="two-sided", rng=np.random.default_rng(20261011))
    pvalues = dict(zip(treatment_names, np.atleast_1d(dunnett.pvalue)))
    statistics = dict(zip(treatment_names, np.atleast_1d(dunnett.statistic)))
    rows = []
    for group in order:
        values = arrays[group]
        rows.append({
            "group": group, "n": len(values), "mean": np.mean(values), "sd": np.std(values, ddof=1),
            "anova_F": anova.statistic, "anova_df_between": len(order) - 1,
            "anova_df_within": sum(map(len, arrays.values())) - len(order), "anova_p": anova.pvalue,
            "dunnett_statistic": statistics.get(group, np.nan),
            "dunnett_adjusted_p": pvalues.get(group, np.nan),
            "significance": "control" if group == control else stars(float(pvalues[group])),
        })
    return arrays, pd.DataFrame(rows)


def plot(arrays, results: pd.DataFrame, value_label: str, output: Path, stem: str):
    order = list(arrays); means = [np.mean(arrays[x]) for x in order]
    sds = [np.std(arrays[x], ddof=1) for x in order]; x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.bar(x, means, yerr=sds, capsize=5, width=0.62, color="#B9D8EA", edgecolor="#4C89A8")
    for index, group in enumerate(order):
        rng = np.random.default_rng(1100 + index)
        ax.scatter(index + rng.uniform(-0.09, 0.09, len(arrays[group])), arrays[group],
                   color="#254F65", edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(x, order, fontweight="bold"); ax.set_ylabel(value_label)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(False)
    fig.savefig(output / f"{stem}.svg", format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(output / f"{stem}.tiff", format="tiff", dpi=300, bbox_inches="tight",
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--group-column", default="Group")
    parser.add_argument("--value-column", required=True)
    parser.add_argument("--control", default="INF")
    parser.add_argument("--order", required=True, help="Comma-separated group order")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--stem", default="anova_dunnett")
    args = parser.parse_args()
    frame = pd.read_excel(args.input) if args.input.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(args.input)
    order = [item.strip() for item in args.order.split(",")]
    args.output.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "pdf.fonttype": 42})
    arrays, results = analyze(frame, args.group_column, args.value_column, args.control, order)
    results.to_csv(args.output / f"{args.stem}_statistics.csv", index=False)
    plot(arrays, results, args.value_column, args.output, args.stem)


if __name__ == "__main__":
    main()
