"""Matched cross-organ Spearman/Pearson correlations and scatter plots.

The annotation table must contain ``Name`` and ``Parentclass``. Input matrices
must contain the same sample-key column followed by numeric metabolite columns.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_excel(path) if path.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(path)


def safe_stem(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_")[:120]


def correlation(x: np.ndarray, y: np.ndarray, method: str):
    keep = np.isfinite(x) & np.isfinite(y)
    n = int(keep.sum())
    if n < 3 or np.unique(x[keep]).size < 2 or np.unique(y[keep]).size < 2:
        return np.nan, np.nan, n
    result = stats.spearmanr(x[keep], y[keep]) if method == "spearman" else stats.pearsonr(x[keep], y[keep])
    return float(result.statistic), float(result.pvalue), n


def fit_line(x: np.ndarray, y: np.ndarray):
    keep = np.isfinite(x) & np.isfinite(y); x, y = x[keep], y[keep]
    slope, intercept, _, _, _ = stats.linregress(x, y)
    grid = np.linspace(x.min(), x.max(), 100); fitted = intercept + slope * grid
    if len(x) <= 2:
        return grid, fitted, fitted, fitted
    residual = y - (intercept + slope * x)
    s_err = np.sqrt(np.sum(residual**2) / (len(x) - 2))
    center = np.sum((x - x.mean()) ** 2)
    se = s_err * np.sqrt(1 / len(x) + (grid - x.mean()) ** 2 / center)
    critical = stats.t.ppf(0.975, len(x) - 2)
    return grid, fitted, fitted - critical * se, fitted + critical * se


def heatmap(table: pd.DataFrame, organ: str, method: str, output: Path):
    coefficient = "Spearman_rho" if method == "spearman" else "Pearson_r"
    p_column = "Spearman_p" if method == "spearman" else "Pearson_p"
    matrix = table.pivot(index="Serum_metabolite", columns="Organ_metabolite", values=coefficient)
    pvalues = table.pivot(index="Serum_metabolite", columns="Organ_metabolite", values=p_column).reindex_like(matrix)
    if matrix.shape[1] > 1:
        prepared = np.nan_to_num(matrix.to_numpy().T, nan=0.0)
        distances = pdist(prepared, metric="euclidean")
        if not np.allclose(distances, 0):
            order = leaves_list(linkage(distances, method="average", optimal_ordering=True))
            matrix = matrix.iloc[:, order]; pvalues = pvalues.iloc[:, order]
    fig_width = max(8, matrix.shape[1] * 0.55)
    fig_height = max(6, matrix.shape[0] * 0.55)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto", interpolation="none")
    ax.set_xticks(np.arange(matrix.shape[1]), matrix.columns, rotation=55, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]), matrix.index)
    ax.set_xlabel(f"{organ} bile acids", fontweight="bold")
    ax.set_ylabel("Serum taurinated bile acids", fontweight="bold")
    ax.set_title(f"{method.title()} correlation", fontweight="bold")
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            p = pvalues.iat[row, col]
            marker = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            if marker:
                ax.text(col, row, marker, ha="center", va="center", color="black", fontweight="bold")
    fig.colorbar(image, ax=ax, label="Correlation coefficient", shrink=0.8)
    stem = safe_stem(f"{organ}_{method}_heatmap")
    fig.savefig(output / f"{stem}.svg", format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(output / f"{stem}.tiff", format="tiff", dpi=300, bbox_inches="tight",
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serum", required=True, type=Path)
    parser.add_argument("--organs", required=True, nargs="+", help="Organ=path pairs")
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--sample-key", default="sample_key")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "pdf.fonttype": 42})

    serum = read_table(args.serum)
    annotations = read_table(args.annotations).dropna(subset=["Name", "Parentclass"])
    serum_names = annotations[
        annotations.Name.isin(serum.columns)
        & annotations.Parentclass.eq("Taurinated bile acids and derivatives")
    ].Name.drop_duplicates().tolist()
    if not serum_names:
        raise ValueError("No serum taurinated bile acids matched the annotations")

    all_rows = []
    for organ_spec in args.organs:
        organ_name, path_text = organ_spec.split("=", 1)
        organ = read_table(Path(path_text))
        organ_names = annotations[
            annotations.Name.isin(organ.columns)
            & annotations.Parentclass.astype(str).str.contains("bile acids", case=False, na=False)
        ].Name.drop_duplicates().tolist()
        merged = serum[[args.sample_key, *serum_names]].merge(
            organ[[args.sample_key, *organ_names]], on=args.sample_key, how="inner", validate="one_to_one"
        )
        for serum_name in serum_names:
            y = pd.to_numeric(merged[serum_name], errors="coerce").to_numpy(float)
            for organ_metabolite in organ_names:
                x = pd.to_numeric(merged[organ_metabolite], errors="coerce").to_numpy(float)
                rho, p_s, n_s = correlation(x, y, "spearman")
                r, p_p, n_p = correlation(x, y, "pearson")
                all_rows.append({"Organ": organ_name, "Serum_metabolite": serum_name,
                                 "Organ_metabolite": organ_metabolite, "n": min(n_s, n_p),
                                 "Spearman_rho": rho, "Spearman_p": p_s,
                                 "Pearson_r": r, "Pearson_p": p_p})
                if (np.isfinite(p_s) and p_s < args.alpha) or (np.isfinite(p_p) and p_p < args.alpha):
                    keep = np.isfinite(x) & np.isfinite(y)
                    grid, fitted, lower, upper = fit_line(x, y)
                    fig, ax = plt.subplots(figsize=(7, 7))
                    ax.scatter(x[keep], y[keep], s=44, color="#527AA3", edgecolor="white", linewidth=0.4)
                    ax.plot(grid, fitted, color="#B2182B", linewidth=2)
                    ax.fill_between(grid, lower, upper, color="#B2182B", alpha=0.15, linewidth=0)
                    ax.set_xlabel(organ_metabolite, fontweight="bold")
                    ax.set_ylabel(serum_name, fontweight="bold")
                    ax.set_title(f"{organ_name}\nSpearman ρ={rho:.2f}, p={p_s:.3g}; Pearson r={r:.2f}, p={p_p:.3g}")
                    ax.spines[["top", "right"]].set_visible(False)
                    stem = safe_stem(f"{organ_name}_{organ_metabolite}_vs_{serum_name}")
                    fig.savefig(args.output / f"{stem}.svg", format="svg", bbox_inches="tight", facecolor="white")
                    fig.savefig(args.output / f"{stem}.tiff", format="tiff", dpi=300, bbox_inches="tight",
                                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
                    plt.close(fig)
    result = pd.DataFrame(all_rows)
    result.to_csv(args.output / "correlation_results.csv", index=False)
    for organ_name, organ_table in result.groupby("Organ", sort=False):
        heatmap(organ_table, organ_name, "spearman", args.output)
        heatmap(organ_table, organ_name, "pearson", args.output)


if __name__ == "__main__":
    main()
