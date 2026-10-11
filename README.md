# Gut–Brain Axis in Alzheimer’s Disease

Analysis code for the multi-organ metabolomics, Caco-2 cell, and western-blot
experiments used in the Alzheimer’s disease gut–brain-axis study.

## Analysis workflow

1. Missing-value filtering and imputation
2. Principal component analysis (PCA)
3. ClassyFire chemical-class summaries, stacked composition bars, and
   organ-level normalization rates
4. Univariate statistics
5. Box plots, bar plots, and volcano plots
6. Spearman/Pearson correlations, heatmaps, and scatter plots
7. Caco-2 cell one-way ANOVA with Dunnett’s comparisons
8. Western-blot one-way ANOVA with Dunnett’s comparisons

The manuscript-ready statistical description is in [METHODS.md](METHODS.md).

## Repository layout

- `R/`: reusable R functions for preprocessing, statistics, PCA, and correlation
- `Scripts/`: original end-to-end R workflows
- `Python/`: publication-figure and ANOVA/Dunnett workflows added for the paper
- `Example data/`: small example inputs only; no raw study data are distributed
- `config/`: example configuration files

## Core definitions

For metabolite selection in the ClassyFire and normalization figures, a feature
qualified when:

```text
(two-sided Mann–Whitney p < 0.05 OR |Hedges' g| >= 0.7)
AND (fold change >= 1.5 OR fold change <= 1/1.5)
```

A normalized feature qualified in both `AD-CT versus WT-CT` and
`AD-LP versus AD-CT` and showed opposite fold-change directions. The organ-level
normalization rate was the number of normalized features divided by the number
of qualifying AD-CT-versus-WT-CT features.

## Requirements

- R >= 4.0
- Python >= 3.10
- R packages listed in `R-packages.txt`
- Python packages listed in `requirements.txt`

## Reproducible use

Run the original R workflows from the repository root:

```r
source("Scripts/MultiOrgan_Total_Split.R")
source("Scripts/Visualization.R")
source("Scripts/Correlation.R")
```

The Python scripts expose command-line help:

```bash
python Python/classyfire_normalization.py --help
python Python/anova_dunnett.py --help
```

Input files must use de-identified sample IDs. Large or sensitive raw data,
local absolute paths, and generated figures should not be committed.

## Citation

Use the repository’s **Cite this repository** menu, which is populated by
`CITATION.cff`. For a permanent manuscript citation, archive a tagged GitHub
release in Zenodo and cite the resulting DOI.

## License

MIT License. See [LICENSE](LICENSE).

