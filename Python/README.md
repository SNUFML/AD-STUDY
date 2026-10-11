# Python paper workflows

All scripts write editable SVG files, 300-dpi LZW-compressed TIFF files, and
machine-readable CSV statistics. They do not contain local absolute paths.

## ClassyFire and normalization figures

Copy `config/paper_pipeline.example.json`, update only the input paths, and run:

```bash
python Python/classyfire_normalization.py --config config/paper_pipeline.json
```

The selection rule is:

```text
(p < 0.05 OR |Hedges' g| >= 0.7) AND (FC >= 1.5 OR FC <= 0.667)
```

## Cross-organ bile-acid correlations

```bash
python Python/correlation_analysis.py \
  --serum data/processed/serum.csv \
  --organs Colon=data/processed/colon.csv Cecum=data/processed/cecum.csv Liver=data/processed/liver.csv \
  --annotations data/processed/class_annotations.csv \
  --sample-key sample_key \
  --output outputs/correlations
```

The script calculates both Spearman and Pearson correlations, writes a complete
results table, generates an organ-specific heatmap for each method, and creates
scatter plots for pairs with raw two-sided p < 0.05 by either method.

## Caco-2 and western-blot ANOVA/Dunnett analysis

The same implementation is used for either experiment. The input is a tidy CSV
or Excel table with one row per replicate.

```bash
python Python/anova_dunnett.py \
  --input data/processed/caco2.csv \
  --group-column Group \
  --value-column Outcome \
  --control INF \
  --order "CON,INF,7KDCA,7KLCA,12KDCA,12KCDCA" \
  --output outputs/caco2 \
  --stem caco2_outcome
```

For multiple outcomes, call the script once per outcome. Dunnett p values are
already adjusted for the family of treatment-versus-control comparisons within
that outcome.
