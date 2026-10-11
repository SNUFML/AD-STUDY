# Statistical analysis

The text below is written for direct adaptation to the Statistical analysis
subsection of the manuscript. Group labels are retained as used in the analysis:
WT-CT, AD-CT, and AD-LP.

## Manuscript-ready text

Metabolomics data were processed separately for each organ or specimen type.
For each experimental group, the proportion of missing observations was
calculated for every metabolite. A metabolite was retained for an organ when it
had no more than 20% missing observations in at least one group from that organ.
After organ-specific filtering, remaining missing values were replaced with
one-fifth of the smallest observed positive abundance in the numeric data block
for that organ. No additional transformation was applied in the original
preprocessing script before the group comparisons described below.

Principal component analysis (PCA) was performed on the metabolite abundance
matrix after mean centering and unit-variance scaling using the R function
`prcomp` with `scale = TRUE`. PCA scores were displayed for the first two
components, and 68% group ellipses were added for visualization. The percentages
of variance explained were calculated from the eigenvalues of the scaled PCA.
When PC-score comparisons were displayed, the same predefined group contrasts
were evaluated on the PC scores.

For univariate metabolite comparisons, two prespecified two-group contrasts
were examined: AD-CT versus WT-CT and AD-LP versus AD-CT. Two-sided
Mann–Whitney U tests were used for the primary nonparametric comparisons.
Hedges’ g was calculated from the pooled standard deviation with the small-sample
correction, and fold change was calculated as the arithmetic mean of the first
group divided by that of the second group. Unless specifically stated otherwise,
the metabolomics screening p values shown in the figures were unadjusted. A
metabolite was included in the ClassyFire and normalization summaries when it
satisfied `(p < 0.05 or |Hedges’ g| >= 0.7)` and also showed an absolute
fold-change threshold of at least 1.5 (fold change >= 1.5 or <= 1/1.5).

Chemical taxonomy was obtained from the ClassyFire superclass, class, subclass,
and parent-level annotations supplied with the metabolite table. Superclass bar
plots show the number of qualifying metabolites in each organ and comparison.
For display, a superclass was retained when at least one organ contained five or
more qualifying metabolites in either comparison. For lipid composition plots,
the five most frequent lipid subclasses across organs and comparisons were shown
individually, remaining subclasses were combined as “Other,” and each bar was
normalized to 100% of qualifying lipid metabolites for that organ and comparison.

A metabolite was defined as showing a normalized pattern when it met the stated
selection criterion in both AD-CT versus WT-CT and AD-LP versus AD-CT and the two
fold changes were in opposite directions. The normalization rate for each organ
was calculated as the number of normalized-pattern metabolites divided by the
number of metabolites qualifying in the AD-CT-versus-WT-CT comparison, multiplied
by 100. Metabolite box plots show the median, interquartile range, whiskers, and
individual observations on the original relative-abundance scale. Statistical
annotations were based on the corresponding two-sided Mann–Whitney U tests; for
nonsignificant comparisons, the exact p value and Hedges’ g were reported.

Associations between serum taurine-conjugated bile acids and bile acids measured
in colon, cecum, and liver were assessed in matched animals. Spearman rank and
Pearson product–moment correlations were calculated using pairwise complete
observations; pairs with fewer than three complete observations or with a
constant variable were not tested. Two-sided raw p values were reported. For
heatmap visualization, columns were hierarchically clustered separately within
each organ using Euclidean distances between correlation profiles and average
linkage. Blue and red indicate negative and positive correlations, respectively.
Scatter plots display individual animals, a fitted linear trend, and its 95%
confidence interval, and report both Spearman and Pearson coefficients and p
values. Because these correlation analyses screened multiple metabolite pairs,
they should be described as exploratory unless a prespecified multiplicity
correction is added.

For Caco-2 cell experiments, each outcome was summarized as mean ± standard
deviation with individual replicate values displayed. Overall group differences
were evaluated using a one-way analysis of variance (ANOVA) assuming equal
within-group variances. When treatment groups were compared with the prespecified
INF control, two-sided Dunnett tests were used, and Dunnett-adjusted p values were
reported. The same analysis was used for western-blot densitometry outcomes,
including the PXR and p-p65/p65 measurements: one-way ANOVA followed by two-sided
Dunnett comparisons of each treatment group with INF. Statistical significance
was defined as adjusted p < 0.05 for these Dunnett comparisons. Bar plots show
mean ± standard deviation and all individual replicate values.

Analyses were implemented in R (version 4.0 or later) and Python (version 3.10 or
later). The R workflow uses `matrixTests` for vectorized two-sample Wilcoxon
tests, base R `prcomp` for PCA, and `factoextra` for PCA visualization. Python
analyses use NumPy, pandas, SciPy, and Matplotlib. All tests were two-sided.

## Reporting notes before submission

- Replace the generic software versions with the exact versions used for the
  final locked analysis environment.
- The metabolomics screening figures use raw p values. If the journal requires
  false-discovery-rate control, add Benjamini–Hochberg q values without silently
  replacing the prespecified selection rule.
- State the biological and technical replicate numbers in each figure legend.
- The original missing-value code uses the minimum of the complete organ-level
  numeric block divided by five. Confirm that all retained abundances are
  positive before applying this rule.
