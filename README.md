# m6APrediction

An R package for predicting N6-methyladenosine (m6A) sites using a
fitted Random Forest classifier. Developed for BIO215 Practicals 6–7.

## Features

- Batch prediction from a feature data frame.
- Single-observation prediction from manually specified features.
- Positive-class probabilities and threshold-based classification.
- Bundled fitted model and a 100-row example input table.

## Installation

Run the following commands in R:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("ziyued58-ctrl/m6APrediction")
library(m6APrediction)
```

The package declares R >= 3.5.0 and imports `randomForest` and `stats`.
Dependencies may require a newer R version; development checks were
performed with R 4.5.1.

## Required input features

| Column | Expected input |
|---|---|
| `gc_content` | Numeric GC proportion between 0 and 1 |
| `RNA_type` | `mRNA`, `lincRNA`, `lncRNA` or `pseudogene` |
| `RNA_region` | `CDS`, `intron`, `3'UTR` or `5'UTR` |
| `exon_length` | Numeric feature on the course dataset's log2 scale |
| `distance_to_junction` | Numeric feature on the course dataset's log2 scale |
| `evolutionary_conservation` | Numeric conservation score between 0 and 1 |
| `DNA_5mer` | Five uppercase nucleotides using A, T, C and G |

Use the same feature scales as the training data. The supplied model
was trained on five-base contexts with A in the central position.
Nucleotide positions are encoded as categorical factors with levels
`A`, `T`, `C`, `G`, rather than one-hot columns.

## Batch prediction

This example loads its own model and input data.

```r
library(m6APrediction)

ml_fit <- readRDS(system.file(
  "extdata", "rf_fit.rds",
  package = "m6APrediction", mustWork = TRUE
))

example_df <- read.csv(system.file(
  "extdata", "m6A_input_example.csv",
  package = "m6APrediction", mustWork = TRUE
))

batch_result <- prediction_multiple(
  ml_fit = ml_fit,
  feature_df = example_df[1:3, ],
  positive_threshold = 0.6
)

batch_result[, c(
  "predicted_m6A_prob",
  "predicted_m6A_status"
)]
```

Expected output with the bundled resources:

```text
  predicted_m6A_prob predicted_m6A_status
1              0.876             Positive
2              0.000             Negative
3              0.498             Negative
```

The returned data frame preserves the original input columns and adds
`predicted_m6A_prob` and `predicted_m6A_status`. Categorical input
columns are converted to factors matching the model.

## Single-observation prediction

This example independently loads its own model.

```r
library(m6APrediction)

ml_fit <- readRDS(system.file(
  "extdata", "rf_fit.rds",
  package = "m6APrediction", mustWork = TRUE
))

prediction_single(
  ml_fit = ml_fit,
  gc_content = 0.6,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 12,
  distance_to_junction = 5,
  evolutionary_conservation = 0.8,
  DNA_5mer = "GGACA",
  positive_threshold = 0.5
)
```

Expected output:

```text
predicted_m6A_prob predicted_m6A_status
           "0.51"           "Positive"
```

The result is a named character vector of length two. Its probability
is converted to character because the vector also contains a text status.

## Decision threshold

A site is classified as `Positive` when its predicted probability is
strictly greater than `positive_threshold`. Otherwise, including equality,
it is classified as `Negative`. The default threshold is 0.5.

The probability is the Random Forest Positive-class vote proportion;
it is not experimental confirmation of an m6A site.

## Documentation

```r
?prediction_multiple
?prediction_single
```

`dna_encoding()` is an internal helper and is not exported.

## Package validation

The local package passed:

```r
devtools::check(args = "--no-manual", cran = FALSE)
```

Result: **0 errors, 0 warnings and 0 notes**.

Both prediction examples ran successfully. PDF manual generation was
excluded from this check. The example data demonstrate software use
and are not an independent test set for assessing prediction accuracy.

## Course resources and license

The fitted model and example input data were supplied through the
BIO215 Practical 6 materials and reused without retraining.

See `LICENSE` and `LICENSE.md` for the package code's MIT license.
The bundled course resources remain subject to their original terms.
