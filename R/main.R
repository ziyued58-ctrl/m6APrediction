#' @title Encode DNA Sequences as Nucleotide Factors
#'
#' @description
#' Split equal-length DNA sequences into nucleotide-position columns.
#' Each column is a categorical factor with levels A, T, C and G,
#' in that order.
#'
#' @param dna_strings A non-empty character vector of equal-length
#'   uppercase DNA sequences containing only A, T, C and G.
#'
#' @return A data frame with one row per sequence and columns named
#'   nt_pos1, nt_pos2 and so on. Each column is a factor.
#'
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "GAACT"))
#'
#' @keywords internal
dna_encoding <- function(dna_strings) {
  stopifnot(
    is.character(dna_strings),
    length(dna_strings) > 0L,
    !anyNA(dna_strings),
    all(grepl("^[ATCG]+$", dna_strings)),
    length(unique(nchar(dna_strings))) == 1L
  )

  nn <- nchar(dna_strings[1])

  seq_m <- matrix(
    unlist(strsplit(dna_strings, "")),
    ncol = nn,
    byrow = TRUE
  )

  colnames(seq_m) <- paste0("nt_pos", seq_len(nn))

  seq_df <- as.data.frame(seq_m)

  seq_df[] <- lapply(
    seq_df,
    factor,
    levels = c("A", "T", "C", "G")
  )

  return(seq_df)
}


#' @title Predict m6A Sites for Multiple Observations
#'
#' @description
#' Use a fitted Random Forest classifier to predict the Positive-class
#' probability and m6A status for each input row. Categorical factor
#' levels are set to match the Practical 5 training procedure.
#'
#' @param ml_fit A fitted randomForest classification model with
#'   classes Negative and Positive.
#' @param feature_df A non-empty data frame containing gc_content,
#'   RNA_type, RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation and DNA_5mer. DNA_5mer must contain
#'   uppercase five-nucleotide sequences. Numeric features must use
#'   the same scales as the training data.
#' @param positive_threshold A single numeric value from 0 to 1.
#'   A row is Positive when its probability is strictly greater
#'   than this threshold. The default is 0.5.
#'
#' @return The input data frame with categorical columns converted
#'   to the required factors and two prediction columns:
#'   predicted_m6A_prob, a numeric Positive-class probability, and
#'   predicted_m6A_status, either Positive or Negative.
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' model_path <- system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' )
#' input_path <- system.file(
#'   "extdata", "m6A_input_example.csv", package = "m6APrediction"
#' )
#' model <- readRDS(model_path)
#' input_data <- utils::read.csv(input_path)
#' predictions <- prediction_multiple(model, input_data)
#' head(predictions)
prediction_multiple <- function(
    ml_fit,
    feature_df,
    positive_threshold = 0.5) {

  required_columns <- c(
    "gc_content",
    "RNA_type",
    "RNA_region",
    "exon_length",
    "distance_to_junction",
    "evolutionary_conservation",
    "DNA_5mer"
  )

  stopifnot(
    is.data.frame(feature_df),
    nrow(feature_df) > 0L,
    all(required_columns %in% colnames(feature_df)),
    is.numeric(positive_threshold),
    length(positive_threshold) == 1L,
    is.finite(positive_threshold),
    positive_threshold >= 0,
    positive_threshold <= 1
  )

  stopifnot(
    !anyNA(feature_df[, required_columns, drop = FALSE]),
    all(grepl(
      "^[ATCG]{5}$",
      as.character(feature_df$DNA_5mer)
    )),
    all(feature_df$RNA_type %in% c(
      "mRNA", "lincRNA", "lncRNA", "pseudogene"
    )),
    all(feature_df$RNA_region %in% c(
      "CDS", "intron", "3'UTR", "5'UTR"
    ))
  )

  feature_df$RNA_type <- factor(
    feature_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )

  feature_df$RNA_region <- factor(
    feature_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )

  model_input <- cbind(
    feature_df[, setdiff(required_columns, "DNA_5mer"),
               drop = FALSE],
    dna_encoding(as.character(feature_df$DNA_5mer))
  )

  inferred_prob <- stats::predict(
    ml_fit,
    newdata = model_input,
    type = "prob"
  )

  feature_df$predicted_m6A_prob <- inferred_prob[, "Positive"]

  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > positive_threshold,
    "Positive",
    "Negative"
  )

  return(feature_df)
}


#' @title Predict an m6A Site for One Observation
#'
#' @description
#' Construct a one-row feature data frame and call
#' prediction_multiple() to predict the m6A probability and status.
#'
#' @param ml_fit A fitted randomForest classification model with
#'   classes Negative and Positive.
#' @param gc_content A single numeric GC-content value from 0 to 1.
#' @param RNA_type One of mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region One of CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length A single numeric exon-length feature on the
#'   same scale as the training data.
#' @param distance_to_junction A single numeric junction-distance
#'   feature on the same scale as the training data.
#' @param evolutionary_conservation A single numeric conservation
#'   value from 0 to 1.
#' @param DNA_5mer A single uppercase DNA sequence of length five,
#'   containing only A, T, C and G.
#' @param positive_threshold A single numeric value from 0 to 1.
#'   Positive is assigned when probability is strictly greater
#'   than this threshold. The default is 0.5.
#'
#' @return A named character vector of length two, with names
#'   predicted_m6A_prob and predicted_m6A_status. The probability
#'   is converted to character because the vector also contains
#'   the character status.
#'
#' @export
#'
#' @examples
#' model_path <- system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' )
#' model <- readRDS(model_path)
#' prediction_single(
#'   ml_fit = model,
#'   gc_content = 0.5,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 10,
#'   distance_to_junction = 8,
#'   evolutionary_conservation = 0.5,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
prediction_single <- function(
    ml_fit,
    gc_content,
    RNA_type,
    RNA_region,
    exon_length,
    distance_to_junction,
    evolutionary_conservation,
    DNA_5mer,
    positive_threshold = 0.5) {

  inputs <- list(
    gc_content,
    RNA_type,
    RNA_region,
    exon_length,
    distance_to_junction,
    evolutionary_conservation,
    DNA_5mer
  )

  stopifnot(all(lengths(inputs) == 1L))

  one_row <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  pred_outcome <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = one_row,
    positive_threshold = positive_threshold
  )

  returned_vector <- c(
    predicted_m6A_prob = pred_outcome$predicted_m6A_prob[1],
    predicted_m6A_status = pred_outcome$predicted_m6A_status[1]
  )

  return(returned_vector)
}
