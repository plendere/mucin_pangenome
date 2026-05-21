#!/usr/bin/env Rscript
# =============================================================================
# Estimate the number of reference genomes required to reach 95%
# leave-one-out genotyping concordance per locus
# =============================================================================
#
# Description
# -----------
# Given leave-one-out genotyping concordance measured at several reference
# panel sizes (a subsampling experiment), this script fits a logistic
# saturation model independently for each locus and each concordance metric,
# and extrapolates the number of reference genomes required to reach 95%
# concordance.
#
# Model
# -----
#   concordance(n) = A / (1 + exp(-k * (n - n0)))
#
#   where
#     A  = asymptotic maximum concordance (the ceiling the curve approaches)
#     k  = growth rate (steepness of the rise)
#     n0 = inflection point (panel size at which concordance rises fastest)
#
# A logistic model was chosen because it (i) is naturally bounded between 0
# and 1, matching a concordance metric, (ii) accommodates an S-shape with a
# slow start followed by saturation, and (iii) provides an explicit asymptote
# (A), which makes it possible to distinguish loci that will plateau below
# 95% concordance from loci that will eventually reach it given more genomes.
#
# For each locus x metric the script reports one of three outcomes:
#   - "Already achieved": current observed concordance >= 0.95.
#   - "Achievable":       fitted asymptote A >= 0.95; n95 is the panel size
#                         at which the predicted curve crosses 0.95.
#   - "Plateau":          fitted asymptote A <  0.95; locus is predicted not
#                         to reach 95% concordance even with more genomes.
#
# Input
# -----
# A long-format TSV/CSV with one row per (locus, metric, panel size, replicate)
# observation. Required columns:
#
#   locus        character  name of the locus / gene
#   metric       character  name of the concordance metric (e.g. "bp",
#                           "haplogroup"). Any number of metrics is allowed;
#                           a separate fit is performed for each.
#   n_genomes    numeric    reference panel size for this observation
#   concordance  numeric    concordance value in [0, 1] for this observation
#
# Replicate observations at the same (locus, metric, n_genomes) are allowed
# and are used as-is by the fit (each contributes one data point).
#
# Output
# ------
# 1. <prefix>_fit_results.tsv : per (locus, metric) fitted parameters,
#                               classification, and n95 estimate.
# 2. <prefix>_predictions.tsv : fitted curves on a dense grid (for plotting).
# 3. <prefix>_curves.pdf      : per-metric facet plot of observations,
#                               fitted curves, and the 95% threshold.
#
# Usage
# -----
#   Rscript estimate_genomes_for_95pct_concordance.R \
#       --input  concordance_long.tsv \
#       --output_prefix results/genomes_for_95pct \
#       [--threshold 0.95] \
#       [--max_predict 2000] \
#       [--current_panel 276]
#
# Dependencies: dplyr, tidyr, ggplot2, optparse, minpack.lm (for nlsLM)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(minpack.lm)
})

# -----------------------------------------------------------------------------
# Command-line arguments
# -----------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Long-format TSV/CSV with columns: locus, metric, n_genomes, concordance"),
  make_option(c("-o", "--output_prefix"), type = "character", default = "genomes_for_95pct",
              help = "Prefix (including directory) for output files [default: %default]"),
  make_option(c("-t", "--threshold"), type = "double", default = 0.95,
              help = "Concordance threshold to extrapolate to [default: %default]"),
  make_option(c("-m", "--max_predict"), type = "double", default = 2000,
              help = "Upper bound (in number of genomes) for the prediction grid [default: %default]"),
  make_option(c("-c", "--current_panel"), type = "double", default = NA,
              help = "Optional: current reference panel size to draw on the plot")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Please provide an input file with --input.")
}

# -----------------------------------------------------------------------------
# Load and validate input
# -----------------------------------------------------------------------------
read_table_auto <- function(path) {
  ext <- tools::file_ext(path)
  if (tolower(ext) == "csv") {
    read.csv(path, stringsAsFactors = FALSE)
  } else {
    read.delim(path, stringsAsFactors = FALSE)
  }
}

dat <- read_table_auto(opt$input)

required_cols <- c("locus", "metric", "n_genomes", "concordance")
missing_cols <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop("Input is missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat <- dat %>%
  mutate(
    locus       = as.character(locus),
    metric      = as.character(metric),
    n_genomes   = as.numeric(n_genomes),
    concordance = as.numeric(concordance)
  ) %>%
  filter(!is.na(n_genomes), !is.na(concordance))

cat("Loaded", nrow(dat), "observations across",
    length(unique(dat$locus)), "loci and",
    length(unique(dat$metric)), "metric(s).\n")

# -----------------------------------------------------------------------------
# Logistic fit per (locus, metric)
# -----------------------------------------------------------------------------
# concordance = A / (1 + exp(-k * (n - n0)))
fit_logistic <- function(df, threshold = 0.95, max_predict = 2000) {

  # Need a reasonable number of distinct panel sizes to fit 3 parameters
  if (length(unique(df$n_genomes)) < 4) {
    return(tibble(
      A = NA_real_, k = NA_real_, n0 = NA_real_,
      n95 = NA_real_, current_max = max(df$concordance, na.rm = TRUE),
      status = "insufficient_data"
    ))
  }

  current_max <- max(df$concordance, na.rm = TRUE)

  if (current_max >= threshold) {
    return(tibble(
      A = NA_real_, k = NA_real_, n0 = NA_real_,
      n95 = NA_real_, current_max = current_max,
      status = "already_achieved"
    ))
  }

  start_A  <- min(max(current_max, 0.5), 1.0)
  start_k  <- 0.01
  start_n0 <- median(df$n_genomes, na.rm = TRUE)

  fit <- tryCatch(
    nlsLM(
      concordance ~ A / (1 + exp(-k * (n_genomes - n0))),
      data    = df,
      start   = list(A = start_A, k = start_k, n0 = start_n0),
      lower   = c(0,   1e-6, 0),
      upper   = c(1,   10,   max_predict * 5),
      control = nls.lm.control(maxiter = 200)
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      A = NA_real_, k = NA_real_, n0 = NA_real_,
      n95 = NA_real_, current_max = current_max,
      status = "fit_failed"
    ))
  }

  cf <- coef(fit)
  A  <- unname(cf["A"])
  k  <- unname(cf["k"])
  n0 <- unname(cf["n0"])

  if (A < threshold) {
    return(tibble(
      A = A, k = k, n0 = n0,
      n95 = NA_real_, current_max = current_max,
      status = "plateau"
    ))
  }

  # Solve A / (1 + exp(-k(n - n0))) = threshold for n
  #     => n = n0 - log(A/threshold - 1) / k
  n95 <- n0 - log(A / threshold - 1) / k

  if (!is.finite(n95) || n95 <= 0) {
    return(tibble(
      A = A, k = k, n0 = n0,
      n95 = NA_real_, current_max = current_max,
      status = "fit_failed"
    ))
  }

  tibble(
    A = A, k = k, n0 = n0,
    n95 = n95, current_max = current_max,
    status = "achievable"
  )
}

fit_results <- dat %>%
  group_by(locus, metric) %>%
  group_modify(~ fit_logistic(.x, threshold = opt$threshold,
                              max_predict = opt$max_predict)) %>%
  ungroup() %>%
  mutate(
    threshold = opt$threshold,
    n95_rounded = ifelse(is.na(n95), NA_integer_, ceiling(n95))
  )

# -----------------------------------------------------------------------------
# Build prediction grid for plotting
# -----------------------------------------------------------------------------
predict_logistic <- function(n, A, k, n0) {
  A / (1 + exp(-k * (n - n0)))
}

grid_n <- seq(1, opt$max_predict, length.out = 400)

pred_list <- list()
for (i in seq_len(nrow(fit_results))) {
  r <- fit_results[i, ]
  if (is.na(r$A) || is.na(r$k) || is.na(r$n0)) next
  pred_list[[length(pred_list) + 1]] <- tibble(
    locus       = r$locus,
    metric      = r$metric,
    n_genomes   = grid_n,
    concordance = predict_logistic(grid_n, r$A, r$k, r$n0),
    status      = r$status
  )
}
pred_df <- if (length(pred_list) > 0) bind_rows(pred_list) else tibble()

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------
out_dir <- dirname(opt$output_prefix)
if (nchar(out_dir) > 0 && !dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

fit_out  <- paste0(opt$output_prefix, "_fit_results.tsv")
pred_out <- paste0(opt$output_prefix, "_predictions.tsv")
pdf_out  <- paste0(opt$output_prefix, "_curves.pdf")

write.table(fit_results, fit_out, sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(pred_df) > 0) {
  write.table(pred_df, pred_out, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("Wrote:", fit_out, "\n")
if (nrow(pred_df) > 0) cat("Wrote:", pred_out, "\n")

# -----------------------------------------------------------------------------
# Diagnostic plot
# -----------------------------------------------------------------------------
p <- ggplot() +
  geom_point(data = dat,
             aes(x = n_genomes, y = concordance, color = locus),
             alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = opt$threshold, linetype = "dashed", color = "grey40") +
  facet_wrap(~ metric, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x     = "Number of reference genomes",
    y     = "Leave-one-out concordance",
    color = "Locus",
    title = sprintf("Logistic extrapolation to %g%% concordance",
                    100 * opt$threshold)
  ) +
  theme_bw() +
  theme(legend.position = "right")

if (nrow(pred_df) > 0) {
  p <- p + geom_line(data = pred_df,
                     aes(x = n_genomes, y = concordance, color = locus),
                     alpha = 0.7)
}

if (!is.na(opt$current_panel)) {
  p <- p + geom_vline(xintercept = opt$current_panel,
                      linetype = "dotted", color = "grey40")
}

ggsave(pdf_out, p, width = 9, height = 5)
cat("Wrote:", pdf_out, "\n")

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------
summary_tbl <- fit_results %>%
  mutate(
    label = case_when(
      status == "already_achieved" ~ sprintf("already >= %.0f%%",  100 * opt$threshold),
      status == "achievable"       ~ sprintf("~%d genomes",        n95_rounded),
      status == "plateau"          ~ sprintf("plateau at ~%.1f%%", 100 * A),
      status == "fit_failed"       ~ "fit failed",
      status == "insufficient_data"~ "insufficient data",
      TRUE                         ~ status
    )
  ) %>%
  select(locus, metric, status, current_max, A, n95_rounded, label)

cat("\nPer-locus summary:\n")
print(as.data.frame(summary_tbl), row.names = FALSE)


