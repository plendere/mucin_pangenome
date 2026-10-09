#!/usr/bin/env Rscript
# =============================================================================
# Association of mucin VNTR variation with a continuous phenotype
# (developed for KNoRMA scores in the Cystic Fibrosis Foundation WGS cohort)
# =============================================================================
# Every model is adjusted for sex and ancestry principal components
# (default PC1-PC10; any PCs you supply, e.g. from somalier, PLINK or EIGENSOFT):
#
#   1. Haplogroup dosage            phenotype ~ dosage + PCs + sex
#        dosage = number of haplotypes (0, 1 or 2) assigned to a haplogroup;
#        fit separately for each mucin and haplogroup
#   2. VNTR length                  phenotype ~ length + PCs + sex
#        length = summed length of the two haplotypes, z-scored per mucin,
#        so betas are per 1 SD of summed length
#   3. Binned length dosage, median split
#        phenotype ~ bin dosage + PCs + sex; bins: short / long (below / above
#        the per-mucin median haplotype length)
#   4. Binned length dosage, extremes
#        phenotype ~ bin dosage + PCs + sex; bins: bottom 20% / middle 60% /
#        top 20% of the per-mucin haplotype length distribution (bin1/bin2/bin3)
#        Bin dosage = number of haplotypes (0, 1 or 2) falling in the bin.
#
# Multiple-testing correction (Bonferroni, Benjamini-Hochberg) is applied within
# each model across all tests it ran, never pooled across models.
#
# USAGE
#   Rscript CF_association_testing_allmodels.R --input=my_table.tsv --outdir=results
#   Rscript CF_association_testing_allmodels.R --help
#
# INPUT: ONE tab-delimited table, one row per sample x mucin, with a header.
#   Required columns (default names; rename with the --col-* options):
#     sample              sample identifier
#     mucin               mucin / locus name (e.g. MUC5B)
#     phenotype           continuous phenotype (numeric)
#     sex                 male/female (also accepts M/F or 1/2)
#     vntr_length_hap1    VNTR length of haplotype 1 (numeric)
#     vntr_length_hap2    VNTR length of haplotype 2 (numeric)
#     haplogroup_hap1     haplogroup label of haplotype 1
#     haplogroup_hap2     haplogroup label of haplotype 2
#     PC1 ... PC<n>       ancestry PCs for that sample (repeated on each row)
#   Missing values (NA or empty) are allowed; affected samples drop out of the
#   tests that need them. If the same sample appears on several mucin rows its
#   sex and PCs should be identical.
#   A small example table is in example_CF_input/example_input.tsv.
#
# OPTIONS
#   --input=FILE              (required) input table
#   --outdir=results          output directory
#   --pcs=10                  number of PCs to adjust for (PC1..PC<n>)
#   --min-samples=30          minimum samples per test
#   --length-mucins=LIST      comma-separated mucins to use for models 2-4
#                             (default: all mucins in the table)
#   --dosage-exclude=LIST     comma-separated mucins to leave out of model 1
#                             (default: none)
#   --col-sample=NAME  --col-mucin=NAME  --col-phenotype=NAME  --col-sex=NAME
#   --col-length1=NAME  --col-length2=NAME  --col-hap1=NAME  --col-hap2=NAME
#                             use these columns instead of the default names
#
# NOTES ON BINNED MODELS
#   VNTR lengths are discrete and often heavily tied, so some quantile edges
#   coincide and neighbouring bins collapse. Zero-width bins are dropped and the
#   remaining bins keep their position (short/long for the median split).
#
# OUTPUT (in --outdir; <tag> is e.g. 10pcs)
#   dosage_association_results_<tag>.csv
#   vntr_length_sum_results_<tag>.csv
#   binned_median_split_results_<tag>.csv
#   binned_bottom20_top20_results_<tag>.csv
#   all_models_summary_<tag>.csv
#   run_log_<tag>.txt
#
# REQUIREMENTS: R (>= 4.0) with dplyr and broom.
# The script uses no random numbers, so results are deterministic.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(broom)
})

# =============================================================================
# ARGUMENT PARSING
# =============================================================================
DEFAULTS <- list(
  input          = NA_character_,
  outdir         = "results",
  pcs            = 10L,
  min_samples    = 30L,
  length_mucins  = NA_character_,
  dosage_exclude = NA_character_,
  col_sample     = "sample",
  col_mucin      = "mucin",
  col_phenotype  = "phenotype",
  col_sex        = "sex",
  col_length1    = "vntr_length_hap1",
  col_length2    = "vntr_length_hap2",
  col_hap1       = "haplogroup_hap1",
  col_hap2       = "haplogroup_hap2"
)

USAGE <- paste(
  "Usage: Rscript CF_association_testing_allmodels.R --input=FILE [options]",
  "",
  "Options:  --outdir=results  --pcs=10  --min-samples=30",
  "          --length-mucins=MUC5AC,MUC5B   (models 2-4; default all mucins)",
  "          --dosage-exclude=MUC2,MUC4     (model 1; default none)",
  "Column names (defaults shown):",
  "          --col-sample=sample  --col-mucin=mucin  --col-phenotype=phenotype",
  "          --col-sex=sex  --col-length1=vntr_length_hap1  --col-length2=vntr_length_hap2",
  "          --col-hap1=haplogroup_hap1  --col-hap2=haplogroup_hap2",
  "The table also needs columns PC1..PC<pcs>. See the header of this script",
  "for the input format and details.",
  sep = "\n")

split_list <- function(x) {
  if (is.na(x)) return(NULL)
  v <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
  v <- v[nzchar(v)]
  if (length(v) == 0) NULL else v
}

parse_args <- function(argv, defaults = DEFAULTS) {
  if (any(argv %in% c("-h", "--help"))) { cat(USAGE, "\n"); quit(status = 0) }
  opts <- defaults
  for (a in argv) {
    if (!grepl("^--[A-Za-z0-9-]+=.+$", a)) stop("Unrecognised argument (use --option=value): ", a, call. = FALSE)
    body <- sub("^--", "", a)
    key  <- gsub("-", "_", sub("=.*$", "", body))
    val  <- sub("^[^=]*=", "", body)
    if (!key %in% names(opts)) stop("Unknown option --", gsub("_", "-", key), call. = FALSE)
    opts[[key]] <- val
  }
  if (is.na(opts$input)) stop("Missing required option --input\n\n", USAGE, call. = FALSE)
  opts$pcs <- suppressWarnings(as.integer(opts$pcs))
  opts$min_samples <- suppressWarnings(as.integer(opts$min_samples))
  if (is.na(opts$pcs) || opts$pcs < 1) stop("--pcs must be a positive integer", call. = FALSE)
  if (is.na(opts$min_samples) || opts$min_samples < 1) stop("--min-samples must be a positive integer", call. = FALSE)
  opts
}

pc_names <- function(n) paste0("PC", seq_len(n))

# Restrict a long-format table (one row per sample x mucin) to a set of mucins.
select_mucins <- function(df, include = NULL, exclude = NULL) {
  if (!is.null(include)) df <- df[df$mucin %in% include, , drop = FALSE]
  if (!is.null(exclude)) df <- df[!df$mucin %in% exclude, , drop = FALSE]
  df
}

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================
read_tsv_plain <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE)
}

# "Male"/"male"/"M"/1 -> "male"; "Female"/"female"/"F"/2 -> "female"; anything else NA.
# Sex is spelled inconsistently across tables ("Male" vs "male"), so it is normalised.
normalize_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("1", "m", "male")]   <- "male"
  x[x %in% c("2", "f", "female")] <- "female"
  x[!x %in% c("male", "female")]  <- NA
  x
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

complete_rows <- function(df, cols) df[stats::complete.cases(df[, cols, drop = FALSE]), , drop = FALSE]


# Read the single input table, check/rename its columns to the internal names
# used by the models, and convert types.
load_input <- function(opts, pc_cols) {
  df <- read_tsv_plain(opts$input)
  map <- c(LoO_sample = opts$col_sample, mucin = opts$col_mucin, knorma = opts$col_phenotype,
           sex = opts$col_sex, pred_vntr_length_hap1 = opts$col_length1,
           pred_vntr_length_hap2 = opts$col_length2,
           pred_haplogroup_hap1 = opts$col_hap1, pred_haplogroup_hap2 = opts$col_hap2)
  miss <- c(setdiff(unname(map), names(df)), setdiff(pc_cols, names(df)))
  if (length(miss) > 0) {
    stop("Input table is missing column(s): ", paste(miss, collapse = ", "),
         "\nColumns found: ", paste(names(df), collapse = ", "),
         "\nUse the --col-* options to name your columns, and supply PC1..PC", length(pc_cols), ".",
         call. = FALSE)
  }
  out <- df[, c(unname(map), pc_cols), drop = FALSE]
  names(out)[seq_along(map)] <- names(map)
  out$LoO_sample <- as.character(out$LoO_sample)
  out$knorma <- to_num(out$knorma)
  out$pred_vntr_length_hap1 <- to_num(out$pred_vntr_length_hap1)
  out$pred_vntr_length_hap2 <- to_num(out$pred_vntr_length_hap2)
  out$sex <- normalize_sex(out$sex)
  for (pc in pc_cols) out[[pc]] <- to_num(out[[pc]])
  if (anyDuplicated(out[, c("LoO_sample", "mucin")]) > 0) {
    stop("The table has more than one row for the same sample and mucin; ",
         "it must have exactly one row per sample x mucin.", call. = FALSE)
  }
  out
}

# =============================================================================
# SHARED MODEL HELPERS
# =============================================================================
covariate_terms <- function(df, pc_cols) {
  c(pc_cols, if (nlevels(droplevels(as.factor(df$sex))) > 1) "sex")
}

fit_lm_terms <- function(sub_df, predictors, pc_cols) {
  sub_df$sex <- droplevels(as.factor(sub_df$sex))
  fml   <- reformulate(c(predictors, covariate_terms(sub_df, pc_cols)), response = "knorma")
  model <- lm(fml, data = sub_df)
  res   <- tidy(model) %>% filter(term %in% predictors, !is.na(estimate), !is.na(p.value))
  list(model = model, res = res)
}

add_corrections <- function(res) {
  res %>% mutate(
    p_bonferroni = p.adjust(p.value, method = "bonferroni"),
    p_fdr_bh     = p.adjust(p.value, method = "BH"),
    p_fdr_by     = p.adjust(p.value, method = "BY")
  )
}

split_by_mucin <- function(df) split(df, df$mucin)

# =============================================================================
# MODEL 1: HAPLOGROUP DOSAGE (phenotype ~ dosage + PCs + sex)
# =============================================================================
calculate_dosage <- function(df, haplogroup) {
  h <- as.character(haplogroup)
  as.numeric(as.character(df$pred_haplogroup_hap1) == h) +
    as.numeric(as.character(df$pred_haplogroup_hap2) == h)
}

test_haplogroup_dosage <- function(mdf, mucin_name, haplogroup, pc_cols,
                                   min_samples, min_dosage_count = 5) {
  # Samples with an unassigned haplogroup on either haplotype are excluded
  sub_df <- mdf[!is.na(mdf$pred_haplogroup_hap1) & !is.na(mdf$pred_haplogroup_hap2), , drop = FALSE]
  if (nrow(sub_df) < min_samples) return(NULL)

  sub_df$dosage <- calculate_dosage(sub_df, haplogroup)
  if (length(unique(sub_df$dosage)) < 2 || sum(sub_df$dosage > 0) < min_dosage_count) return(NULL)

  sub_df <- complete_rows(sub_df, c(pc_cols, "sex", "knorma"))
  if (nrow(sub_df) < min_samples) return(NULL)
  if (length(unique(sub_df$dosage)) < 2 || sum(sub_df$dosage > 0) < min_dosage_count) return(NULL)

  tryCatch({
    fit <- fit_lm_terms(sub_df, "dosage", pc_cols)
    if (nrow(fit$res) == 0) return(NULL)
    gl <- glance(fit$model)
    fit$res %>% mutate(
      mucin = mucin_name, haplogroup = as.character(haplogroup), n_pcs = length(pc_cols),
      n_samples = nrow(sub_df),
      n_dosage_0 = sum(sub_df$dosage == 0), n_dosage_1 = sum(sub_df$dosage == 1),
      n_dosage_2 = sum(sub_df$dosage == 2),
      mean_phenotype_dosage0 = mean(sub_df$knorma[sub_df$dosage == 0]),
      mean_phenotype_dosage1 = mean(sub_df$knorma[sub_df$dosage == 1]),
      mean_phenotype_dosage2 = mean(sub_df$knorma[sub_df$dosage == 2]),
      model_r_squared = gl$r.squared, model_adj_r_squared = gl$adj.r.squared
    )
  }, error = function(e) {
    warning("Dosage model failed for ", mucin_name, " haplogroup ", haplogroup, ": ", conditionMessage(e))
    NULL
  })
}

run_haplogroup_dosage <- function(df, pc_cols, min_samples) {
  haplogroups <- unique(c(as.character(df$pred_haplogroup_hap1), as.character(df$pred_haplogroup_hap2)))
  haplogroups <- sort(haplogroups[!is.na(haplogroups) & haplogroups != ""])
  by_mucin <- split_by_mucin(df)
  cat("Mucins:", length(by_mucin), " Haplogroups:", length(haplogroups),
      " Candidate tests:", length(by_mucin) * length(haplogroups), "\n")

  out <- list()
  for (m in names(by_mucin)) for (h in haplogroups) {
    r <- test_haplogroup_dosage(by_mucin[[m]], m, h, pc_cols, min_samples)
    if (!is.null(r)) out[[length(out) + 1]] <- r
  }
  if (length(out) == 0) return(NULL)

  bind_rows(out) %>% add_corrections() %>%
    select(mucin, haplogroup, n_pcs, estimate, std.error, statistic, p.value,
           p_bonferroni, p_fdr_bh, p_fdr_by, n_samples, n_dosage_0, n_dosage_1, n_dosage_2,
           mean_phenotype_dosage0, mean_phenotype_dosage1, mean_phenotype_dosage2,
           model_r_squared, model_adj_r_squared) %>%
    arrange(p.value)
}

# =============================================================================
# MODEL 2: VNTR LENGTH (phenotype ~ length + PCs + sex)
# =============================================================================
test_vntr_length <- function(mdf, mucin_name, length_model, pc_cols, min_samples,
                             scale_length = TRUE) {
  sub_df <- mdf[!is.na(mdf$pred_vntr_length_hap1) & !is.na(mdf$pred_vntr_length_hap2), , drop = FALSE]
  sub_df <- complete_rows(sub_df, c(pc_cols, "sex", "knorma"))
  if (nrow(sub_df) < min_samples) return(NULL)

  h1 <- sub_df$pred_vntr_length_hap1
  h2 <- sub_df$pred_vntr_length_hap2
  lens <- list(
    length_sum = h1 + h2, length_mean = (h1 + h2) / 2,
    length_min = pmin(h1, h2), length_max = pmax(h1, h2),
    length_hap1 = h1, length_hap2 = h2
  )
  for (nm in names(lens)) {
    v <- lens[[nm]]
    sub_df[[nm]] <- if (scale_length) { if (sd(v) > 0) as.numeric(scale(v)) else NA_real_ } else v
  }

  terms <- switch(length_model,
    sum = "length_sum", mean = "length_mean", min = "length_min", max = "length_max",
    both = c("length_hap1", "length_hap2"),
    stop("Invalid length_model: ", length_model, call. = FALSE))
  if (any(is.na(sub_df[[terms[1]]]))) return(NULL)   # no variation in length

  tryCatch({
    fit <- fit_lm_terms(sub_df, terms, pc_cols)
    if (nrow(fit$res) == 0) return(NULL)
    gl <- glance(fit$model)
    all1 <- mdf$pred_vntr_length_hap1; all2 <- mdf$pred_vntr_length_hap2
    fit$res %>% mutate(
      mucin = mucin_name, length_model = length_model, n_pcs = length(pc_cols),
      n_samples = nrow(sub_df), scaled = scale_length,
      mean_length_hap1 = mean(all1, na.rm = TRUE), sd_length_hap1 = sd(all1, na.rm = TRUE),
      mean_length_hap2 = mean(all2, na.rm = TRUE), sd_length_hap2 = sd(all2, na.rm = TRUE),
      length_range = paste0(round(min(c(all1, all2), na.rm = TRUE)), "-",
                            round(max(c(all1, all2), na.rm = TRUE))),
      model_r_squared = gl$r.squared, model_adj_r_squared = gl$adj.r.squared
    )
  }, error = function(e) {
    warning("Length model failed for ", mucin_name, ": ", conditionMessage(e))
    NULL
  })
}

run_vntr_length <- function(df, pc_cols, min_samples, length_model = "sum", scale_length = TRUE) {
  by_mucin <- split_by_mucin(df)
  out <- list()
  for (m in names(by_mucin)) {
    r <- test_vntr_length(by_mucin[[m]], m, length_model, pc_cols, min_samples, scale_length)
    if (!is.null(r)) out[[length(out) + 1]] <- r
  }
  if (length(out) == 0) return(NULL)

  bind_rows(out) %>%
    mutate(p_bonferroni = p.adjust(p.value, method = "bonferroni"),
           p_fdr_bh     = p.adjust(p.value, method = "BH")) %>%
    select(mucin, term, length_model, n_pcs, estimate, std.error, statistic, p.value,
           p_bonferroni, p_fdr_bh, n_samples, length_range,
           mean_length_hap1, sd_length_hap1, mean_length_hap2, sd_length_hap2,
           model_r_squared, model_adj_r_squared, scaled) %>%
    arrange(p.value)
}

# =============================================================================
# MODEL 3: BINNED VNTR-LENGTH DOSAGE
# =============================================================================
# When quantile edges coincide (discrete, tied lengths) zero-width bins are
# dropped and the remaining bins are relabelled by order.
collapse_bin_labels <- function(method, labels, keep) {
  nb <- sum(keep)
  if (method %in% c("tertiles", "quartiles", "median")) {
    return(switch(as.character(nb),
                  "1" = "short",
                  "2" = c("short", "long"),
                  "3" = c("short", "medium", "long"),
                  labels[seq_len(nb)]))
  }
  labels[keep]
}

get_length_bins <- function(lengths, method = "tertiles", custom_probs = NULL) {
  lengths <- lengths[!is.na(lengths)]
  if (length(lengths) == 0) return(NULL)

  if (method == "tertiles") {
    breaks <- quantile(lengths, probs = c(0, 1/3, 2/3, 1)); labels <- c("short", "medium", "long")
  } else if (method == "quartiles") {
    breaks <- quantile(lengths, probs = c(0, 0.25, 0.5, 0.75, 1)); labels <- c("short", "med_short", "med_long", "long")
  } else if (method == "median") {
    breaks <- quantile(lengths, probs = c(0, 0.5, 1)); labels <- c("short", "long")
  } else if (method == "custom" && !is.null(custom_probs)) {
    breaks <- quantile(lengths, probs = custom_probs); labels <- paste0("bin", seq_len(length(custom_probs) - 1))
  } else stop("Invalid bin method: ", method, call. = FALSE)

  breaks <- unname(breaks)
  if (any(duplicated(breaks))) {          # ties (discrete VNTR lengths): collapse bins
    keep   <- diff(breaks) > 0
    breaks <- unique(breaks)
    labels <- collapse_bin_labels(method, labels, keep)
  }
  if (length(breaks) < 2) return(NULL)    # all lengths identical, nothing to bin
  list(breaks = breaks, labels = labels)
}

assign_length_bins <- function(df, bin_info) {
  df$length_bin_hap1 <- cut(df$pred_vntr_length_hap1, breaks = bin_info$breaks,
                            labels = bin_info$labels, include.lowest = TRUE)
  df$length_bin_hap2 <- cut(df$pred_vntr_length_hap2, breaks = bin_info$breaks,
                            labels = bin_info$labels, include.lowest = TRUE)
  df
}

test_bin_dosage <- function(mdf, mucin_name, bin_category, bin_info, pc_cols, min_samples,
                            min_dosage_count = 5) {
  sub_df <- mdf[!is.na(mdf$pred_vntr_length_hap1) & !is.na(mdf$pred_vntr_length_hap2), , drop = FALSE]
  if (nrow(sub_df) == 0) return(NULL)
  sub_df <- assign_length_bins(sub_df, bin_info)
  sub_df$dosage <- as.numeric(sub_df$length_bin_hap1 == bin_category) +
                   as.numeric(sub_df$length_bin_hap2 == bin_category)

  sub_df <- complete_rows(sub_df, c(pc_cols, "sex", "knorma", "dosage"))
  if (nrow(sub_df) < min_samples) return(NULL)
  if (length(unique(sub_df$dosage)) < 2 || sum(sub_df$dosage > 0) < min_dosage_count) return(NULL)

  tryCatch({
    fit <- fit_lm_terms(sub_df, "dosage", pc_cols)
    if (nrow(fit$res) == 0) return(NULL)
    gl <- glance(fit$model)
    idx <- which(bin_info$labels == bin_category)
    fit$res %>% mutate(
      mucin = mucin_name, length_bin = as.character(bin_category), n_pcs = length(pc_cols),
      n_samples = nrow(sub_df),
      n_dosage_0 = sum(sub_df$dosage == 0), n_dosage_1 = sum(sub_df$dosage == 1),
      n_dosage_2 = sum(sub_df$dosage == 2),
      mean_pheno_dose0 = mean(sub_df$knorma[sub_df$dosage == 0]),
      mean_pheno_dose1 = mean(sub_df$knorma[sub_df$dosage == 1]),
      mean_pheno_dose2 = mean(sub_df$knorma[sub_df$dosage == 2]),
      bin_range_min = bin_info$breaks[idx], bin_range_max = bin_info$breaks[idx + 1],
      model_r_squared = gl$r.squared
    )
  }, error = function(e) {
    warning("Binned model failed for ", mucin_name, " ", bin_category, ": ", conditionMessage(e))
    NULL
  })
}

run_binned_dosage <- function(df, pc_cols, min_samples, bin_method = "tertiles",
                              custom_probs = NULL) {
  by_mucin <- split_by_mucin(df)
  out <- list()
  for (m in names(by_mucin)) {
    mdf <- by_mucin[[m]]
    # Bin edges come from every haplotype length for the mucin, independent of covariates
    bin_info <- get_length_bins(c(mdf$pred_vntr_length_hap1, mdf$pred_vntr_length_hap2),
                                method = bin_method, custom_probs = custom_probs)
    if (is.null(bin_info)) next
    for (b in bin_info$labels) {
      r <- test_bin_dosage(mdf, m, b, bin_info, pc_cols, min_samples)
      if (!is.null(r)) out[[length(out) + 1]] <- r
    }
  }
  if (length(out) == 0) return(NULL)

  bind_rows(out) %>%
    mutate(p_bonferroni = p.adjust(p.value, method = "bonferroni"),
           p_fdr_bh     = p.adjust(p.value, method = "BH")) %>%
    select(mucin, length_bin, n_pcs, bin_range_min, bin_range_max,
           estimate, std.error, statistic, p.value, p_bonferroni, p_fdr_bh,
           n_samples, n_dosage_0, n_dosage_1, n_dosage_2,
           mean_pheno_dose0, mean_pheno_dose1, mean_pheno_dose2, model_r_squared) %>%
    arrange(p.value)
}

# =============================================================================
# REPORTING
# =============================================================================
describe_results <- function(res, title, label_fn, effect_fn, alpha = 0.05) {
  cat("\n=============================================================\n")
  cat(" ", title, "\n")
  cat("=============================================================\n")
  if (is.null(res) || nrow(res) == 0) { cat("No tests completed.\n"); return(invisible(NULL)) }

  cat(sprintf("Tests: %d | nominal p<%.2f: %d | Bonferroni<%.2f: %d | FDR(BH)<%.2f: %d | FDR(BH)<0.10: %d\n",
              nrow(res), alpha, sum(res$p.value < alpha), alpha, sum(res$p_bonferroni < alpha),
              alpha, sum(res$p_fdr_bh < alpha), sum(res$p_fdr_bh < 0.10)))
  cat("Top associations:\n")
  top <- head(res, 10)
  for (i in seq_len(nrow(top))) {
    r <- top[i, ]
    cat(sprintf("  %2d. %-28s %-34s p=%.2e  FDR=%.2e\n", i, label_fn(r), effect_fn(r), r$p.value, r$p_fdr_bh))
  }
}

summary_row <- function(res, analysis, label_fn) {
  if (is.null(res) || nrow(res) == 0) {
    return(data.frame(analysis = analysis, n_tests = 0L, n_nominal = NA, n_bonferroni = NA,
                      n_fdr_0.05 = NA, n_fdr_0.10 = NA, min_p = NA, top_hit = NA, top_hit_fdr = NA,
                      stringsAsFactors = FALSE))
  }
  data.frame(analysis = analysis, n_tests = nrow(res),
             n_nominal = sum(res$p.value < 0.05), n_bonferroni = sum(res$p_bonferroni < 0.05),
             n_fdr_0.05 = sum(res$p_fdr_bh < 0.05), n_fdr_0.10 = sum(res$p_fdr_bh < 0.10),
             min_p = min(res$p.value), top_hit = label_fn(res[1, ]), top_hit_fdr = res$p_fdr_bh[1],
             stringsAsFactors = FALSE)
}

lab_hap  <- function(r) paste0(r$mucin, " hap ", r$haplogroup)
lab_len  <- function(r) paste0(r$mucin, " ", r$term)
lab_bin  <- function(r) sprintf("%s %s (%d-%d bp)", r$mucin, r$length_bin, round(r$bin_range_min), round(r$bin_range_max))
eff_beta <- function(r) sprintf("beta=%.4f (SE %.4f)", r$estimate, r$std.error)

save_csv <- function(res, outdir, name, tag) {
  if (is.null(res)) return(invisible(NULL))
  path <- file.path(outdir, sprintf("%s_%s.csv", name, tag))
  write.csv(res, path, row.names = FALSE)
  cat("  wrote", path, "\n")
}

log_environment <- function(argv) {
  pk <- c("dplyr", "broom")
  ver <- vapply(pk, function(p) {
    if (requireNamespace(p, quietly = TRUE)) paste0(p, " ", as.character(utils::packageVersion(p)))
    else paste0(p, " (not installed)")
  }, character(1))
  cat(R.version.string, "\n")
  cat("Packages:", paste(ver, collapse = ", "), "\n")
  cat("Arguments:", if (length(argv)) paste(argv, collapse = " ") else "(defaults)", "\n")
}

# =============================================================================
# MAIN
# =============================================================================
main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opts    <- parse_args(argv)
  n_pcs   <- opts$pcs
  pc_cols <- pc_names(n_pcs)
  tag     <- sprintf("%dpcs", n_pcs)

  dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)
  log_con <- file(file.path(opts$outdir, sprintf("run_log_%s.txt", tag)), open = "wt")
  sink(log_con, split = TRUE)
  on.exit({ sink(); close(log_con) }, add = TRUE)

  cat("Association models with", n_pcs, "PCs + sex  (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")\n")
  log_environment(argv)
  cat("Output directory:", opts$outdir, "\n\n")

  # ---- load the input table --------------------------------------------------
  cff <- load_input(opts, pc_cols)
  cat(sprintf("Input: %s\nrows: %d, samples: %d, mucins: %d, sex: %s\n\n", opts$input, nrow(cff),
              length(unique(cff$LoO_sample)), length(unique(cff$mucin)),
              paste(names(table(cff$sex, useNA = "ifany")), table(cff$sex, useNA = "ifany"), sep = "=", collapse = ", ")))

  # ---- mucin sets per model (default: all mucins in the table) ----------------
  len_include <- split_list(opts$length_mucins)
  dos_exclude <- split_list(opts$dosage_exclude)
  cff_len <- select_mucins(cff, include = len_include)   # models 2-4
  cff_dos <- select_mucins(cff, exclude = dos_exclude)   # model 1
  if (nrow(cff_len) == 0 || nrow(cff_dos) == 0) stop("No rows left after applying the mucin selection; check --length-mucins / --dosage-exclude.", call. = FALSE)
  cat("Mucins in models 2-4:", paste(sort(unique(cff_len$mucin)), collapse = ", "), "\n")
  cat("Mucins in model 1:", paste(sort(unique(cff_dos$mucin)), collapse = ", "), "\n\n")

  summaries <- list()

  # ---- 1. haplogroup dosage --------------------------------------------------
  cat("Running haplogroup dosage ...\n")
  res_dosage <- run_haplogroup_dosage(cff_dos, pc_cols, opts$min_samples)
  describe_results(res_dosage, sprintf("1. HAPLOGROUP DOSAGE (%d PCs + sex)", n_pcs), lab_hap, eff_beta)
  save_csv(res_dosage, opts$outdir, "dosage_association_results", tag)
  summaries[["haplogroup_dosage"]] <- summary_row(res_dosage, "haplogroup_dosage", lab_hap)

  # ---- 2. VNTR length --------------------------------------------------------
  for (lm_name in "sum") {   # summed haplotype length (the reported length model)
    cat("\nRunning VNTR length model:", lm_name, "...\n")
    res <- run_vntr_length(cff_len, pc_cols, opts$min_samples, length_model = lm_name)
    describe_results(res, sprintf("2. VNTR LENGTH [%s] (%d PCs + sex; beta per 1 SD)", lm_name, n_pcs), lab_len, eff_beta)
    save_csv(res, opts$outdir, paste0("vntr_length_", lm_name, "_results"), tag)
    summaries[[paste0("vntr_length_", lm_name)]] <- summary_row(res, paste0("vntr_length_", lm_name), lab_len)
  }

  # ---- 3. binned length dosage: median split ----------------------------------
  cat("\nRunning binned length dosage, median split ...\n")
  res_med <- run_binned_dosage(cff_len, pc_cols, opts$min_samples, "median")
  describe_results(res_med, sprintf("3. BINNED MEDIAN SPLIT (%d PCs + sex)", n_pcs), lab_bin, eff_beta)
  save_csv(res_med, opts$outdir, "binned_median_split_results", tag)
  summaries[["binned_median_split"]] <- summary_row(res_med, "binned_median_split", lab_bin)

  # ---- 4. binned length dosage: extremes (bottom 20% / middle 60% / top 20%) ---
  cat("\nRunning binned length dosage, bottom 20% / top 20% ...\n")
  res_ext <- run_binned_dosage(cff_len, pc_cols, opts$min_samples, "custom",
                               custom_probs = c(0, 0.2, 0.8, 1))
  describe_results(res_ext, sprintf("4. BINNED BOTTOM20 / TOP20 (%d PCs + sex)", n_pcs), lab_bin, eff_beta)
  save_csv(res_ext, opts$outdir, "binned_bottom20_top20_results", tag)
  summaries[["binned_bottom20_top20"]] <- summary_row(res_ext, "binned_bottom20_top20", lab_bin)

  # ---- overall summary -------------------------------------------------------
  overall <- bind_rows(summaries)
  cat("\n=============================================================\n")
  cat(" SUMMARY ACROSS ALL MODELS (", n_pcs, " PCs + sex)\n", sep = "")
  cat("=============================================================\n")
  print(as.data.frame(overall), row.names = FALSE, digits = 3)
  save_csv(overall, opts$outdir, "all_models_summary", tag)

  cat("\nINTERPRETATION\n")
  cat("  Linear models: positive beta = higher phenotype with more copies / longer VNTR.\n")
  cat("  VNTR length betas are per 1 SD of length. Binned betas are per copy of that bin.\n")
  cat("  FDR/Bonferroni are applied across the included mucin loci within each analysis, not across analyses.\n")
  cat("\nAll outputs in:", opts$outdir, "\n")
  invisible(overall)
}

if (!interactive()) main()

