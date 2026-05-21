## Mucin VNTR motif finding in protein sequences

library(Biostrings)

analyze_repeats <- function(protein_seq, canonical_repeat, split_blocks = FALSE, max_gap = NULL) {
  # Guard against NA, NULL, empty, or non-character input
  if (is.null(protein_seq) || length(protein_seq) == 0 ||
      is.na(protein_seq) || !is.character(protein_seq) ||
      nchar(trimws(protein_seq)) == 0) {
    return(NULL)
  }
  if (is.null(canonical_repeat) || length(canonical_repeat) == 0 ||
      is.na(canonical_repeat) || !is.character(canonical_repeat) ||
      nchar(trimws(canonical_repeat)) == 0) {
    return(NULL)
  }

  # Strip any whitespace or asterisks (stop codons) from protein sequence
  protein_seq <- gsub("[\\s*]", "", protein_seq, perl = TRUE)
  if (nchar(protein_seq) == 0) return(NULL)

  repeat_len <- nchar(canonical_repeat)

  if (is.null(max_gap)) {
    max_gap <- repeat_len * 2
  }

  # --- Stage 1: Find all individual repeat matches ---
  doubled <- paste0(canonical_repeat, canonical_repeat)

  aln_double <- matchPattern(
    AAString(doubled),
    AAString(protein_seq),
    max.mismatch = floor(repeat_len * 0.8),
    with.indels = FALSE
  )

  aln_single <- matchPattern(
    AAString(canonical_repeat),
    AAString(protein_seq),
    max.mismatch = floor(repeat_len * 0.4),
    with.indels = FALSE
  )

  all_starts <- c(start(aln_double), start(aln_single))
  all_ends <- c(end(aln_double), end(aln_single))

  if (length(all_starts) == 0) {
    return(NULL)
  }

  ord <- order(all_starts)
  all_starts <- all_starts[ord]
  all_ends <- all_ends[ord]

  # --- Stage 2: Define repeat blocks ---
  if (split_blocks) {
    # Cluster matches into separate blocks based on max_gap
    blocks <- list()
    block_start <- all_starts[1]
    block_end <- all_ends[1]

    for (j in seq_along(all_starts)) {
      if (j == 1) {
        next
      }
      gap <- all_starts[j] - block_end

      if (gap > max_gap) {
        blocks[[length(blocks) + 1]] <- c(block_start, block_end)
        block_start <- all_starts[j]
        block_end <- all_ends[j]
      } else {
        block_end <- max(block_end, all_ends[j])
      }
    }
    blocks[[length(blocks) + 1]] <- c(block_start, block_end)
  } else {
    # Single block spanning from first match to last match
    blocks <- list(c(min(all_starts), max(all_ends)))
  }

  # --- Stage 3: For each block, determine phase and tile ---
  all_results <- list()
  global_repeat_num <- 1

  for (b in seq_along(blocks)) {
    region_start <- blocks[[b]][1]
    region_end <- blocks[[b]][2]
    repeat_region <- substr(protein_seq, region_start, region_end)

    # Phase detection
    best_offset <- 0
    best_score <- -Inf

    for (offset in 0:(repeat_len - 1)) {
      rotated <- paste0(
        substr(canonical_repeat, offset + 1, repeat_len),
        substr(canonical_repeat, 1, offset)
      )
      compare_len <- min(repeat_len, nchar(repeat_region))
      rotated_aa <- strsplit(substr(rotated, 1, compare_len), "")[[1]]
      region_aa <- strsplit(substr(repeat_region, 1, compare_len), "")[[1]]

      matches <- sum(rotated_aa == region_aa)
      if (matches > best_score) {
        best_score <- matches
        best_offset <- offset
      }
    }

    phased_canonical <- paste0(
      substr(canonical_repeat, best_offset + 1, repeat_len),
      substr(canonical_repeat, 1, best_offset)
    )
    canon_aa <- strsplit(phased_canonical, "")[[1]]

    region_len <- nchar(repeat_region)
    n_full <- floor(region_len / repeat_len)

    # Tile full repeats
    for (i in seq_len(n_full)) {
      start_pos <- (i - 1) * repeat_len + 1
      end_pos <- i * repeat_len
      chunk <- substr(repeat_region, start_pos, end_pos)
      chunk_aa <- strsplit(chunk, "")[[1]]
      diffs <- canon_aa != chunk_aa

      all_results[[length(all_results) + 1]] <- data.frame(
        block = b,
        repeat_num = global_repeat_num,
        abs_start = region_start + start_pos - 1,
        abs_end = region_start + end_pos - 1,
        mismatches = sum(diffs),
        pct_identity = round(100 * sum(!diffs) / repeat_len, 1),
        pct_variability = round(100 * sum(diffs) / repeat_len, 1),
        phase_offset = best_offset,
        repeat_seq = chunk,
        stringsAsFactors = FALSE
      )
      global_repeat_num <- global_repeat_num + 1
    }

    # Handle partial repeat at end of block
    leftover <- region_len - (n_full * repeat_len)
    if (leftover > 0) {
      start_pos <- n_full * repeat_len + 1
      chunk <- substr(repeat_region, start_pos, region_len)
      chunk_aa <- strsplit(chunk, "")[[1]]
      partial_canon <- canon_aa[seq_len(leftover)]
      diffs <- partial_canon != chunk_aa

      all_results[[length(all_results) + 1]] <- data.frame(
        block = b,
        repeat_num = global_repeat_num,
        abs_start = region_start + start_pos - 1,
        abs_end = region_end,
        mismatches = sum(diffs),
        pct_identity = round(100 * sum(!diffs) / leftover, 1),
        pct_variability = round(100 * sum(diffs) / leftover, 1),
        phase_offset = best_offset,
        repeat_seq = chunk,
        stringsAsFactors = FALSE
      )
      global_repeat_num <- global_repeat_num + 1
    }
  }

  if (length(all_results) == 0) {
    return(NULL)
  }
  do.call(rbind, all_results)
}

# --- Apply across your dataframe ---
# canonical <- "AHGVTSAPDTRPAPGSTAPP"  # your canonical repeat
#
# all_results <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
#   res <- analyze_repeats(df$protein_seq[i], canonical, split_blocks = FALSE)
#   if (!is.null(res)) {
#     res$haplotype <- df$haplotype[i]
#   }
#   res
# }))



