###function to infer haplogroups from a phylogenetic tree with annotations for VNTR length


infer_haplogroups <- function(tree, meta, vntr_col = "vntr_length",
                              min_cluster_size = 3,
                              bootstrap_cutoff = 70) {

  #--------------------------------
  # 1. Extract bootstrap per node
  #--------------------------------
  gt <- ggtree(tree)
  node_support <- gt$data %>% select(node, label)

  # Internal nodes only
  internal_support <- node_support %>%
    filter(node > length(tree$tip.label)) %>%
    mutate(label = suppressWarnings(as.numeric(label))) %>%
    mutate(label = ifelse(is.na(label), 0, label))
  
  #--------------------------------
  # 2. Compute patristic distances
  #--------------------------------
  pd <- cophenetic.phylo(tree)
  hc <- hclust(as.dist(pd), method = "average")
  tree_height <- max(hc$height)

  thresholds <- seq(tree_height * 0.01, tree_height * 0.9, length.out = 200)

  results <- list()

  #--------------------------------
  # Helper: check if a cluster is supported by bootstrap >= 70
  #--------------------------------
  cluster_supported <- function(tips_in_cluster) {
    mrca <- getMRCA(tree, tips_in_cluster)
    if (is.null(mrca)) return(FALSE)
    support <- internal_support %>% filter(node == mrca) %>% pull(label)
    if (length(support) == 0) support <- 0
    return(support >= bootstrap_cutoff)
  }

  #--------------------------------
  # 3. Loop through thresholds, evaluate R² (VNTR ~ cluster)
  #--------------------------------
  for (h in thresholds) {
    cl <- cutree(hc, h = h)
    df <- data.frame(tip = names(cl), cluster = cl) %>%
            left_join(meta, by = "tip")

    # Remove tiny clusters (< min_cluster_size)
    small <- names(which(table(df$cluster) < min_cluster_size))
    df <- df %>% filter(!cluster %in% small)

    if (n_distinct(df$cluster) < 2) next

    # Apply bootstrap filter
    valid <- logical(nrow(df))
    for (k in unique(df$cluster)) {
      tips_k <- df$tip[df$cluster == k]
      ok <- cluster_supported(tips_k)
      valid[df$cluster == k] <- ok
    }
    df <- df[valid, ]
    if (n_distinct(df$cluster) < 2) next

    # Compute ANOVA / R²
    a <- aov(df[[vntr_col]] ~ as.factor(df$cluster))
    ss <- summary(a)[[1]]
    results[[as.character(h)]] <- list(
      height = h,
      clusters = n_distinct(df$cluster),
      R2 = ss$`Sum Sq`[1] / sum(ss$`Sum Sq`),
      F  = ss$`F value`[1],
      p  = ss$`Pr(>F)`[1]
    )
  }

  if (length(results) == 0) stop("No bootstrap-supported clusterings found.")
  res_df <- do.call(rbind, lapply(results, as.data.frame))
  best <- res_df[which.max(res_df$R2), ]
  best_height <- best$height

  #----------------------------------------
  # Recover final cluster assignments
  #----------------------------------------
  cl_final <- cutree(hc, h = best_height)
  df_final <- data.frame(tip = names(cl_final), cluster = cl_final) %>%
    left_join(meta, by = "tip")

  small2 <- names(which(table(df_final$cluster) < min_cluster_size))
  df_final <- df_final %>% filter(!cluster %in% small2)

  valid2 <- logical(nrow(df_final))
  for (k in unique(df_final$cluster)) {
    tips_k <- df_final$tip[df_final$cluster == k]
    valid2[df_final$cluster == k] <- cluster_supported(tips_k)
  }
  df_final <- df_final[valid2, ]

  #----------------------------------------
  # Assign orphan tips (NA tips) to nearest surviving tip's cluster safely
  #----------------------------------------
  rownames(pd) <- tree$tip.label
  colnames(pd) <- tree$tip.label

  all_tips <- meta$tip
  clustered_tips <- df_final$tip
  orphan_tips <- setdiff(all_tips, clustered_tips)

  if (length(orphan_tips) > 0) {
    for (tip in orphan_tips) {
      if (!(tip %in% rownames(pd))) {
        warning(sprintf("Tip %s not found in distance matrix; assigning NA cluster.", tip))
        next
      }
      # use only clustered_tips that exist in pd
      valid_clustered <- clustered_tips[clustered_tips %in% colnames(pd)]
      if (length(valid_clustered) == 0) {
        warning("No tips left to assign orphans to. Assigning NA cluster.")
        next
      }
      dists <- pd[tip, valid_clustered, drop = TRUE]
      nearest_tip <- valid_clustered[which.min(dists)]
      new_row <- data.frame(
        tip = tip,
        cluster = df_final$cluster[df_final$tip == nearest_tip],
        vntr_length = meta[[vntr_col]][meta$tip == tip]
      )
      df_final <- rbind(df_final, new_row)
    }
  }

  #----------------------------------------
  # Recode clusters to sequential numbers starting from 1
  #----------------------------------------
  df_final$cluster <- as.integer(factor(df_final$cluster, 
                                        levels = sort(unique(df_final$cluster))))

  #----------------------------------------
  # Visualization
  #----------------------------------------
  pdf(paste0("./", curr_muc, "_haplogroup_predictions.pdf"))
  p <- ggtree(tree) %<+% df_final +
       geom_tippoint(aes(color = factor(cluster)), size = 2) +
       theme_tree2() +
       ggtitle(sprintf("Bootstrap >= %d | Optimal height = %.3f", bootstrap_cutoff, best_height))
  plot(p)
  dev.off()

  return(list(
    best_threshold = best,
    clusters = df_final,
    plot = p,
    results = res_df
  ))
}

