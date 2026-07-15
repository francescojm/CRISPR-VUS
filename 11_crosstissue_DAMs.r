# -------- ADDITIONAL RESULTS EXPLORATION, CROSS-TISSUE DAM ANALYSIS -------- #

library(ggplot2)
library(ggvenn)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(stringr)

# Load per-tissue DAM analysis results from RData files
load_results <- function(home, resultPath) {
  setwd(resultPath)

  tissues <- gsub("_results.RData", "", list.files(pattern="results.RData"))

  results <- list()
  ind_iter <- 0
  for (ctiss in tissues) {
    ind_iter <- ind_iter + 1
    load(paste(ctiss, "_results.RData", sep = ""))
    results[[ind_iter]] <- RESTOT
  }
  names(results) <- tissues

  setwd(home)

  return(list(tissues = tissues, results = results))
}

# Plot diagnostic distributions of p-values, rank ratios and medFitEff
# ADJUST IF USED! DEP_THRESHOLD, RR_TH, PRODUCE_PLOTS ...
# make_diagnostic_distributions <- function(home, resultPath, results, tissues) {

#   # p-distribution 
#   p <- c()
#   for (ctiss in tissues) {
#     p <- c(p, results[[ctiss]]$pval_rand[results[[ctiss]]$rank_ratio < 1.6 & results[[ctiss]]$medFitEff < -0.5])
#   }
#   df_p <- data.frame(pvalue = p)
#   pdf(file.path(home, resultPath, "_figures_source", "p_distribution.pdf"), 6, 4)
#   print(
#     ggplot(df_p, aes(x = pvalue)) +
#       geom_histogram(bins = 100) +
#       theme_classic() +
#       geom_vline(xintercept = 0.2, linetype = "dashed", color = "red", size = 0.5)
#   )
#   dev.off()

#   # rank ratio distribution
#   rr <- c()
#   for (ctiss in tissues) {
#     rr <- c(rr, results[[ctiss]]$rank_ratio)
#   }
#   df_rr <- data.frame(RankRatio = rr)
#   pdf(file.path(home, resultPath, "_figures_source", "rr_distribution.pdf"), 6, 4)
#   print(
#     ggplot(df_rr, aes(x = RankRatio)) +
#       geom_histogram(bins = 100) +
#       theme_classic() +
#       geom_vline(xintercept = 1.6, linetype = "dashed", color = "red", size = 0.5)
#   )
#   dev.off()

#   # medFitEff distribution
#   mfe <- c()
#   for (ctiss in tissues) {
#     mfe <- c(mfe, results[[ctiss]]$medFitEff)
#   }
#   df_mfe <- data.frame(medFitEff = mfe)
#   pdf(file.path(home, resultPath, "_figures_source", "medFitEff_distribution.pdf"), 6, 4)
#   print(
#     ggplot(df_mfe, aes(x = medFitEff)) +
#       geom_histogram(bins = 100) +
#       theme_classic() +
#       geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", size = 0.5)
#   )
#   dev.off()

#   invisible(TRUE)
# }

# Compute mutational burden per tissue using annotations and variant data
compute_mut_burden <- function(tissues, CMP_annot, cl_variants) {

  # Average annotated mutational burden per tissue
  mut_burden_anno <- c()
  for (ctiss in tissues) {
    mut_burden_anno <- c(mut_burden_anno, mean(CMP_annot$mutational_burden[CMP_annot$cancer_type == ctiss], na.rm = TRUE))
  }
  names(mut_burden_anno) <- tissues

  # Recompute mutational burden as the mean number of mutated genes per model
  mut_burden_new <- c()
  for (ctiss in tissues) {

    # Select variants belonging to models of the current tissue
    cl_variants_tmp <- cl_variants[cl_variants$model_id %in% CMP_annot$model_id[CMP_annot$cancer_type == ctiss], ]
    num_mut <- c()

    # Count unique mutated genes in each model of the current tissue
    for (model in unique(cl_variants_tmp$model_id)) {
      num_mut <- c(num_mut, length(unique(cl_variants_tmp$gene_symbol_2023[cl_variants_tmp$model_id == model])))
    }

    # Average number of mutated genes across models in the tissue
    mut_burden_new <- c(mut_burden_new, mean(num_mut, na.rm = TRUE))
  }
  names(mut_burden_new) <- tissues

  return(list(mut_burden_anno = mut_burden_anno,
              mut_burden_new  = mut_burden_new))
}

# Count the number of unique DAM-bearing genes in each tissue.
compute_hits_by_tissue <- function(results, tissues, dep_threshold, RR_th) {
  num_hits <- c()
  for (ctiss in tissues) {
    num_hits <- c(num_hits, length(unique(
      results[[ctiss]]$GENE[
        results[[ctiss]]$rank_ratio < RR_th &
          results[[ctiss]]$medFitEff < dep_threshold &
          results[[ctiss]]$hypTest_p < 0.2 &
          results[[ctiss]]$empPval < 0.2
      ]
    )))
  }
  names(num_hits) <- tissues
  return(num_hits)
}

# Report the tissues with the highest and lowest normalized hit counts (after correcting for mutational burden)
report_hits_by_tissue <- function(num_hits_norm) {

  max_tiss <- names(num_hits_norm)[which.max(num_hits_norm)]
  min_tiss <- names(num_hits_norm)[which.min(num_hits_norm)]

  cat("the tissue with the highest percentage of DAM-bearing genes (considering mutational burden) is ",
      max_tiss, "with", num_hits_norm[max_tiss], "DAM-bearing genes\n", sep = " ")

  cat("the tissue with the lowest percentage DAM-bearing genes (considering mutational burden) is ",
      min_tiss, "with", num_hits_norm[min_tiss], "DAM-bearing genes\n", sep = " ")

  invisible(list(max = max_tiss, min = min_tiss))
}

# Collect all DAM-bearing genes across tissues
collect_hits <- function(results, tissues, dep_threshold, RR_th) {
  hits <- c()
  for (ctiss in tissues) {
    hits <- c(hits, unique(
      results[[ctiss]]$GENE[
        results[[ctiss]]$rank_ratio < RR_th &
          results[[ctiss]]$medFitEff < dep_threshold &
          results[[ctiss]]$hypTest_p < 0.2 &
          results[[ctiss]]$empPval < 0.2
      ]
    ))
  }
  return(hits)
}

# Print summary statistics for DAM-bearing genes across tissues,
# including recurrence across tissues and overlap with known drivers.
report_hit_summary <- function(hits, driver_genes) {

  cat("Number of unique DAM-bearing genes:", length(unique(hits)), "\n")

  hits_table <- table(hits)
  # Keep genes observed in more than one cancer type
  hits_2ormore <- sort(hits_table, decreasing = TRUE)[1:sum(hits_table>1)]

  cat("Number of DAM-bearing genes in at least two cancer types:", length(hits_2ormore), "\n")
  cat("of which not known to be drivers:", length(setdiff(names(hits_2ormore), driver_genes)), "\n")
  cat("Number of DAM-bearing genes known as driver:", length(intersect(hits, driver_genes)), "\n")
  cat("Number of DAM-bearing genes not known as driver:", length(setdiff(hits, driver_genes)), "\n")

  invisible(TRUE)
}

# Plot how many tissues each hit gene appears in.
# Produces one plot for all hits and one restricted to non-driver hits.
plot_hit_frequency <- function(hits, driver_genes, figuresPath, produce_plots = TRUE) {
  
  # Distribution of the number of tissues in which each hit gene is detected
  df <- data.frame(counts = table(table(hits)))
  df$counts.Freq <- log10(df$counts.Freq + 1)
  df$counts.Var1 <- as.numeric(as.character(df$counts.Var1))
  
  if (produce_plots) {
    # not assigned figure number 
    pdf(file.path(figuresPath, "num_tissues_perhit_log.pdf"), 5, 5)
    print(
      ggplot(df, aes(x = counts.Var1, y = counts.Freq)) +
        geom_bar(stat = "identity") +
        theme_classic() +
        ylab("log10(counts+1)") +
        xlab("Number of cancer types") +
        scale_x_continuous(breaks = c(1:15))
    )
    dev.off()
  }
  
  # Same distribution restricted to genes not annotated as drivers
  df2 <- data.frame(counts = table(table(hits)[setdiff(hits, driver_genes)]))
  df2$counts.Freq <- log10(df2$counts.Freq + 1)
  df2$counts.Var1 <- as.numeric(as.character(df2$counts.Var1))

  if (produce_plots) {
    # not assigned figure number 
    pdf(file.path(figuresPath, "num_tissues_perhit_nondrivers_log.pdf"), 2, 5)
    print(
      ggplot(df2, aes(x = counts.Var1, y = counts.Freq)) +
        geom_bar(stat = "identity") +
        theme_classic() +
        ylab("log10(counts+1)") +
        xlab("Number of cancer types") +
        scale_x_continuous(breaks = c(1:10))
    )
    dev.off()
  }

  invisible(TRUE)
}

# Plot overlap between detected hits and known driver genes.
# Returns the fraction of hits that are known drivers.
plot_driver_hits_venn <- function(figuresPath, driver_genes, hits, produce_plots = TRUE) {
  
  set.seed(245374)

  # FIGURE 3A
  if (produce_plots) {
    pdf(file.path(figuresPath, "Venn_driver_hits.pdf"), 5, 5)
    print(ggvenn(
      data = list("Driver genes" = driver_genes, "Hits" = hits),
      text_size = 5, fill_color = c("#F0E442", "#0072B2"), show_percentage = FALSE
    ))
    dev.off()
  }

  perc_int <- length(intersect(hits, driver_genes)) / length(hits)
  cat("Intersection ", perc_int, "\n", sep = "")
  return(perc_int)
}

# Compare detected hits with several external benchmark driver lists.
# For each list:
# - test enrichment with Fisher's exact test
# - save a Venn diagram showing overlap with hits
plot_other_driver_venns_and_tests <- function(benchmark, hits, cl_variants, figuresPath, produce_plots = TRUE) {

  for (col in 1:8) {
    dg1 <- na.omit(benchmark[-1, col])

    # Counts used for Fisher enrichment test
    a <- length(intersect(hits, dg1))
    b <- length(setdiff(hits, dg1))
    c <- length(setdiff(dg1, hits))
    # Background gene universe
    d <- length(unique(cl_variants$gene_symbol))  # or gene_symbol_2023 ??

    cat(as.character(benchmark[1, col]), " ",
        fisher.test(matrix(c(a, b, c, d), ncol = 2), alternative = "greater")$p.value,
        "\n", sep = "")
    
    # SUPPLEMENTARY FIGURE 4
    if (produce_plots) {
      pdf(file.path(figuresPath, paste0("Venn_driver_", col, ".pdf")), 5, 5)
      print(
        ggvenn(
          data = list("Driver genes" = dg1, "Hits" = hits),
          text_size = 5, fill_color = c("#F0E442", "#0072B2"), show_percentage = FALSE
        )
      )
      dev.off()
    }
  }

  invisible(TRUE)
}

# Build matrices summarizing driver-gene hits across tissues.
# Each driver gene / tissue pair is labeled as:
# - "Known Found": known driver in that cancer type and detected here
# - "Novel": detected here but not previously known in that cancer type
# - "Known Not Found": known in that cancer type but not detected here
compute_driver_summary_matrices <- function(driver_genes, tissues, results, cancer_match_long_CMP, cancer_match_long_into, 
                                            intogen_drivers, dep_threshold, RR_th) {

  summary_drivers <- matrix(nrow = length(driver_genes), ncol = length(tissues))
  rownames(summary_drivers) <- driver_genes
  colnames(summary_drivers) <- tissues

  for (dg in driver_genes) {
    for (ctiss in tissues) {
      
      # Map CMP tissue name to the corresponding IntOGen cancer type(s)
      ct_into <- cancer_match_long_into[which(cancer_match_long_CMP == ctiss)]
      
      # Check whether the driver gene is detected as a DAM hit in this tissue
      if (dg %in% results[[ctiss]]$GENE[
        results[[ctiss]]$rank_ratio < RR_th &
          results[[ctiss]]$medFitEff < dep_threshold &
          results[[ctiss]]$hypTest_p < 0.2 &
          results[[ctiss]]$empPval < 0.2
      ]) {
        
        # Detected as a hit, initially classified as novel
        summary_drivers[dg, ctiss] <- "Novel"
        
        # Reclassify as known if the gene is an IntOGen driver in the matched cancer type
        if (length(intersect(ct_into, intogen_drivers$CANCER_TYPE[intogen_drivers$SYMBOL == dg])) > 0) {
          summary_drivers[dg, ctiss] <- "Known Found"
        }

      } else if (length(intersect(ct_into, intogen_drivers$CANCER_TYPE[intogen_drivers$SYMBOL == dg])) > 0) {
        # Known driver for this cancer type, but not detected in current DAM analysis
        summary_drivers[dg, ctiss] <- "Known Not Found"
      }
    }
  }
  
  # Convert categorical matrix to binary form for heatmap:
  # 1 = Known Found, 0 = Novel
  summary_drivers_bin <- summary_drivers
  summary_drivers_bin[summary_drivers_bin == "Known Found"] <- 1
  summary_drivers_bin[summary_drivers_bin == "Novel"] <- 0
  class(summary_drivers_bin) <- "numeric"
  
  # Keep drivers observed in more than one tissue
  selhits2 <- names(which(rowSums(summary_drivers == ("Known Found") | summary_drivers == ("Novel"), na.rm = TRUE) > 1))

  return(list(
    summary_drivers = summary_drivers,
    summary_drivers_bin = summary_drivers_bin,
    selhits2 = selhits2
  ))
}

# Plot a heatmap of recurrent driver-gene hits across tissues.
# Rows are drivers and columns are tissues.
plot_driver_heatmap <- function(figuresPath, summary_drivers_bin, selhits2, produce_plots = TRUE) {

  toplot <- summary_drivers_bin[selhits2, , drop = FALSE]

  toplot <- toplot[order(rowSums(!is.na(toplot)), decreasing = TRUE),
    order(colSums(!is.na(toplot)), decreasing = TRUE), drop = FALSE]
  
  # SUPPLEMENTARY FIGURE 6 pt.2
  if (produce_plots) {
    outfile <- file.path(figuresPath, "Drivers_pheat_novelvsknown.pdf")
    ph <- pheatmap::pheatmap(toplot, cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 15, cellheight = 15, na_col = "grey90", silent = TRUE)
    pdf(outfile, width = 15, height = 15)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()
  }

  invisible(TRUE)
}

# Plot, for recurrent driver genes, the number of cancer types in which they are:
# - found in a known cancer type
# - found in a novel cancer type
# - known in that cancer type but not found
plot_driver_barplot <- function(figuresPath, summary_drivers, tissues, produce_plots = TRUE) {
  
  selhits <- names(which(rowSums(summary_drivers == c("Known Found") | summary_drivers == c("Novel"), na.rm = TRUE) > 1))

  df <- data.frame(
    counts = c(rowSums(summary_drivers == "Known Found", na.rm = TRUE),
               rowSums(summary_drivers == "Known Not Found", na.rm = TRUE),
               rowSums(summary_drivers == "Novel", na.rm = TRUE)),
    type = rep(c("Found in a Known Cancer Type", "Not Found in a Known Cancer Type", "Found in a Novel Cancer Type"),
               each = nrow(summary_drivers)),
    gene = c(rownames(summary_drivers), rownames(summary_drivers), rownames(summary_drivers))
  )

  df$gene <- factor(df$gene, levels = c(df$gene[df$type == "Found in a Known Cancer Type"][
    order(df$counts[df$type == "Found in a Known Cancer Type"] + df$counts[df$type == "Found in a Novel Cancer Type"])
  ]))

  df$counts[df$type == "Not Found in a Known Cancer Type"] <- -df$counts[df$type == "Not Found in a Known Cancer Type"]
  df$type <- factor(df$type, levels = c("Found in a Known Cancer Type", "Found in a Novel Cancer Type", "Not Found in a Known Cancer Type"))

  # FIGURE 3B
  if (produce_plots) {
    pdf(file.path(figuresPath, "Drivers_barplot_thr2.pdf"), 10, 5)
    print(
      ggplot(subset(df, (gene %in% selhits) & (type != "Not Found in a Known Cancer Type")),
            aes(x = gene, y = counts, fill = type)) +
        geom_bar(stat = "identity") +
        theme_classic() +
        ylab("Number of cancer types") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 90, size = 10)) +
        scale_fill_manual(values = c("#0072B2", "#56B4E9"), labels = c("Known", "Novel")) +
        guides(fill = guide_legend(title = "Cancer Type"))
    )
    dev.off()
  }

  invisible(TRUE)
}

# Build a binary matrix showing in which tissues each non-driver hit is detected.
# Rows are genes, columns are tissues, and entries are 1 if the gene is a hit.
compute_nondriver_summary_matrix <- function(hits_nodriver, tissues, results, dep_threshold, RR_th) {

  summary_nodrivers <- matrix(0, nrow = length(hits_nodriver), ncol = length(tissues))
  rownames(summary_nodrivers) <- hits_nodriver
  colnames(summary_nodrivers) <- tissues

  for (dg in hits_nodriver) {
    for (ctiss in tissues) {
      # Mark presence if the non-driver gene is detected as a hit in this tissue
      if (dg %in% results[[ctiss]]$GENE[
        results[[ctiss]]$rank_ratio < RR_th &
          results[[ctiss]]$medFitEff < dep_threshold &
          results[[ctiss]]$hypTest_p < 0.2 &
          results[[ctiss]]$empPval < 0.2
      ]) {
        summary_nodrivers[dg, ctiss] <- 1
      }
    }
  }

  return(summary_nodrivers)
}

# Plot the number of tissues in which each non-driver hit appears,
# restricted to genes observed in more than two tissues.
plot_nondriver_ntissues <- function(figuresPath, summary_nodrivers, produce_plots = TRUE) {
  
  # Sort genes by the number of tissues in which they are detected
  df <- data.frame(DAMbg = names(sort(rowSums(summary_nodrivers), decreasing = TRUE)),
                   ntissues = sort(rowSums(summary_nodrivers), decreasing = TRUE))
  df$DAMbg <- factor(df$DAMbg, levels = c(names(sort(rowSums(summary_nodrivers), decreasing = TRUE))))
  
  # SUPPLEMENTARY FIGURE 6 pt.1
  if (produce_plots) {
    pdf(file.path(figuresPath, "DAMbgs_unreported_ntiss.pdf"), 12, 4)
    print(
      ggplot(subset(df, ntissues > 2), aes(x = DAMbg, y = ntissues)) +
        geom_bar(stat = "identity") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, size = 10)) +
        scale_fill_manual(values = c("#0072B2", "#E69F00")) +
        guides(fill = guide_legend(title = "Hit"))
    )
    dev.off()
  }

  invisible(TRUE)
}

# Plot stacked normalized hit counts per tissue, split into:
# - known driver hits
# - novel driver hits
# - non-driver hits
# Only the top tissues by normalized hit burden are shown.
plot_hits_tissue_stack <- function(figuresPath, summary_drivers, summary_nodrivers,
                                  tissues, mut_burden_new, num_hits_norm, produce_plots = TRUE) {

  df <- data.frame(
    counts = c(colSums(summary_drivers == "Known Found", na.rm = TRUE),
               colSums(summary_drivers == "Novel", na.rm = TRUE),
               colSums(summary_nodrivers == 1)),
    type = rep(c("Known Driver", "Known Driver", "Not Known as Driver"), each = length(tissues)),
    tissue = tissues
  )
  
  # Normalize hit counts by mutational burden
  df$counts <- df$counts / rep(mut_burden_new, 3)
  df$tissue <- factor(df$tissue, levels = c(df$tissue[order(num_hits_norm, decreasing = TRUE)]))
  df_sel <- df[df$tissue %in% df$tissue[order(num_hits_norm, decreasing = TRUE)][1:10], ]

  if (produce_plots) {
    # not assigned figure number 
    pdf(file.path(figuresPath, "Hits_tissue_stack_simpl.pdf"), 5, 5)
    print(
      ggplot(df_sel, aes(x = tissue, y = counts, fill = type)) +
        geom_bar(stat = "identity") +
        theme_classic() +
        ylab("% of hits") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 90, size = 10)) +
        scale_fill_manual(values = c("#0072B2", "#E69F00")) +
        guides(fill = guide_legend(title = "Hit"))
    )
    dev.off()
  }

  invisible(TRUE)
}

# Summarize how DAMs are distributed across cell lines and cancer types.
# Produces:
# - histogram of number of DAMs across cell lines
# - stacked barplot of DAM counts across cell lines by cancer type
# - barplot of percentage of cell lines with DAMs in each cancer type
# SUPPLEMENTARY FIGURE 8
plot_DAM_summary <- function(figuresPath, allDAMs, incl_cl_annot, produce_plots = TRUE) {

  ncell <- nrow(incl_cl_annot)
  ctiss <- unique(allDAMs$ctype)
  # Extract all cell line IDs listed in allDAMs$ps_cl
  ucls <- unique(unlist(str_split(allDAMs$ps_cl, ', ')))
  
  # number of DAMs observed per cell line
  nDAMsPerCCL <- unlist(lapply(ucls, function(U){
    length(grep(U, allDAMs$ps_cl))
  }))
  names(nDAMsPerCCL) <- ucls

  cat('DAMs are observed in ', length(ucls), ' cell lines (', 100 * length(ucls) / ncell, '%)\n', sep = '')

  # Plot histogram of number of DAMs per cell line
  # SUPPLEMENTARY FIGURE 8A
  if (produce_plots) {
    pdf(file.path(figuresPath,'n.DAMs across n. CCLs harbouring them.pdf'),6,6)
    hist(nDAMsPerCCL,100,
        main = 'n. DAMs across n. of CCLs in which they are observed',
        xlab='n. DAMs',
        ylab='n. CCLs')
    dev.off()
  }

  cat(length(which(nDAMsPerCCL < 4)) / length(nDAMsPerCCL), ' of CCLs with at least one DAM and no more than 3 DAMs\n')

  # Unique cancer types present in allDAMs
  utiss <- unique(allDAMs$ctype)
  
  # For each cancer type, build a vector of DAM counts across all its cell lines
  nDAMsAcrossCLS <- lapply(utiss, function(ut){

    # Cell lines belonging to the current cancer type
    currentCLS <- incl_cl_annot$model_id[
      which(incl_cl_annot$cancer_type == ut)
    ]

    nocc <- rep(0, length(currentCLS))
    names(nocc) <- currentCLS
    
    # Find which current cell lines have at least one DAM
    ii <- intersect(currentCLS, names(nDAMsPerCCL))
    nocc[ii] <- nDAMsPerCCL[ii]

    sort(nocc, decreasing = TRUE)
  })

  names(nDAMsAcrossCLS) <- utiss
  
  # Order cancer types by total DAM burden across their cell lines
  oo <- order(
    unlist(lapply(nDAMsAcrossCLS, 'sum')),
    decreasing = TRUE
  )

  nDAMsAcrossCLS <- nDAMsAcrossCLS[oo]

  # Pad with zeros to same length
  max_len <- max(lengths(nDAMsAcrossCLS))
  mat <- sapply(nDAMsAcrossCLS, function(x){
    c(x, rep(0, max_len - length(x)))
  })
  mat <- as.matrix(mat)
  
  # Count how many cell lines per cancer type have at least one DAM
  nCellLinesWithDAMs <- unlist(
    lapply(nDAMsAcrossCLS, function(x){
      length(which(x > 0))
    })
  )
  
  # Count how many total cell lines were screened in each cancer type
  nclsAcrossCtype <- unlist(
    lapply(names(nDAMsAcrossCLS), function(dd){
      length(which(incl_cl_annot$cancer_type == dd))
    })
  )
  
  # Stacked barplot: each bar is a cancer type, and the stack represents DAM counts across its cell lines
  # SUPPLEMENTARY FIGURE 8B
  if (produce_plots) {
    pdf(file.path(figuresPath,'n.DAMs_across_cell_lines_in_each_ctype.pdf'),15,8)
    par(mar=c(17,4,2,2))
    tt <- barplot(mat,
                  beside = FALSE,
                  col = c("skyblue", "salmon", "lightgreen", "orange"),
                  border = FALSE,
                  las = 2,
                  ylab='n. DAMs across cell lines',
                  ylim=c(0,280))

    text(tt, unlist(lapply(nDAMsAcrossCLS,sum)) + 5, nCellLinesWithDAMs)
    dev.off()
  }
  
  # Percentage of cell lines with DAMs in each cancer type
  # SUPPLEMENTARY FIGURE 8C
  if (produce_plots) {
    pdf(file.path(figuresPath,'perc_cell_lines_with_DAMs_in_each_ctype.pdf'),15,8)
    par(mar=c(17,4,2,2))
    tt <- barplot(100 * nCellLinesWithDAMs / nclsAcrossCtype,
                  las=2,
                  ylim=c(0,100),
                  border=FALSE,
                  ylab='% of cell lines with DAMs')

    text(tt, 100 * nCellLinesWithDAMs / nclsAcrossCtype + 2, nclsAcrossCtype)

    dev.off()
  }

}