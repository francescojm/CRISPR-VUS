# -------- CO-OCCURRENCE ANALYSIS -------- #

library(tidyverse)

# Compute co-occurrence per cell line
# For each cell line carrying at least one DAM in an unreported DAM background, this function summarizes:
# - the DAMs present in that cell line
# - which of those DAMs are druggable
# - whether the cell line also carries tissue-specific activating driver mutations
# - whether those activating drivers are essential
compute_DAM_oncogenic_cooccurrence <- function(allDAMs, bdep, CMP_annot, cl_variants, driver_genes,
                                               intOGen_drivers, mapping, tractable_targets) {
  
  # DAM-bearing lines
  DAM_bearing_lines <- unique(unlist(strsplit(allDAMs$ps_cl[-which(allDAMs$GENE %in% driver_genes)], ", ")))
  known_DAM_bearing_lines <- unique(unlist(strsplit(allDAMs$ps_cl[which(allDAMs$GENE %in% driver_genes)], ", ")))
  
  RES <- do.call(rbind, lapply(DAM_bearing_lines, function(l) {
    
    # Essential genes in this cell line
    ess_genes <- names(which(bdep[, l] == 1))
    
    # Cancer type and tissue-specific activating driver genes
    ct <- CMP_annot$cancer_type[CMP_annot$model_id == l]
    ct_into <- unlist(strsplit(mapping[ct, 1], " \\| "))
    act_genes <- intOGen_drivers$SYMBOL[which(intOGen_drivers$ROLE == "Act" & is.element(intOGen_drivers$CANCER_TYPE, ct_into))]
    
    # Mutated genes in this cell line
    all_Mutated_Genes <- cl_variants$gene_symbol_2023[which(cl_variants$model_id == l)]
    # Keep activating drivers
    all_Mutated_act_Genes <- intersect(all_Mutated_Genes, act_genes)
    # Keep essential genes
    Act_Ess_Drivers <- intersect(all_Mutated_act_Genes, ess_genes)
    
    # Observed DAMs in this line
    idcell <- grep(l, allDAMs$ps_cl)
    gg <- allDAMs$GENE[idcell]
    vars <- allDAMs$var[idcell]
    observed_DAMs <- paste(gg, vars, sep = '')
    
    # Split DAMs by driver/tractability
    DAMs_in_unreported_DAMbgs <- paste(observed_DAMs[which(!is.element(allDAMs$GENE[idcell], driver_genes))], collapse = ', ')
    druggable_DAMs_in_unreported_DAMbgs <- paste(intersect(setdiff(gg, driver_genes), tractable_targets), collapse = ', ')
    DAMs_in_known_DAMbgs <- paste(observed_DAMs[which(is.element(allDAMs$GENE[idcell], act_genes))], collapse = ', ')
    druggable_known_DAMs <- paste(intersect(intersect(gg, act_genes), tractable_targets), collapse = ', ')
    druggable_Act_Ess_Drivers <- paste(intersect(Act_Ess_Drivers, tractable_targets), collapse = ', ')
    Act_Ess_Drivers <- paste(Act_Ess_Drivers, collapse = ', ')
    
    c(ct, l, DAMs_in_unreported_DAMbgs, druggable_DAMs_in_unreported_DAMbgs,
      DAMs_in_known_DAMbgs, druggable_known_DAMs, Act_Ess_Drivers, druggable_Act_Ess_Drivers)
    
  }))
  
  colnames(RES) <- c('cancer_type','cell_line','DAMs_in_unreported_DAMbgs','druggable_unreported',
                     'DAMs_in_known_DAMbgs','druggable_known','GoF_mutated_essential_drivers','druggable_Act_Ess_Drivers')
  
  return(RES)
}

# Compute summary statistics and co-occurrence by tissues.
# This function summarizes how often cell lines with DAMs in unreported DAM
# backgrounds also show evidence of co-occurring oncogenic addiction.
summarize_cooccurrence <- function(RES, figuresPath, plot_prefix = "unreportedDAMs_oncogenicAddiction", produce_plots = TRUE) {
  
  # Percentages of missing DAMs in known drivers
  nn1 <- 100*length(which(RES[,5] == ''))/nrow(RES)
  cat(sprintf("%.2f %% of CCLs with an unreported DAM lack DAMs in tissue-specific known GoF drivers\n", nn1))
  
  nn2 <- 100*length(which(RES[,7] == ''))/nrow(RES)
  cat(sprintf("%.2f %% of CCLs with an unreported DAM lack tissue-specific essential GoF drivers\n", nn2))
  
  # Compute co-occurrence matrix
  occs <- unlist(lapply(RES[,5], function(x){x==''})) + 0
  tissues <- unique(RES[,1])
  
  # For each tissue, count:
  # - o: number of cell lines with co-occurring oncogenic addiction
  # - s: number of cell lines without it
  coc_results <- do.call(rbind, lapply(tissues, function(x){
    temp <- occs[which(RES[,1] == x)]
    n <- length(temp)
    s <- sum(temp)
    o <- n - s
    c(o, s)
  }))
  rownames(coc_results) <- tissues
  colnames(coc_results) <- c('co-occurring_oncAdd','other')
  
  cat(sprintf("Median across cancer types = %.2f\n", median(100 * coc_results[,2] / rowSums(coc_results))))
  
  coc_results <- coc_results[order(rowSums(coc_results), decreasing = TRUE), ]
  
  # Plot barplots
  # SUPPLEMENTARY FIGURE 13A
  if (produce_plots) {
    pdf(file.path(figuresPath, paste0(plot_prefix, "_coOcc.pdf")), 11,7)
    par(mar=c(16,4,2,2))
    barplot(t(coc_results), border=FALSE, las=2,
            main = paste(sum(c(coc_results)),'cell lines with DAM in unreported DAMbgs'),
            ylab = 'n. cell lines', col = c('darkgray','blue'))
    legend('topright', c('co-occurrent oncogenic addiction','others'),
          fill = c('darkgray','blue'), border = FALSE)
    dev.off()
  }
  
  # Repeat for relaxed criterion (column 7)
  occs <- unlist(lapply(RES[,7], function(x){x==''})) + 0
  coc_results <- do.call(rbind, lapply(tissues, function(x){
    temp <- occs[which(RES[,1] == x)]
    n <- length(temp)
    s <- sum(temp)
    o <- n - s
    c(o, s)
  }))
  rownames(coc_results) <- tissues
  colnames(coc_results) <- c('co-occurring_oncAdd','other')

  cat(sprintf("Median across cancer types = %.2f\n", median(100 * coc_results[,2] / rowSums(coc_results))))
  coc_results <- coc_results[order(rowSums(coc_results), decreasing = TRUE), ]
  
  # SUPPLEMENTARY FIGURE 13B
  if (produce_plots) {
    pdf(file.path(figuresPath, paste0(plot_prefix, "_coOcc_relaxed_criterion.pdf")), 11,7)
    par(mar=c(16,4,2,2))
    barplot(t(coc_results), border=FALSE, las=2,
            main = paste(sum(c(coc_results)),'cell lines with DAM in unreported DAMbgs'),
            ylab = 'n. cell lines', col = c('darkgray','blue'))
    legend('topright', c('co-occurrent oncogenic addiction','others'),
          fill = c('darkgray','blue'), border = FALSE)
    dev.off()
  }
  
  return(coc_results)
}