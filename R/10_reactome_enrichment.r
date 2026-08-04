# -------- REACTOME ENRICHMENT ANALYSIS -------- #

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(VennDiagram)
  library(ReactomePA)
  library(reactome.db)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
library(biomaRt)

# Performs a hypergeometric test
my.hypTest <- function(x,k,n,N) {  
  
  PVALS <- phyper(x-1,n,N-n,k,lower.tail=FALSE)
  
  return(PVALS)
}

# Converts HGNC gene symbols to Entrez Gene IDs
# using a provided biomaRt Ensembl object.
# Returns unique non-NA Entrez IDs as character vector.
convert_symbols_to_entrez <- function(symbols, ensembl) { 

  options(timeout = 300) 
  res <- getBM(attributes = c('hgnc_symbol', 'entrezgene_id'), 
            filters = 'hgnc_symbol', 
            values = symbols, 
            mart = ensembl) 

  out <- as.character(res$entrezgene_id) 

  return(out) 
}

convert_symbols_to_entrez2 <- function(symbols, ensembl = NULL) {

  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols)]
  symbols <- symbols[symbols != ""]

  res <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("ENTREZID")
  )

  out <- unique(as.character(res$ENTREZID))
  out <- out[!is.na(out)]
  out <- out[out != ""]

  return(out)
}

# Performs Reactome pathway enrichment analysis using enrichPathway.
run_reactome_enrichment <- function(gene_entrez, universe_entrez, prefix, resultPath, figuresPath, pcut = 0.05, showCategory = 10, produce_plots = TRUE) {

  enr <- enrichPathway(gene = gene_entrez, pvalueCutoff = pcut, readable = TRUE, universe = universe_entrez)
  #enr_df <- as.data.frame(enr) # ??? not used 

  # write table
  out_file <- file.path(resultPath, paste0('_pws_enriched_in_', prefix, '_DAMbearing_genes.txt'))
  write.table(enr, quote = FALSE, sep = '\t', dec = ',', row.names = FALSE, file = out_file) 

  # dotplot
  # SUPPLEMENTARY FIGURE 9 
  if (produce_plots && nrow(as.data.frame(enr)) > 0) {
    pdf(file.path(figuresPath, paste0('REACTOME_enrich_', prefix, 'DAMbearing.pdf')), 
      width = 7, height = 6)
    print(dotplot(enr, showCategory = showCategory))
    dev.off()
  }

  return(enr)
}

# Compares Reactome enrichment results from all DAM-bearing genes and known DAM-bearing genes.
# Generates a Venn diagram of overlapping enriched pathways and writes two output tables:
# (1) conserved pathways shared between both analyses with genes split into known vs unreported,
# (2) pathways enriched only in the full (all genes) set with genes classified as IntOGen drivers vs others.
compare_enrichments_and_write <- function(all_enr, known_enr, resultPath, figuresPath, intOGen_drivers, produce_plots = TRUE) {

  # venn plot 
  # FIGURE 3C
  if (produce_plots) {
    pdf(file.path(figuresPath, "REACTOME_allDAMEnrich_vs_knownDAMEnrich.pdf"), 7, 7)
    venn.plot <- draw.pairwise.venn(
      area1 = length(known_enr$Description),
      area2 = length(all_enr$Description),
      cross.area = length(intersect(known_enr$Description, all_enr$Description)),
      category = c("kDAMbEnr", "aDAMbEnr"),
      fill = c("purple", "orange"),
      lty = "blank",
      cex = 2,
      cat.cex = 2,
      cat.pos = c(-20, 20)
    )
    dev.off()
  }

  conserved <- intersect(rownames(all_enr), rownames(known_enr))
  newOnly <- setdiff(rownames(all_enr), rownames(known_enr))

  # table for conserved enrichments: split genes into known/unreported
  if (length(conserved) > 0) {
    id1 <- match(conserved, all_enr$ID)
    id2 <- match(conserved, known_enr$ID)

    allG <- strsplit(all_enr[id1, 'geneID'], '/')
    knownG <- strsplit(known_enr[id2, 'geneID'], '/')
    unreportedG <- lapply(1:length(allG),function(x){setdiff(allG[[x]],knownG[[x]])})

    res <- cbind(all_enr[id1,1:3],
          unlist(lapply(1:length(knownG),function(x){paste(knownG[[x]],collapse=', ')})),
          unlist(lapply(1:length(unreportedG),function(x){paste(unreportedG[[x]],collapse=', ')})))

    colnames(res) <- c('pathway id','id','description','known DAM-bearing genes','unreported DAM-bearing genes')

    write.table(res,sep='\t',quote=FALSE,
                file=file.path(resultPath,'_pws_enriched_conserved_in_all_vs_known_DAMbearing_genes.txt'))
  }

  # table for newOnly pathways with gene breakdown
  if (length(newOnly) > 0) {
    id1 <- match(newOnly, all_enr$ID)
    id2 <- match(newOnly, known_enr$ID)

    allG <- strsplit(all_enr[id1,'geneID'],'/')
    knownG <- lapply(1:length(allG),function(x){intersect(allG[[x]], intOGen_drivers)})
    unreportedG <- lapply(1:length(allG),function(x){setdiff(allG[[x]],knownG[[x]])})

    res2 <- cbind(all_enr[id1,1:3],
                unlist(lapply(1:length(knownG),function(x){paste(knownG[[x]],collapse=', ')})),
                unlist(lapply(1:length(unreportedG),function(x){paste(unreportedG[[x]],collapse=', ')})))

    colnames(res2) <- c('pathway id','id','description','known DAM-bearing genes','unreported DAM-bearing genes')

    write.table(res2,sep='\t',quote=FALSE,file=file.path(resultPath,'_pws_enriched_setDiff_of_enrichments_all_vs_known_DAMbearing_genes.txt'))
  }

  # return the two pathway sets invisibly
  return(invisible(list(known_DAM_enrichments = known_enr,
                        all_DAM_enrichments = all_enr,
                        conserved = conserved,
                        newOnly = newOnly)))
}

# Prepares background genes and pathway mapping required for downstream enrichment simulations and co-occurrence analyses.
# Retrieves all Entrez genes, defines new background (excluding drivers),
# and returns gene universe and pathway-to-gene mappings.
prepare_background_data <- function(background, intOGen_drivers, ensembl) {
  
  all_genes <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), mart = ensembl)
  
  all_genes <- all_genes$entrezgene_id
  all_genes <- all_genes[!is.na(all_genes)]
  
  N <- length(all_genes)
  
  new_background <- setdiff(background, intOGen_drivers)
  
  new_background_entrez <- getBM(
    attributes = c("hgnc_symbol", "entrezgene_id"),
    filters = "hgnc_symbol",
    values = new_background,
    mart = ensembl
  )
  
  new_background_entrez <- new_background_entrez$entrezgene_id
  
  pathway_to_genes <- as.list(reactomePATHID2EXTID)
  
  return(list(all_genes = all_genes,
              N = N,
              new_background_entrez = new_background_entrez,
              pathway_to_genes = pathway_to_genes))
}

prepare_background_data2 <- function(background, intOGen_drivers, ensembl = NULL) {

  all_genes <- keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID")

  all_genes <- unique(as.character(all_genes))
  all_genes <- all_genes[!is.na(all_genes)]
  all_genes <- all_genes[all_genes != ""]

  N <- length(all_genes)

  new_background <- setdiff(background, intOGen_drivers)

  new_background_entrez <- convert_symbols_to_entrez2(new_background)

  pathway_to_genes <- as.list(reactome.db::reactomePATHID2EXTID)

  pathway_to_genes <- lapply(pathway_to_genes, function(x) {unique(as.character(x))})

  return(list(
    all_genes = all_genes,
    N = N,
    new_background_entrez = new_background_entrez,
    pathway_to_genes = pathway_to_genes
  ))
}

# Performs permutation-based Reactome enrichment simulation.
# Randomly samples genes preserving pathway membership proportions,
# computes overlap with known enrichments, and saves null distributions.
run_random_enrichment_simulation <- function(new_DAM_bearing_entrez, known_DAM_bearing_entrez, new_background_entrez, all_genes, nsim = 1000,
  pathway_to_genes, known_DAM_enrichments, background_entrez, OnlyInAllPaths, resultPath) {
  
  cleng_ <- NULL
  eleng_ <- NULL
  presencePath <- NULL
  
  # how many new DAM genes are present in Reactome
  n_new_DAM_entrez_in_Reactome <- length(intersect(new_DAM_bearing_entrez, unlist(pathway_to_genes)))
  # how many new DAM genes are not present in Reactome
  n_new_DAM_entrez_out_Reactome <- length(setdiff(new_DAM_bearing_entrez, unlist(pathway_to_genes)))
  
  # split all genes into those in Reactome vs not
  allGenes_in_reactome <- intersect(all_genes, unlist(pathway_to_genes))
  allGenes_out_reactome <- setdiff(all_genes, unlist(pathway_to_genes))
  
  set.seed(1234)
  
  for (i in 1:nsim) {
    
    # randomly sample genes that are in Reactome pathways
    randomGenesInPathways <- as.character(sample(intersect(new_background_entrez, allGenes_in_reactome), n_new_DAM_entrez_in_Reactome))
    # randomly sample genes that are not in Reactome pathways
    randomGenesOutPathways <- as.character(sample(intersect(new_background_entrez, allGenes_out_reactome), n_new_DAM_entrez_out_Reactome))
    
    randomGenes <- c(randomGenesInPathways, randomGenesOutPathways)
    
    random_DAM_enrichments <- enrichPathway(
      gene = c(known_DAM_bearing_entrez, randomGenes),
      pvalueCutoff = 0.05,
      readable = TRUE,
      universe = background_entrez)
    
    # compare pathways vs the known enrichment result
    commonRandpaths <- intersect(known_DAM_enrichments$Description, random_DAM_enrichments$Description)
    onlyNewRandPaths <- setdiff(random_DAM_enrichments$Description, known_DAM_enrichments$Description)
    
    cleng_[i] <- length(commonRandpaths)
    eleng_[i] <- length(onlyNewRandPaths)
    
    presencePath <- cbind(presencePath, unlist(lapply(OnlyInAllPaths, function(x){is.element(x, onlyNewRandPaths)})) + 0)
  }
  
  rownames(presencePath) <- OnlyInAllPaths

  return(list(cleng_ = cleng_,
              eleng_ = eleng_,
              presencePath = presencePath))
  
}

# Computes empirical p-values for conserved and new-only enrichments
# using permutation distributions.
# Generates histograms and pathway-level empirical frequency plots.
run_pathway_enrichment_empirical_tests <- function(cleng_, eleng_, Conserved_paths, newOnly, new_DAM_bearing_entrez,
    presencePath, figuresPath, produce_plots = TRUE) {
  
  # conserved
  # SUPPLEMENTARY FIGURE 11A 
  if (produce_plots) {
    pdf(file.path(figuresPath, "distr_of_conserved_pathEnrichments.pdf"), 5, 4)
    hist(cleng_,
        main = "Conserved pathway enrichments\n when adding randomly selected genes\nto the known DAM-bearing",
        xlab = "n. conserved enrichments")
    abline(v = length(Conserved_paths), col = "red")
    dev.off()
  }
  
  pval_conserved <- length(which(cleng_>=length(Conserved_paths)))/length(cleng_)
  
  cat("probability of having", length(Conserved_paths),
  "conserved enriched pathways when adding to the DAM bearing genes known to be cancer driver genes an additional",
  length(new_DAM_bearing_entrez),
  "genes not known to be driver and selected by random chance (preserving the ratio of genes in reactome pathways or out) =",
  pval_conserved, "\n")
  
  # new only
  # SUPPLEMENTARY FIGURE 11B
  if (produce_plots) {
    pdf(file.path(figuresPath, "distr_of_newOnly_pathEnrichments.pdf"), 5, 4)
    hist(eleng_,
        main = "New pathway enrichments\n when adding randomly selected genes\nto the known DAM-bearing",
        xlab = "n. of new enrichments")
    abline(v = length(newOnly), col = "red")
    dev.off()
  }
  
  pval_newOnly <- length(which(eleng_>=length(newOnly)))/length(eleng_)
  
  cat(
  "probability of having", length(newOnly),
  "newly enriched pathways when adding to the DAM bearing genes known to be cancer driver genes an additional",
  length(new_DAM_bearing_entrez),
  "genes not known to be driver and selected by random chance (preserving the ratio of genes in reactome pathways or out) =",
  pval_newOnly, "\n")
  
  # empirical frequency plot
  # SUPPLEMENTARY FIGURE 11C
  if (produce_plots) {
    pdf(file.path(figuresPath, "newOnly_pathEnrichments_empPval.pdf"), 11, 5)
    par(mar = c(6, 28, 0, 3))
    barplot(rowSums(presencePath) / ncol(presencePath),
            horiz = TRUE, las = 2, xlab = "p")
    dev.off()
  }
  
  return(list(
    pval_conserved = pval_conserved,
    pval_newOnly = pval_newOnly
  ))
}

# Tests whether new DAM-bearing genes co-occur in Reactome pathways
# with known cancer driver genes.
# Computes hypergeometric p-value, permutation null distribution,
# and generates summary plots.
run_driver_cooccurrence_analysis <- function(pathway_to_genes, IntOGen_Drivers_entrez, new_DAM_bearing_entrez, all_genes,
    new_gene_symbols, figuresPath, nperm = 1000, produce_plots = TRUE) {

  if (length(new_DAM_bearing_entrez) == 0) {
    cat("No unreported DAM-bearing genes available for Reactome co-occurrence analysis\n")
    return(invisible(NULL))
  }
  
  # identify pathways that contain at least one known driver
  has_driver <- unlist(lapply(pathway_to_genes, function(x){length(intersect(x,IntOGen_Drivers_entrez))>0}))
  
  # collect all genes that are in pathways containing a driver
  genes_with_driver <- unique(unlist(pathway_to_genes[which(has_driver)]))
  genes_with_driver <- intersect(genes_with_driver,all_genes)

  # hypergeometric test setup
  x <- length(intersect(new_DAM_bearing_entrez, genes_with_driver))
  k <- length(new_DAM_bearing_entrez)
  n <- length(genes_with_driver)
  N <- length(all_genes)
  
  res <- NULL
  set.seed(1234567)
  # null distribution by permutation
  for (i in 1:nperm){
    nn <- sample(all_genes,length(new_DAM_bearing_entrez))  
    res[i] <- length(intersect(nn,genes_with_driver))
  }
  
  observed_p <- my.hypTest(x, k, n, N)
  expectation <- mean(res)

  cat(x, " of the new DAM-bearing genes (",
  round(100 * x / k), "%) co-occurred in at least one pathway with known cancer driver gene\n",
  sep = "")
  cat("p =", observed_p, "\n")
  cat("expected value =", expectation, "\n")
  
  # plot
  # SUPPLEMENTARY FIGURE 12
  if (produce_plots) {
    pdf(file.path(figuresPath, "CoOcc_with_known_drivers_in_REACTOME_pathways.pdf"), 5,3)
    hist(res,main=paste('co-occurrences in at least one pathway\nwith a cancer driver gene'),
        xlab=paste('n. out of ',length(new_gene_symbols),'randomly selected genes (1,000 simulations)'),xlim = range(c(res, x)))
    abline(v=x,col='red')
    dev.off()
    
    simulated_log_p <- unlist(lapply(res, function(x) {-log10(my.hypTest(x, k, n, N))}))
    if (observed_p == 0) {
      observed_log_p <- max(simulated_log_p, na.rm = TRUE) + 1
      observed_label <- expression(-log[10](p[observed]) == Inf)
    } else {
      observed_log_p <- -log10(observed_p)
      observed_label <- paste0("observed = ", round(observed_log_p, 2))}

    pdf(file.path(figuresPath, "CoOcc_with_known_drivers_in_REACTOME_pathways_pvals.pdf"), 5,3)
    hist(simulated_log_p,
        main='empirical p-values across 1,000 simulations',xlab='-log10(p)',xlim = range(c(simulated_log_p, observed_log_p)))
    abline(v = observed_log_p, col = "red")
    text(x = observed_log_p, y = par("usr")[4] * 0.9, labels = observed_label, col = "red", pos = 2, cex = 0.7)
    dev.off()
  }
  
  # FIGURE 3H
  if (produce_plots) {
    pdf(file.path(figuresPath, "n.newDAMbgs_in_a_pathways_with_a_known_driver.pdf"), 6,6)
    pie(c(x,length(new_DAM_bearing_entrez)-x),col=c('#273c97','white'),main=paste(length(new_gene_symbols),'unreported DAMbgs'),labels = c(c(x,length(new_DAM_bearing_entrez)-x)))
    dev.off()
  }
  
  return(list(
    observed = x,
    expectation = expectation,
    pval = observed_p,
    null_distribution = res
  ))
}

# Quantifies pathway coverage increase after adding new DAM-bearing genes.
# Generates plots for conserved pathway coverage expansion
# and coverage of newly enriched pathways.
run_pathway_coverage_analysis <- function(known_DAM_enrichments, all_DAM_enrichments, Conserved_paths, newOnly, known_DAM_bearing_entrez,
    new_DAM_bearing_entrez, reactome.db, figuresPath, produce_plots = TRUE) {
  
  # two-column matrix of gene counts for conserved pathways
  Path_increasedCoverage <- cbind(
    known_DAM_enrichments[
      match(Conserved_paths, known_DAM_enrichments$Description), "Count"],
    all_DAM_enrichments[
      match(Conserved_paths, all_DAM_enrichments$Description), "Count"]
  )
  
  rownames(Path_increasedCoverage) <- Conserved_paths
  # keeps only pathways where the count increased when going from known to all
  Path_increasedCoverage <- Path_increasedCoverage[Path_increasedCoverage[,2] - Path_increasedCoverage[,1] > 0, , drop = FALSE]
  
  # plot 1: all pathways with increased coverage
  # SUPPLEMENTARY FIGURE 10
  oo <- order(Path_increasedCoverage[,2] - Path_increasedCoverage[,1])
  if (produce_plots && nrow(Path_increasedCoverage) > 0) {
    pdf(file.path(figuresPath, "path_increasedcoverage.pdf"), 10,15)
    par(mar=c(4,26,0,0.5))
    barplot(t(cbind(Path_increasedCoverage[oo,2]-Path_increasedCoverage[oo,1],Path_increasedCoverage[oo,1])),
            beside = FALSE,horiz = TRUE,las=2,xlab='n. genes',xlim=c(0,80),cex.names=0.5,border=FALSE,
            col=c('orange','purple'))
    dev.off()
  }

  # plot 2: top 10 pathways by increased coverage
  # FIGURE 3D
  oo <- order(Path_increasedCoverage[,2]-Path_increasedCoverage[,1],decreasing = TRUE)
  if (produce_plots && nrow(Path_increasedCoverage) > 0) {
    top_n <- min(10, length(oo))
    top_idx <- rev(oo[seq_len(top_n)])
    pdf(file.path(figuresPath, "path_increasedcoverage_top10.pdf"),10,4)
    par(mar=c(4,10,0,0.5))
    barplot(t(cbind(Path_increasedCoverage[top_idx, 1],
                    Path_increasedCoverage[top_idx, 2] - Path_increasedCoverage[top_idx, 1])),
            beside = FALSE,horiz = TRUE,las=2,xlab='n. genes',xlim=c(0,75),cex.names=0.5,border=FALSE,
            col=c('purple','orange'))
    dev.off()
  }
  
  # compute gene composition for “newOnly” pathways
  if (length(newOnly) > 0) {
    res <- do.call(rbind,lapply(1:length(newOnly),function(x) {
      pathToEntrez <- suppressMessages(AnnotationDbi::select(reactome.db,
                          keys = newOnly[x],
                          keytype = "PATHID",
                          columns = c("ENTREZID")))

      knownG <- length(intersect(known_DAM_bearing_entrez,pathToEntrez$ENTREZID))
      newG <- length(intersect(new_DAM_bearing_entrez,pathToEntrez$ENTREZID))
      return(c(knownG,newG))
      }))

    rownames(res) <- all_DAM_enrichments[match(newOnly,all_DAM_enrichments$ID),'Description']

    oo <- order(res[,2])

    # plot 3: “NewOnly” pathway gene counts
    # FIGURE 3J
    if (produce_plots) {
      pdf(file.path(figuresPath, "path_increasedcoverage_NewOnly.pdf"),10,4)
      par(mar=c(4,10,0,0.5))
      barplot(t(res[oo,]),
              beside = FALSE,horiz = TRUE,las=2,xlab='n. genes',xlim=c(0,60),cex.names=0.5,border=FALSE,
              col=c('purple','orange'))
      dev.off()
    }
  }
  
  return(Path_increasedCoverage)
}
