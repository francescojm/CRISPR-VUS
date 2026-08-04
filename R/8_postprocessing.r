# -------- POST-PROCESSING, SUMMARY, AND ENRICHMENT ANALYSES -------- #

library(stringr)
library(openxlsx)
library(ggplot2)

# Performs a hypergeometric test
my.hypTest <- function(x,k,n,N) {  
  
  PVALS <- phyper(x-1,n,N-n,k,lower.tail=FALSE)
  
  return(PVALS)
}

# split multi-variant hits into one row per variant
decoupleMultipleHits <- function(hitTable) {
  vars <- hitTable$var
  iimultiple <- grep(' | ', vars)

  # new data frame containing only rows that do not have multiple variants
  finalHit <- hitTable[setdiff(1:nrow(hitTable), iimultiple), ]

  finalHit <- rbind(finalHit,
                    do.call('rbind', lapply(iimultiple, function(x) {
                      # splits the variant string into individual variants and extract the resulting vector 
                      indivar <- setdiff(strsplit(vars[x], ' | ')[[1]], '|')
                      tmpHit <- NULL
                      # for each individual variant duplicate the original row and append to tmpHit
                      for (i in 1:length(indivar)) {
                        tmpHit <- rbind(tmpHit, hitTable[x, ])
                      }
                      # replaces the var column in the duplicated rows with the individual variant values
                      tmpHit$var <- indivar
                      return(tmpHit)
                    })))
  return(finalHit)
}

# summarize included cell lines across cancer types
# reports basic statistics on the number of cell lines per cancer type
# and generates a horizontal barplot of cell line counts
summarize_cell_lines <- function(incl_cl_annot, figuresPath, produce_plots = TRUE) {
  cat(nrow(incl_cl_annot), "cell lines included in the analysis\n")
  cat(length(unique(incl_cl_annot$cancer_type)), "considered cancer types\n")

  # number of cell lines per cancer type 
  ss <- summary(as.factor(incl_cl_annot$cancer_type))
  cat("median number of cell lines per cancer type =", median(ss), "\n")
  cat("min =", min(ss), "for", names(sort(ss))[1], "\n")
  cat("max =", max(ss), "for", names(sort(ss, decreasing=TRUE))[1], "\n")

  # SUPPLEMENTARY FIGURE 1A
  if (produce_plots) {
    pdf(file.path(figuresPath, 'Ncelllines.pdf'), 7,10)
    par(mar=c(4,16,0,2))
    barplot(sort(ss), horiz = TRUE, las=2, border=FALSE, col='blue', xlab='n. cell lines')
    abline(v=median(ss), lty=2)
    dev.off()
  }

  return(ss)
}

# collect and aggregate all tested variants across cancer types
# loads all *_testedVariants.RData files in the results directory, merges
# them into a single table, and collapses identical variants by reporting
# all positive cell lines for each variant
# the resulting table is saved as both .RData and .tsv
collect_tested_variants <- function(resultPath) {
  fc <- dir(resultPath)
  fc <- grep('testedVariants.RData', fc, value=TRUE)
  totalTestedVariants <- NULL
  for (i in 1:length(fc)) {
    # extracts the cancer type name from the filename
    cty <- strsplit(fc[i],'_testedVariants.RData')[[1]]
    load(paste(resultPath,'/',fc[i],sep=''))

    # creates a unique identifier for each variant combining gene_symbol and protein_mutation
    variant_UMIs <- paste(ts_cl_variants$gene_symbol, ts_cl_variants$protein_mutation)
    unique_variants <- unique(variant_UMIs)

    curRes <- do.call('rbind', lapply(unique_variants, function(x) {
       # rows corresponding to a specific variant
       idxs <- which(variant_UMIs == x)
       # collects all model_id where this variant appears and concatenates them into a comma-separated string
       positiveCls <- paste(ts_cl_variants$model_id[idxs], collapse=', ')
       # summarized row 
       res <- cbind(cty, x, ts_cl_variants[idxs[1], c(1,2,4,5,6,7,8,9,10,11,12)], positiveCls)
       return(res)
    }))
    totalTestedVariants <- rbind(totalTestedVariants, curRes)
  }

  return(totalTestedVariants)
}

# summarize tested variants and generate diagnostic plots
# computes:
# - number of tested variant x cancer type combinations
# - number of unique variants
# - number of tested genes
# - number of cell lines carrying each variant
# also produces histograms of the number of cell lines per variant
summarize_tested_variants <- function(totalTestedVariants, figuresPath, tissues, produce_plots = TRUE) {

  # each row is a variant–cancer type combination
  ntestedVarCtypeCombos <- nrow(totalTestedVariants)
  cat(ntestedVarCtypeCombos, "tested variants x cancer type combos\n")

  # unique variants across cancer types
  ntestedIndividualVariants <- length(unique(totalTestedVariants$x))
  cat(ntestedIndividualVariants, "tested individual variants\n")

  # unique tested genes 
  ntestedGenes <- length(unique(totalTestedVariants$gene_symbol))
  cat("involving", ntestedGenes, "genes\n")

  # number of cell lines per variant 
  var_nlines <- str_count(totalTestedVariants$positiveCls, ',') + 1
  totalTestedVariants$var_nlines <- var_nlines
  
  # SUPPLEMENTARY FIGURE 1B
  if (produce_plots) {
    pdf(file.path(figuresPath, 'var_nlines.pdf'), 7,5)
    print(ggplot(totalTestedVariants, aes(x = var_nlines)) +
      geom_histogram(bins = 30, color = '#000000', fill = 'steelblue') +
      theme_classic() +
      labs(x = '# of cell lines bearing the variant', y = '# of variants'))
    dev.off()

    pdf(file.path(figuresPath, 'var_nlines_log10.pdf'), 7,5)
    print(ggplot(totalTestedVariants, aes(x = var_nlines)) +
      geom_histogram(bins = 30, color = '#000000', fill = 'steelblue') +
      theme_classic() +
      scale_y_log10() +
      labs(x = '# of cell lines bearing the variant', y = '# of variants (log10)'))
    dev.off()
  }

  # how many tested variants belong to a specific cancer type
  ct_tested_variants <- unlist(lapply(tissues, function(x) {
    length(which(totalTestedVariants$cty == x))
  }))
  names(ct_tested_variants) <- tissues
  return(list(totalTestedVariants = totalTestedVariants,
              ct_tested_variants = ct_tested_variants,
              ntestedGenes = ntestedGenes))
}

# collect all DAMs and hits across tissues
# For each *_results.RData file:
# - selects significant hits based on predefined thresholds
# - decouples multi-variant hits
# - harmonizes positive cell lines with totalTestedVariants
# - builds a gene-level table of DAM-bearing genes with IntOGen info
# Saves:
# - allHits (.RData and .tsv)
# - allDAMs (.RData and .tsv)
# - allDAM_bearing_genes (.RData and .tsv)
collect_all_DAMs <- function(resultPath, totalTestedVariants, RR_th, dep_threshold, empPval = 0.20, hypTest_p = 0.20) {
  fc <- dir(resultPath)
  fc <- grep('_results.RData', fc, value=TRUE)
  allDAMs <- NULL
  allHits <- NULL

  # loop over result files
  for (i in 1:length(fc)){
    # extract cancer type
    cty <- strsplit(fc[i],'_results.RData')[[1]]
    load(paste(resultPath, "/", fc[i], sep = ""))

    # select significant hits
    hitsIdxs <- which(RESTOT$medFitEff < dep_threshold & RESTOT$rank_ratio < RR_th &
                      RESTOT$empPval < empPval & RESTOT$hypTest_p < hypTest_p)
    currHits <- RESTOT[hitsIdxs,]
    currDAMs <- decoupleMultipleHits(hitTable = currHits)

    # where ps_cl contains multiple cell lines
    idxs <- grep(', ', currDAMs$ps_cl)
    positiveCLs <- lapply(idxs, function(iii){
      # extract gene and variant
      currGene <- currDAMs$GENE[iii]
      currVar  <- currDAMs$var[iii]
      # extract listed positive cell lines 
      currPsCls <- unlist(str_split(currDAMs$ps_cl[iii], ', '))

      # find corresponding tested variant entry
      subi <- which(totalTestedVariants$gene_symbol == currGene &
                    totalTestedVariants$protein_mutation == currVar)
      # keeps only cell lines that are positive in current hit
      PSCl <- unlist(str_split(totalTestedVariants$positiveCls[subi], ', '))
      PSCl <- paste(intersect(PSCl, currPsCls), collapse=', ')
      return(PSCl)
    })

    currDAMs$ps_cl[idxs] <- unlist(positiveCLs)
    rownames(currHits) <- NULL
    rownames(currDAMs) <- NULL

    currHits <- currHits[order(currHits$GENE),]
    currDAMs <- currDAMs[,c(1:3,11)]
    currDAMs <- currDAMs[order(currDAMs$GENE),]
    allDAMs <- rbind(allDAMs, currDAMs)
    allHits <- rbind(allHits, currHits)
  }

  # get unique genes with significant hits
  allDAM_bearing_genes <- sort(unique(allHits$GENE))
  id <- match(allDAM_bearing_genes, allHits$GENE)
  # for each gene collect cancer types
  DAMbearing_in <- unlist(lapply(allDAM_bearing_genes, function(x) {
    cid <- which(allHits$GENE == x)
    analyses <- paste(sort(allHits$ctype[cid]), collapse=' | ')
  }))
  
  # combine into summary table
  allDAM_bearing_genes <- cbind(allDAM_bearing_genes,
    allHits[id, c('intOGen_driver_for','Act','LoF','ambigous')],
    DAMbearing_in)

  return(list(allHits = allHits,
              allDAMs = allDAMs,
              allDAM_bearing_genes = allDAM_bearing_genes))
}

# global summary of DAMs and hits
# prints descriptive statistics on:
# - total number of significant hits
# - total number of DAMs
# - number of unique variants
# - fraction of tested cases that are DAMs
# - number and fraction of genes affected
summarize_DAMs <- function(allHits, allDAMs, totalTestedVariants, ntestedGenes) {
  cat(nrow(allHits), "hits with rankratio, medFitness effect, emp pvalue, and HG pvalue below the thresholds (default: 1.71, -0.5, 0.2, and 0.2 respectively), across cancer types\n")
  cat(nrow(allDAMs), "individual DAMs across cancer types\n")
  cat("corresponding to", length(unique(paste(allDAMs$GENE, allDAMs$var))), "individual variants\n")
  cat(100 * length(paste(allDAMs$GENE, allDAMs$var)) / nrow(totalTestedVariants),
      "% of tested cases (considering variants tested multiple times across different cancer types)\n")
  cat("involving", length(unique(allDAMs$GENE)), "genes\n")
  cat(100 * length(unique(allDAMs$GENE)) / ntestedGenes, "% of tested cases\n")
}

# plot DAM counts by cancer type
# Generates two barplots:
# - number of DAMs per cancer type
# - number of DAM-bearing genes per cancer type
plot_DAMs_by_ctype <- function(allDAMs, clc, figuresPath, produce_plots = TRUE) {

  # SUPPLEMENTARY FIGURE 2A 
  if (produce_plots) {
    pdf(file.path(figuresPath, 'DAMs_ctype.pdf'), 12,10)
    par(mar=c(25,5,2,5))
    barplot(table(allDAMs$ctype)[order(table(allDAMs$ctype), decreasing=TRUE)], las=2,
            ylab='# of DAMs',
            col=clc[names(sort(table(allDAMs$ctype), decreasing=TRUE)),1],
            border=FALSE, cex.names = 1.5)
    invisible(dev.off())
    
    # SUPPLEMENTARY FIGURE 2B
    pdf(file.path(figuresPath, 'DAMsbearing_ctype.pdf'), 12,10)
    par(mar=c(25,5,2,5))
    allDAMbearing <- allDAMs[!duplicated(allDAMs[,c('ctype','GENE')]),]
    barplot(table(allDAMbearing$ctype)[order(table(allDAMbearing$ctype), decreasing=TRUE)], las=2,
            ylab='# of DAM-bearing genes',
            col=clc[names(sort(table(allDAMbearing$ctype), decreasing=TRUE)),1],
            border=FALSE, cex.names = 1.5)
    invisible(dev.off())
  }
}

# compute per-tissue DAM statistics
# For each cancer type, calculates:
# - number of unique DAM-bearing genes
# - total number of DAMs
# Also prints median, min, and max across tissues
compute_DAM_stats <- function(allDAMs, tissues) {
  nDAMbearing <- c()
  for(t in tissues){
    allDAMs_sel <- allDAMs[allDAMs$ctype == t,]
    nDAMbearing <- c(nDAMbearing, length(unique(allDAMs_sel$GENE)))
  }
  names(nDAMbearing) <- tissues

  cat(median(nDAMbearing), "median number of DAM bearing genes across cancer types\n")
  cat("min =", min(nDAMbearing), "for", names(sort(nDAMbearing))[1], "\n")
  cat("max =", max(nDAMbearing), "for", names(sort(nDAMbearing, decreasing=TRUE))[1], "\n")

  nDAMs <- c()
  for(t in tissues) {
    allDAMs_sel2 <- allDAMs[allDAMs$ctype == t,]
    nDAMs <- c(nDAMs, length(allDAMs_sel2$GENE))
  }
  names(nDAMs) <- tissues

  cat(median(nDAMs), "median number of DAMs across cancer types\n")
  cat("min =", min(nDAMs), "for", names(sort(nDAMs))[1], "\n")
  cat("max =", max(nDAMs), "for", names(sort(nDAMs, decreasing=TRUE))[1], "\n")

  return(list(nDAMbearing = nDAMbearing, nDAMs = nDAMs))
}

# global enrichment of IntOGen cancer drivers among DAM-bearing genes
# performs a hypergeometric test to assess whether DAM-bearing genes
# are enriched for known cancer drivers from IntOGen
intOGen_global_enrichment <- function(intogen_drivers, allDAM_bearing_genes, ntestedGenes) {

  N <- ntestedGenes
  n <- length(intogen_drivers)
  k <- length(unique(allDAM_bearing_genes$allDAM_bearing_genes))
  x <- length(intersect(intogen_drivers, allDAM_bearing_genes$allDAM_bearing_genes))

  cat("Of the", k, "DAM bearing genes,", x, "have been previously reported as cancer drivers (intOGen) (", round(100 * x / k, 2), "%)\n")
  cat("Hypergeometric test p-value:\n")
  pval <- my.hypTest(x, k, n, N)
  cat(pval, "\n")

  return(list(x = x, k = k, N = N, n = n, pval = pval))
}

# role-specific IntOGen enrichment (oncogenes vs TSG vs ambiguous)
# Computes the predominant role (Act/LoF/ambiguous) for each IntOGen gene
# Classifies genes as:
#  - oncogenes (Act)
#  - tumor suppressors (LoF)
#  - ambiguous
# Tests enrichment of each class among DAM-bearing genes
# Generates a pie chart summarizing the composition
intOGen_role_enrichment <- function(intogen_drivers, allDAM_bearing_genes, k, N, x_global, figuresPath, produce_plots = TRUE) {

  ug <- sort(unique(intogen_drivers$SYMBOL))

  # for each gene compute proportion of roles
  roles <- do.call('rbind', lapply(ug, function(x) {
    ids <- which(intogen_drivers$SYMBOL == x)
    intogen_Role_Perc <- c(0,0,0)
    names(intogen_Role_Perc) <- c('Act','LoF','ambiguous')
    rolSums <- summary(as.factor(intogen_drivers$ROLE[ids]))
    intogen_Role_Perc[names(rolSums)] <- round(100*rolSums/sum(rolSums), 2)
    return(intogen_Role_Perc)
  }))

  rownames(roles) <- ug
  putRole <- rep(NA, length(ug))
  # for each gene, pick the role with the highest percentage and store 
  for (i in 1:nrow(roles)) {
    putRole[i] <- colnames(roles)[order(roles[i,], decreasing = TRUE)[1]]  
  }
  names(putRole) <- ug

  # classify genes in three categories
  oncogenes <- names(which(roles[,'Act'] == 100))
  tsg <- names(which(roles[,'LoF'] == 100))
  ambig <- setdiff(rownames(roles), c(oncogenes, tsg))

  # Act (oncogenes)
  n <- length(oncogenes)
  x <- length(intersect(oncogenes, allDAM_bearing_genes$allDAM_bearing_genes))
  pct <- if (k == 0) NA else round(100 * x / k, 2)
  cat("Of the", k, "DAM bearing genes,", x, "have been previously reported as 100% oncoGenes (intOGen) (", pct, "%)\n")
  cat("Hypergeometric test p-value:\n")
  pAct <- my.hypTest(x, k, n, N)
  cat(pAct, "\n")
  nAct <- x

  # LoF (TSG)  
  n <- length(tsg)
  x <- length(intersect(tsg, allDAM_bearing_genes$allDAM_bearing_genes))
  pct <- if (k == 0) NA else round(100 * x / k, 2)
  cat("Of the", k, "DAM bearing genes,", x, "have been previously reported as 100% TSG (intOGen) (", pct, "%)\n")
  cat("Hypergeometric test p-value:\n")
  pTsg <- my.hypTest(x, k, n, N)
  cat(pTsg, "\n")
  nTsg <- x

  # Ambiguous
  n <- length(ambig)
  x <- length(intersect(ambig, allDAM_bearing_genes$allDAM_bearing_genes))
  pct <- if (k == 0) NA else round(100 * x / k, 2)
  cat("Of the", k, "DAM bearing genes,", x, "have been previously reported as TSGs or oncoGenes (intOGen) (", pct, "%)\n")
  cat("Hypergeometric test p-value:\n")
  pAmb <- my.hypTest(x, k, n, N)
  cat(pAmb, "\n")
  nAmb <- x

  # SUPPLEMENTARY FIGURE 5A
  if (produce_plots && sum(c(nTsg, nAmb, nAct)) > 0) {
    pdf(file.path(figuresPath, "pie_drivers.pdf"), 7,8)
    pie(c(nTsg,nAmb,nAct),
        col=c('#004add','#d5aff3','red'),
        border=FALSE,
        main=paste(x_global, paste('DAM-bearing genes known as cancer driver\nacross cancer types')),
        labels = c('Frequent LoF','Ambigous','Frequent GoF'),
        cex=1.2)
    dev.off()
  }

  return(list(nAct = nAct, pAct = pAct,
              nTsg = nTsg, pTsg = pTsg,
              nAmb = nAmb, pAmb = pAmb,
              roles = roles,
              putRole = putRole))
}

# tissue-specific composition of DAM-bearing genes by IntOGen role
# for each cancer type:
# - maps to corresponding IntOGen cancer types
# - classifies DAM-bearing genes as TSG, ambiguous, oncogene, or novel
# - performs tissue-specific hypergeometric enrichment tests
# Generates summary barplots and enrichment scatter plots
composition_analysis <- function(intogen_drivers, ct_mapping, figuresPath, tissues,
                                 allHits, totalTestedVariants, produce_plots = TRUE) {

  COMPOSITION <- NULL
  COMPOSITIONp <- NULL
  CT_DAM_Bearing <- NULL

  for (i in 1:length(tissues)) {
    # extract the DAM-bearing gene symbols for each tissue and store
    ctiss <- tissues[i]
    current_Dam_bearing <- allHits$GENE[which(allHits$ctype == ctiss)]
    CT_DAM_Bearing[[i]] <- current_Dam_bearing

    intoTypes <- setdiff(unlist(strsplit(ct_mapping[ctiss,1],' | ')),'|')   #intogen cancer type identifier
    currentOG <- unique(intogen_drivers[which(is.element(intogen_drivers$CANCER_TYPE,intoTypes) & intogen_drivers$ROLE=='Act'),'SYMBOL'])   #oncogenes from Intogen in the specific cancer type
    currentTSG <- unique(intogen_drivers[which(is.element(intogen_drivers$CANCER_TYPE,intoTypes) & intogen_drivers$ROLE=='LoF'),'SYMBOL'])  #tumor suppressor genes from Intogen in the specific cancer type
    both <- intersect(currentOG,currentTSG)                                 #ambiguous OG or TSG from Intogen in the specific cancer type
    currentOG <- setdiff(currentOG,both)                                    #oncogenes (non ambiguous) from Intogen in the specific cancer type
    currentTSG <- setdiff(currentTSG,both)                                  #TSG (non ambiguous) from Intogen in the specific cancer type

    # classify the DAM-bearing genes into TSG / Amb / OG / Novel and store counts
    classOG <- intersect(current_Dam_bearing, currentOG)
    classAMB <- intersect(current_Dam_bearing, both)
    classTSG <- intersect(current_Dam_bearing, currentTSG)
    novel <- setdiff(current_Dam_bearing, c(classOG, classAMB, classTSG))
    RES <- c(length(classTSG), length(classAMB), length(classOG), length(novel))
 
    bg <- length(unique(totalTestedVariants$gene_symbol[totalTestedVariants$cty == ctiss]))  #number of tested genes in the specific tissue
 
    # run hypergeometric tests for enrichment of each role
    OGp  <- my.hypTest(length(classOG),  length(current_Dam_bearing), length(currentOG), bg)
    AMBp <- my.hypTest(length(classAMB), length(current_Dam_bearing), length(both), bg)
    TSGp <- my.hypTest(length(classTSG), length(current_Dam_bearing), length(currentTSG), bg)
    RESp <- c(TSGp, AMBp, OGp, NA)

    COMPOSITION  <- cbind(COMPOSITION, RES)
    COMPOSITIONp <- cbind(COMPOSITIONp, RESp)
  }

  names(CT_DAM_Bearing) <- tissues

  # ordering by total counts 
  oo <- order(colSums(COMPOSITION))
  colnames(COMPOSITION) <- tissues
  colnames(COMPOSITIONp) <- tissues
  rownames(COMPOSITION) <- c('TSG','Amb','OG','Novel')
  rownames(COMPOSITIONp) <- c('TSG','Amb','OG','Novel')
  
  # SUPPLEMENTARY FIGURE 5B
  if (produce_plots) {
    pdf(file.path(figuresPath, 'ActLoFenrichment.pdf'), 12,9)
    par(mfrow=c(1,3))
    par(mar=c(4,12,1,2))
    barplot(colSums(COMPOSITION)[oo], horiz=TRUE, las=2,
            xlab='n. DAM-bearing genes', xlim=c(0,200), cex.names=0.6, border=NA,
            col='#75b4d9')
    abline(v= median(colSums(COMPOSITION)), lty=2, col='darkgray')

    barplot(100*COMPOSITION[1:3,oo,drop = FALSE]/t(matrix(rep(colSums(COMPOSITION[,oo,drop = FALSE]),3),
                                            ncol(COMPOSITION), 3)),
            horiz=TRUE, las=2, xlab='% DAM-bearing genes', xlim=c(0,25),
            cex.names=0.6, border=NA,
            col=c('#004add','#d5aff3','red'))

    COMPOSITIONp[COMPOSITION == 0] <- NA
    par(mar=c(4,10,1.5,2))
    plot(-log10(COMPOSITIONp[1,oo]), 1:ncol(COMPOSITIONp), col='#004add', pch=16,
        xlim=c(0,6), xlab='enrichment -log10 pval', yaxt='n',
        frame.plot=FALSE, ylab='')
    points(-log10(COMPOSITIONp[2,oo]), 1:ncol(COMPOSITIONp), col='#d5aff3', pch=16)
    points(-log10(COMPOSITIONp[3,oo]), 1:ncol(COMPOSITIONp), col='red', pch=16)
    abline(v= -log10(0.05), lty=2)
    dev.off()
  }

  n_tissues <- length(tissues)
  cat("DAM-bearing genes enriched for cancer type specific TSGs for ",
      length(which(COMPOSITIONp[1, ] < 0.05)), " cancer types (",
      round(100 * length(which(COMPOSITIONp[1, ] < 0.05)) / n_tissues, 2), "%)\n", sep = "")

  cat("DAM-bearing genes enriched for cancer type specific ambiguous drivers for ",
      length(which(COMPOSITIONp[2, ] < 0.05)), " cancer types (",
      round(100 * length(which(COMPOSITIONp[2, ] < 0.05)) / n_tissues, 2), "%)\n", sep = "")

  cat("DAM-bearing genes enriched for cancer type specific OG drivers for ",
      length(which(COMPOSITIONp[3, ] < 0.05)), " cancer types (",
      round(100 * length(which(COMPOSITIONp[3, ] < 0.05)) / n_tissues, 2), "%)\n", sep = "")

  return(list(COMPOSITION = COMPOSITION,
              COMPOSITIONp = COMPOSITIONp,
              CT_DAM_Bearing = CT_DAM_Bearing))
}
