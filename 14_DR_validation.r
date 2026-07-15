# -------- DR VALIDATION PIPELINE -------- #

library(stringr)

# Expands rows where a hit contains multiple variants in the same entry,
# so that each gene-variant pair is represented as a separate row.
decoupleMultipleHits <- function(hitTable) {
  vars <- hitTable$Hit

  # separate gene and variant columns
  genes <- str_split_fixed(vars,' ',2)
  vars <- genes[,2]
  genes <- genes[,1]
  
  # Identify entries containing multiple variants
  iimultiple <- grep(' | ',vars)
  # Build extended table with explicit gene and variant columns
  ehitTable <- cbind(hitTable[,1],genes,vars,hitTable[,c(1,3:ncol(hitTable))])
  # Keep rows with only one variant unchanged
  finalHit <- ehitTable[setdiff(1:nrow(hitTable),iimultiple),]
  
  # Expand rows with multiple variants into one row per variant
  finalHit <- rbind(finalHit,
                  
                  do.call('rbind',lapply(iimultiple,function(x) {
                    
                    indivar <- setdiff(strsplit(vars[x],' | ')[[1]],'|')
                    tmpHit <- NULL
                    for (i in 1:length(indivar)) {
                      tmpHit <- rbind(tmpHit,ehitTable[x,])
                    }
                    
                    tmpHit$vars <- indivar
                    
                    return(tmpHit)
                  })))
  
  colnames(finalHit)[1] <- 'cancer_type'
  return(finalHit)
}

# For each tissue/cancer type, loads screening results, filters significant hits,
# matches them to available drug-response files, and builds a combined table
# of drug-response validations across tissues.
build_dr_validations <- function(tissues, resultPath, DRplotsPath, RR_th=1.71, dep_threshold=-0.5) {
  
  allDRvalidations <- NULL
  tissueDRvalidations <- list()
  
  for (ctiss in tissues) {
    
    load(file.path(resultPath, paste0(ctiss, "_results.RData")))
    
    # Select hits passing effect size and significance thresholds
    hits <- RESTOT[
      RESTOT$medFitEff < dep_threshold &
      RESTOT$rank_ratio < RR_th &
      RESTOT$hypTest_p < 0.20 &
      RESTOT$empPval < 0.20, ]
    
    allTar <- hits$GENE
    allVar <- hits$var
    
    plot_path <- file.path(DRplotsPath, ctiss)
    allfnames <- list.files(pattern = "RData", path = plot_path)
    fnames <- unlist(lapply(str_split(allfnames,' _ '),function(x){x[1]}))
    
    # Keep only hits for which a corresponding result file exists
    allVar <- allVar[which(allTar %in% fnames)]
    allTar <- intersect(allTar, fnames)
    
    if (length(allTar) > 0) {
      
      rRES <- do.call(rbind, lapply(1:length(allTar), function(i) {
        
        x <- allTar[i]
        
        # matching result file for the current target gene
        current_fn <- allfnames[
          match(x, unlist(lapply(str_split(allfnames,' _ '), function(x){x[[1]][1]})))
        ]
        
        load(file.path(plot_path, current_fn))
        
        nnd <- nrow(SCREENdata$screenInfo)
        
        # Combine tissue, hit, and screen-level information
        RES <- cbind(
          rep(ctiss, nnd),
          rep(paste(x, allVar[i]), nnd),
          SCREENdata$screenInfo
        )
        
        colnames(RES)[c(1,2)] <- c('ctype','Hit')
        rownames(RES) <- NULL
        
        return(RES)
      }))
      
      tissueDRvalidations[[ctiss]] <- rRES
      # Append to global validation table
      allDRvalidations <- rbind(allDRvalidations, rRES)
    }
  }
  
  allDRvalidations <- allDRvalidations[!is.na(allDRvalidations$validated), ]
  
  return(list(
    allDRvalidations = allDRvalidations,
    tissueDRvalidations = tissueDRvalidations
  ))
}

# Summarises the number of cancer-type-specific hits, the subset that can be
# assessed using available drug-response data, and how many of those are validated.
summarise_dr_validations <- function(allDRvalidations, allHits) {
  
  # Unique cancer type + hit combinations with DR information
  sigs <- paste(allDRvalidations$ctype, allDRvalidations$Hit)
  
  # Total number of cancer-type-specific hits in the original hit table
  ncts_hits <- length(unique(paste(allHits$ctype, allHits$GENE, allHits$var)))
  
  nvalidable <- length(unique(sigs))
  
  # Count how many validation entries for each hit
  validated <- unlist(lapply(unique(sigs), function(x) {
    ii <- which(sigs == x)
    sum(allDRvalidations[ii, ]$validated)
  }))
  
  # Number of unique hits with at least one successful validation
  nvalidated <- length(unique(sigs)[which(validated > 0)])
  
  cat("Of the", ncts_hits, "cancer-type-specific hits (DAMs or DAMs combinations with optimal RankRatio, fitness effect and pvalues),", nvalidable,
      "involve a DAM-bearing gene that is druggable and targeted by a compound with available cancer-type matching drug-response data on GDSC.\n")
  cat("Of these,", nvalidated, "are validated.\n")
  
  return(list(
    ncts_hits = ncts_hits,
    nvalidable = nvalidable,
    nvalidated = nvalidated
  ))
}

# Builds a table of validated SAMs, expands multiple-hit
# entries into individual gene-variant pairs, and annotates each SAM according
# to whether it involves a known cancer driver and/or a cancer-type-specific driver.
# Also labels potentially repurposable drug targets.
build_annotated_SAMs <- function(allDRvalidations, intOGen_drivers, intOGen_drivers_s, ctypeMapping) {
  
  # Keep only validated drug-response associations
  allDRvalidations <- allDRvalidations[allDRvalidations$validated == TRUE, ]
  
  allSAMs <- decoupleMultipleHits(hitTable = allDRvalidations)
  allSAMs <- allSAMs[order(allSAMs$cancer_type), ]

  cat(paste("Encompassing", length(unique(paste(allSAMs$cancer_type,allSAMs$genes,allSAMs$vars))),
            "individual cancer-type-specific drug validated DAMs (SAMs)\n"))
  
  # Subset SAMs involving known IntOGen driver genes
  ii <- which(is.element(allSAMs$genes, intOGen_drivers_s))
  kallSAMs <- allSAMs[ii, ]
  
  cat(paste("Encompassing", length(unique(paste(kallSAMs$cancer_type,kallSAMs$genes,kallSAMs$vars))),
            "individual cancer-type-specific drug validated DAMs (SAMs) involving known cancer drivers\n"))
  
  nSAMs <- nrow(allSAMs)
  
  isAcancerDriver <- rep(NA, nSAMs)
  isA_ct_specific_Driver <- rep(NA, nSAMs)
  
  for (i in 1:nSAMs) {
    # General cancer driver annotation
    isAcancerDriver[i] <- is.element(allSAMs$genes[i],intOGen_drivers$SYMBOL)
    intTypes <- str_trim(unlist(str_split(ctypeMapping[allSAMs$cancer_type[i],1],'\\|')))
    # Cancer-type-specific driver annotation
    isA_ct_specific_Driver[i] <- is.element(allSAMs$genes[i],intOGen_drivers$SYMBOL[which(is.element(intOGen_drivers$CANCER_TYPE,intTypes))])
  }
  
  # A SAM is considered repurposable if the gene is not a known driver
  # or not a driver in that specific cancer type
  repurposableDrug <- (!isAcancerDriver | !isA_ct_specific_Driver)
  
  allSAMs <- cbind(allSAMs, isAcancerDriver, isA_ct_specific_Driver, repurposableDrug)
  
  return(allSAMs)
}

