# -------- DRUG RESPONSE ANALYSIS -------- #

library(stringr)

# Performs a hypergeometric test.
my.hypTest <- function(x,k,n,N) {  # mutated in the top k, top k most dependent, number mutated cell lines, tot cell lines
  
  PVALS <- phyper(x-1,n,N-n,k,lower.tail=FALSE)
  
  return(PVALS)
}

# Identifies drug response associations for genes within a tissue.
# For genes showing significant dependency, retrieves drugs
# targeting the gene and tests whether mutant cell lines are more sensitive.
# Computes rank-based enrichment, hypergeometric and empirical p-values,
# defines significant gene–drug associations, generates diagnostic plots,
# and saves detailed screening results per gene.
# 
# Diagnostic plots include:
# Ranked ln(IC50) values across cell lines with mutant lines highlighted
# and drug concentration thresholds indicated.
analyze_drugs <- function(RESTOT, CMP_annot, drugTargetInfo, gdscAll, DRplotsPath, ctiss, RR_th=1.71, dep_threshold=-0.5, hypTest_p=0.2, 
                 empPval=0.2, drug_RR=1.5, drug_HG_pval=0.2, drug_EMP_pval=0.2, n_rantrials = 1000, display=TRUE, tissue_idx=NULL, ntiss=NULL) {
  
  # select tissue-specific cell lines
  clTiss <- CMP_annot$model_id[CMP_annot$cancer_type == ctiss]
  
  # retains strongest DAMs
  RESTOT <- RESTOT[which(RESTOT$rank_ratio < RR_th & 
                            RESTOT$medFitEff < dep_threshold & 
                            RESTOT$hypTest_p < hypTest_p & 
                            RESTOT$empPval < empPval),]
  
  RES <- lapply(1:nrow(RESTOT), function(x) {
    
    # extract gene, mutation, mutant cell lines
    target <- RESTOT$GENE[x]
    variant <- RESTOT$var[x]
    cellLines <- unlist(str_split(RESTOT$ps_cl[x],', '))

    # identifies drugs targeting that gene
    ids <- which(drugTargetInfo$Gene.Target==target)
    drug_ids <- drugTargetInfo$Drug.ID[ids]
    drug_names <- drugTargetInfo$Name[ids]
    
    drug_ids <- sort(unique(gdscAll$uDRUG_ID[is.element(gdscAll$DRUG_ID,drug_ids)]))
    drug_names <- gdscAll$DRUG_NAME[match(drug_ids,gdscAll$uDRUG_ID)]
    
    if (length(ids)>0) {
      
      cat("\r", strrep(" ", 120), "\r")
      cat(sprintf("SAM analysis: processing tissue %3d/%3d (%-15s)", tissue_idx, ntiss, ctiss))
      flush.console()
      
      # subset drug response data 
      drug_response_data <- gdscAll[which(is.element(gdscAll$uDRUG_ID,drug_ids) & 
                                           is.element(gdscAll$SANGER_MODEL_ID,clTiss)),]

      # extracts unique drug IDs and model IDs
      udid <- unique(drug_response_data$uDRUG_ID)
      ucl <- unique(drug_response_data$SANGER_MODEL_ID)
      
      # initialize matrices 
      zscores <- matrix(NA,nrow = length(udid),
                        ncol = length(ucl),
                        dimnames = list(udid,ucl))
      
      lnIC50 <- matrix(NA,nrow = length(udid),
                       ncol = length(ucl),
                       dimnames = list(udid,ucl))
      
      # fill matrices 
      if(nrow(drug_response_data)>0) {
        for (i in 1:nrow(drug_response_data)){
          zscores[as.character(drug_response_data$uDRUG_ID)[i],
                  drug_response_data$SANGER_MODEL_ID[i]] <- drug_response_data$Z_SCORE[i]
          
          lnIC50[as.character(drug_response_data$uDRUG_ID)[i],
                 drug_response_data$SANGER_MODEL_ID[i]] <- drug_response_data$LN_IC50[i]
        }
      }
      
      # create drug annotation table
      additionalInfos <- as.data.frame(matrix(NA,nrow=length(udid),ncol=6))
      colnames(additionalInfos) <- c('screen','drug_id','drug_name','put_target','min_conc','max_conc')
      rownames(additionalInfos) <- sort(udid)
      
      # fill annotation info
      ii <- match(rownames(additionalInfos),drug_response_data$uDRUG_ID)
      additionalInfos[,'screen'] <- drug_response_data$DATASET[ii]
      additionalInfos[,'drug_id'] <- as.character(drug_response_data$DRUG_ID[ii])
      additionalInfos[,'drug_name'] <- drug_response_data$DRUG_NAME[ii]
      additionalInfos[,'put_target'] <- drug_response_data$PUTATIVE_TARGET[ii]
      additionalInfos[,'min_conc'] <- drug_response_data$MIN_CONC[ii]
      additionalInfos[,'max_conc'] <- drug_response_data$MAX_CONC[ii]
      
      # keep only mutant lines with IC50 data
      MutantCellLines <- intersect(cellLines,colnames(lnIC50)) 
      
      if(length(MutantCellLines)>0) {
        
        ndrugs <- nrow(additionalInfos)
        
        # loop over drugs
        SAMscores <- do.call('rbind', lapply(1:ndrugs, function(y) {
          
          # sort IC50 values
          drpattern <- sort(lnIC50[rownames(additionalInfos)[y],])
          # remove NAs
          drpattern <- drpattern[!is.na(drpattern)]
          
          # median mutant sensitivity
          medianLnIC50ofMutant <- median(drpattern[MutantCellLines],na.rm = TRUE)
          
          # find ranks of mutant cell lines and compute RR
          hits <- match(MutantCellLines,names(drpattern)) 
          hits <- hits[!is.na(hits)]
          nhits <- length(hits)
          rankRatio <- sum(hits)/(nhits*(nhits+1)/2)
          
          # hypergeometric test to compute enrichment p-value
          k <- max(hits)
          x <- length(hits)
          n <- length(hits)
          N <- length(which(!is.na(drpattern)))
          HG_pval <- my.hypTest(x,k,n,N)
          
          # permutation test
          NT <- n_rantrials
          randRankRatio <- rep(NA,NT)
          for (rt in 1:NT) {
            hits <- match(MutantCellLines,names(drpattern[sample(length(drpattern))]))
            hits <- hits[!is.na(hits)]
            nhits <- length(hits)
            randRankRatio[rt] <- sum(hits)/(nhits*(nhits+1)/2)
          }
          
          # empirical p-value
          EMP_pval <- (1+(length(which(randRankRatio<=rankRatio))))/(1+NT)
          
          res <- c(medianLnIC50ofMutant,rankRatio,HG_pval,EMP_pval)
        }))
        
        colnames(SAMscores) <- c('medLn50_mutCLs','rankRatio','HG_pval','EMP_pval')
        
        additionalInfos <- cbind(additionalInfos,SAMscores)  
        
        # defines a significant gene-drug association 
        validated <- (additionalInfos$medLn50_mutCLs < log(additionalInfos$max_conc)) & 
          (additionalInfos$rankRatio <= drug_RR) &
          (additionalInfos$HG_pval < drug_HG_pval) &
          (additionalInfos$EMP_pval < drug_EMP_pval)
        
        additionalInfos <- cbind(additionalInfos,validated)
        
        # keep validated drugs
        SAMs <- additionalInfos[which(additionalInfos$validated),]
        
        # plotting
        if(display && nrow(SAMs)>0) {
          lapply(1:nrow(SAMs),function(p) {
            
            allPattern <- lnIC50[rownames(SAMs)[p],]
            bg <- rep(rgb(0,0,255,alpha = 110,maxColorValue = 255),length(allPattern))
            names(bg) <- names(allPattern)
            bg[MutantCellLines] <- 'red'
            
            dname <- SAMs$drug_name[p]
            dname <- str_replace_all(dname,'/','|')
            
            pdf_file <- paste0(target, " - ", SAMs$screen[p], "_", SAMs$drug_id[p], "_", dname, ".pdf")
            pdf_file <- gsub("\\|", "", pdf_file)
            pdf(file.path(DRplotsPath, ctiss, pdf_file), width = 8.85, height = 4.20)
            
            par(mfrow=c(1,2))
            plot(sort(allPattern),
                 col=bg[order(allPattern)],xlab='cell lines',ylab='ln IC50',pch=16,
                 main = paste(dname,' [',target,']\nrankRatio = ',round(SAMs$rankRatio[p],2),
                              '\nHGp = ',formatC(SAMs$HG_pval[p], format = "e", digits = 2),
                              '\nEMPp =',formatC(SAMs$EMP_pval[p], format = "e", digits = 2),sep=''),
                 cex.main = 0.8)
            
            abline(h=log(SAMs$max_conc[p]),col='gray',lty=2)
            abline(h=SAMs$medLn50_mutCLs[p],col='pink',lty=1)
            
            plot(0,0,frame.plot = FALSE,axes = FALSE,col=NA,xlab='',ylab='')
            legend('left',pch=16,col=c('red',"#0000FF6E"),
                   legend=c(variant,'others'),
                   bg = 'white',title=paste(target,'status'),cex=0.8)
            dev.off()
            
          })
        }
        
        # creating the SCREENdata object with full drug-response data for later reuse
        SCREENdata <- list(screenInfo=additionalInfos,
                           lnIC50=lnIC50,
                           Zscores=zscores,
                           MutantCellLines=MutantCellLines)
        
        save_file <- paste0(target, " _ ", paste(gsub("\\?|\\*|!|>|\\|", "", variant), collapse = "AND"), "_screenRes.RData")

        list(
          save_file = save_file,
          SCREENdata = SCREENdata
        )

      }
      
    }
  })
  
  RES <- Filter(Negate(is.null), RES)
  return(RES)
}

