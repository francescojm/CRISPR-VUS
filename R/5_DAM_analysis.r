### -------- CRISPR TISSUE-SPECIFIC ANALYSIS -------- ##

suppressPackageStartupMessages(library(sets))
library(binaryLogic) 
suppressPackageStartupMessages(library(tidyverse))
library(CoRe) #installed from github

# Extracts the set of unique, non-ambiguous protein mutations observed in a
# given gene across cell lines. Returns a data frame containing one row per
# distinct mutation with associated annotation fields
variantSpectrum <- function(cl_var,gene) { 
  
  # extract variants in that gene
  tmp <- cl_var[cl_var$gene_symbol_2023==gene,c('protein_mutation', 'cdna_mutation', 'gene_symbol_2023')] 
  
  # list of mutations
  aa <- sort(unique(tmp$protein_mutation))
  aa <- setdiff(aa, c("-", "p.?"))
  
  aa <- tmp[match(aa,tmp$protein_mutation),]
  
  # df with mutations in that gene
  return(aa)
}

# Returns cell lines carrying at least one of the specified protein variants in the given gene
positives <- function(cl_var,gene,variants) {
  return(sort(unique(cl_var$model_id[cl_var$gene_symbol_2023==gene & is.element(cl_var$protein_mutation,variants)])))
}

# Tests whether a set of mutations in a gene is associated with increased CRISPR dependency
# Computes:
# - Rank ratio of mutated cell lines in the dependency ranking
# - Median fitness effect among mutated cell lines
# - Jaccard index between mutated and essential cell lines
# - Expression statistics
#
# Plotting:
# Generates a two-panel diagnostic plot for the tested gene–mutation set:
# (1) Ranked CRISPR fitness effects across cell lines with mutated lines highlighted
# (2) Histogram of fitness effect distribution with mutated lines overlaid
clasTest <- function(ess_scores,bess_scores,cl_var,gene,vs,vs_cds,produce_plots=TRUE,save=NULL,dep_threshold) {
  
  ps_cl <- positives(cl_var,gene,variants=vs)    # cell lines with mutations
  ps_cl <- intersect(ps_cl,colnames(ess_scores)) # and with CRISPR data
  
  # finds ranking of mutated cell lines (in the list of CRISPR effect values of a gene)
  FCp <- unlist(ess_scores[gene,]) 
  hits <- match(ps_cl,names(FCp)[order(FCp)])    # increasing order
  
  # number of mutated cell lines
  n <- length(hits)
  
  if (sum(hits)>0) {
    dependent_cls <- names(which(bess_scores[gene,]>0)) # cell lines where the gene is essential
    
    JI_essMut <- length(intersect(dependent_cls,ps_cl))/length(union(dependent_cls,ps_cl))
    rankRatio <- sum(hits)/((n*(n+1))/2)
    medEff <- median(FCp[ps_cl]) # median fitness of the mutated cell lines
    
    if(length(ps_cl) != length(FCp)) {
      # most dependent cell line, not mutated (baseline)
      mostDepNegative <- sort(FCp[setdiff(names(FCp),ps_cl)])[1] 
    } else { # perfect overlap
      mostDepNegative <- NA
    }
    
    if (is.element(gene,rownames(basal_exp))) {
      expPattern <- basal_exp[gene,intersect(colnames(basal_exp),names(FCp))]
      avgEXP_fpkm_in_positive_cl <- mean(expPattern[ps_cl])
      
      if(sum(!is.na(expPattern)) >= 10) {
        cdf <- ecdf(expPattern)
        perc_basal_exp_of_Avg_positive_cl <- round(100*cdf(mean(expPattern[ps_cl]))) 
      } else {
        
        DD <- sort(expPattern)
        testval <- mean(expPattern[ps_cl])
        # percentile of expression wrt all cell lines
        perc_basal_exp_of_Avg_positive_cl <- 100*max(which(DD<=testval))/length(DD)
      }
    } else {
      avgEXP_fpkm_in_positive_cl <- NA
      perc_basal_exp_of_Avg_positive_cl <- NA
    }
    
    RES <- data.frame(GENE=gene,var=c(paste(vs,collapse=' | ')),var_cds=c(paste(vs_cds,collapse=' | ')),
      npos=length(ps_cl),                                        # number of mutated cell lines for this gene and variant set
      medFitEff=medEff,                                          # median CRISPR fitness effect among mutated cell lines only
      leap=medEff-mostDepNegative,                               # how much the mutation explains the dependency
      totcl=length(FCp),                                         # total number of cell lines with CRISPR data for this gene
      rank_ratio = rankRatio,                                    # rank ratio
      matching=JI_essMut,                                        # Jaccard index
      ps_cl=paste(ps_cl,collapse=', '),                          # list of cell lines carrying at least one of these variants
      highest_rank=max(hits),                                    # the worst (largest) rank among the mutated cell lines
      avgBasalEXP_fpkm_in_ps_cl=avgEXP_fpkm_in_positive_cl,      # average basal expression of this gene in mutated cell lines
      percBasalEXP_of_ps_cl=perc_basal_exp_of_Avg_positive_cl,   # percentile of the average expression in mutated cell lines, compared to all cell lines
      stringsAsFactors = FALSE)
    
    # plotting
    if (produce_plots) {
      FCp <- sort(FCp)
      
      pcpatt <- rep(16,length(FCp))
      pcpatt[which(bess_scores[gene,names(FCp)]>0)] <- 18
      
      bg <- rep(rgb(0,0,255,alpha = 110,maxColorValue = 255),length(FCp))
      names(bg) <- names(FCp)
      
      bg[ps_cl] <- 'red'
      
    # if(rankRatio==1 & JI_essMut>0 & medEff< -1){
        if(length(save)>0) {
          pdf(paste(save,'/',gene,'_',paste(gsub("\\?|\\*|!|>", "", vs),collapse='AND'),'.pdf',sep=''),7.60,7.20)
        }
        
        par(mfrow=c(2,1))
        plot(FCp,col=bg,pch=pcpatt,xlab='cell lines',
             ylab=paste(gene,'scaled fitness effect'),
             main=paste('Rank ratio =',format(sum(hits)/((n*(n+1))/2),digits=3),', mut/ess match = ',
                                                                 format(100*JI_essMut,digits=3),'%'))
        abline(h= -1,col='gray')
        abline(h= dep_threshold,col='gray',lty=2)
        legend('bottomright',pch=c(15,18,16),col=c('red','gray','gray'),
               legend = c(paste(vs,collapse=' | '),'significant fitness effect','non significant fitness effect'),bg = 'white')
        
        hist(FCp,100,col=rgb(0,0,255,alpha = 110,maxColorValue = 255),border=NA,xlab=paste(gene,'fitness effect'),main='')
        
        points(FCp[ps_cl],rep(0,length(ps_cl)),col='red',pch=15,cex=1.5)
        abline(v= -1,col='gray')
        abline(v= dep_threshold,col='gray',lty=2)
        
        if(length(save)>0) {
          dev.off()
        } 
      }    
    # }
    
  }
  
  return(RES)
  
}

# Identifies the subset of mutations within a gene that maximizes the CRISPR dependency signal.
# Evaluates all non-empty subsets of variants, selects the combination that
# optimizes rank ratio (subject to effect size threshold).
Optimal_clasTest <- function(ess_scores,bess_scores,cl_var,gene,vs,vs_cds,produce_plots=TRUE,save=NULL,RR_th,dep_threshold) {
  
  if(length(vs)>1) {
    # generate subsets
    vs_subset_idx <- do.call(rbind,as.binary(1:(2^length(vs)-1),n = length(vs)))+0  # k mutations, 2^k - 1 non empty subsets
    
    # loop over all subsets 
    res <- do.call(rbind,lapply(1:nrow(vs_subset_idx),function(x) {         
      
      # select mutations included in this subset
      curr_var <- vs[which(vs_subset_idx[x,]>0)]
      curr_var_cds <- vs_cds[which(vs_subset_idx[x,]>0)]

      # test association
      clasTest(ess_scores = ess_scores,
               bess_scores = bess_scores,
               cl_var = cl_var,
               gene = gene,
               vs = curr_var,
               vs_cds = curr_var_cds,
               produce_plots = FALSE,
               dep_threshold = dep_threshold)
    }))
    
    # res$rank_ratio[res$rank_ratio<RR_th]<-1 #serviva perché altrimenti si selezionano sempre varianti presenti in una sola linea
#   if(length(which(res$rank_ratio<RR_th))==0)
    if(length(which(res$rank_ratio<RR_th & res$medFitEff< dep_threshold))==0) {       # no subset passes thresholds
      res <- res[which(res$rank_ratio==min(res$rank_ratio)),]               # smallest RR
    } else {
      res <- res[which(res$rank_ratio<RR_th & res$medFitEff< dep_threshold),]         # keep only subsets passing both thresholds      
    }
    
    # if multiple valid subsets remain, choose the subset containing the largest number of mutations
    if(nrow(res)>1) {
      ind <- which.max(unlist(lapply(res$var, function(x){length(unlist(str_split(x," \\| ")))})))
      res <- res[ind,]
    }
  } else {
    vs_subset_idx <- 1
    
    res <-
      clasTest(ess_scores = ess_scores,
               bess_scores = bess_scores,
               cl_var = cl_var,
               gene = gene,
               vs = vs,
               vs_cds = vs_cds,
               produce_plots = FALSE,
               dep_threshold = dep_threshold)
    }
  
  if(res$rank_ratio < RR_th & res$medFitEff < dep_threshold) {
    # reconstruct mutation vectors
    curr_var <- setdiff(unlist(str_split(res$var,'[ | ]')),'')
    curr_var_cds <- setdiff(unlist(str_split(res$var_cds,'[ | ]')),'')
    
    clasTest(ess_scores = ess_scores,
             bess_scores = bess_scores,
             cl_var = cl_var,
             gene = gene,
             vs = curr_var,
             vs_cds = curr_var_cds,
             produce_plots = produce_plots,
             save = save,
             dep_threshold = dep_threshold)
  }
  
  return(res) 
}

# Retrieves mutation frequencies for a specific CDS/AA mutation from COSMIC,
# parsing both targeted and WGS data across tumor types.
# Returns a data frame with per-disease mutation counts, tested samples,
# percentages, and annotation details, sorted by decreasing frequency.
getCOSMICfreqs <- function(COSMICf,GENE,Mutation_CDS,Mutation_AA) {
  # filters the COSMIC table to keep only the row matching gene and mutations
  RES <- COSMICf[which(COSMICf$GENE_NAME==GENE & COSMICf$`Mutation CDS`==Mutation_CDS & COSMICf$`Mutation AA`==Mutation_AA),]

  # extracts disease names, number of positives, and number of tested
  tmp <- unlist(str_split(RES$DISEASE,';'))
  tmp <- str_split(tmp,'=')
  disease <- unlist(lapply(tmp,function(x){x[1]}))
  n_positives <- unlist(lapply(str_split(unlist(lapply(tmp,function(x){x[2]})),'/'),function(x){as.numeric(x[1])}))
  n_tested <- unlist(lapply(str_split(unlist(lapply(tmp,function(x){x[2]})),'/'),function(x){as.numeric(x[2])}))
  
  # computes percentages
  percs <- round(100*n_positives/n_tested,2)
  
  # same applied to WGS
  tmp_wgs <- unlist(str_split(RES$WGS_DISEASE,';'))
  tmp_wgs <- str_split(tmp_wgs,'=')
  disease_wgs <- unlist(lapply(tmp_wgs,function(x){x[1]}))
  n_positives_wgs <- unlist(lapply(str_split(unlist(lapply(tmp_wgs,function(x){x[2]})),'/'),function(x){as.numeric(x[1])}))
  n_tested_wgs <- unlist(lapply(str_split(unlist(lapply(tmp_wgs,function(x){x[2]})),'/'),function(x){as.numeric(x[2])}))
  
  percs_wgs <- round(100*n_positives_wgs/n_tested_wgs,2)
  
  # merge targeted + WGS
  if(!is.na(tmp_wgs[1])) {
    disease <- c(disease,disease_wgs)
    n_positives <- c(n_positives,n_positives_wgs)
    n_tested <- c(n_tested,n_tested_wgs)
    percs <- c(percs,percs_wgs)
    
  }
  
  ncas <- length(percs)
  
  RES <- data.frame(
    GENE=rep(GENE,ncas),
    ROLE = rep(RES$ONC_TSG, ncas),
    CDS = rep(RES$`Mutation CDS`, ncas),
    AA = rep(RES$`Mutation AA`, ncas),
    desc_CDS = rep(RES$`Mutation Description CDS`, ncas),
    desc_AA = rep(RES$`Mutation Description AA`, ncas),
    N_POS_TOTAL = rep(RES$COSMIC_SAMPLE_MUTATED,ncas),
    N_TESTED_TOTAL = rep(RES$COSMIC_SAMPLE_TESTED,ncas),
    DISEASE = disease,
    nMutant = n_positives,
    nTested = n_tested,
    perc = percs, stringsAsFactors = FALSE)
  
  RES <- RES[order(RES$perc,decreasing=TRUE),]
  
  return(RES)
  
}

# Tests all mutated genes within a tissue for mutation-specific dependency.
# For genes showing dependency signals, performs permutation testing to 
# compute empirical p-values and hypergeometric enrichment tests. 
# Adjusts p-values (FDR) and returns a ranked results data frame.
multipleGeneClasTests <- function(ess_scores, bess_scores, cl_var, genes, path, RR_th, dep_threshold, n_rantrials, produce_plots=TRUE, tissue_idx=NULL, ntiss=NULL, tissue_name=NULL) {
  
  nn <- ncol(ess_scores)
  
  genes <- intersect(genes,cl_var$gene_symbol_2023)
  ngenes <- length(genes) # genes to test (with mutations)

  resTOT <- NULL
  for (i in 1:ngenes) {
    cat("\r", strrep(" ", 120), "\r")
    cat(sprintf("CRISPR analysis: processing tissue %3d/%3d (%-15s) | gene %4d/%4d (%-15s)", tissue_idx, ntiss, tissue_name, i, ngenes, genes[i]))
    flush.console()

    # extracts mutations in gene
    vs <- variantSpectrum(cl_var = cl_var,gene = genes[i]) 
    vs_cds <- vs$cdna_mutation
    vs <- vs$protein_mutation
    
    # run optimal classification test (find best mutation combination and compute statistics)
    curRES <- Optimal_clasTest(ess_scores = ess_scores, 
                             bess_scores = bess_scores,
                             cl_var = cl_var,
                             gene = genes[i],
                             vs = vs,
                             vs_cds = vs_cds,
                             produce_plots = produce_plots,
                             save = path,
                             RR_th = RR_th,
                             dep_threshold = dep_threshold)
    
    if(curRES$rank_ratio<RR_th & curRES$medFitEff < dep_threshold) { 
      
      # prepare permutation test
      raness_scores <- ess_scores
      ranbess_scores <- bess_scores
      ranRANKRATIO <- rep(0,n_rantrials)

      # loop of permutations
      for (RR in 1:n_rantrials) { 
        ransamp <- sample(1:nn)   # generates permutation of cell lines
        # shuffle gene scores, null hypothesis 
        raness_scores[genes[i],] <- raness_scores[genes[i],ransamp]
        ranbess_scores[genes[i],] <- ranbess_scores[genes[i],ransamp]
        
        # computes RR under null
        ranRES <- Optimal_clasTest(ess_scores = raness_scores,
                                 bess_scores = ranbess_scores,
                                 cl_var = cl_var,
                                 gene = genes[i],
                                 vs = vs,
                                 vs_cds = vs_cds,
                                 save = path,
                                 produce_plots = FALSE,
                                 RR_th = RR_th,
                                 dep_threshold = dep_threshold)
        
        ranRANKRATIO[RR] <- ranRES$rank_ratio
        
      }
      
      # computes empirical permutation-based p-value
      empPval <- (1+(sum(ranRANKRATIO<=curRES$rank_ratio)))/(1+n_rantrials)
      # hypergeometric test
      hypTest_p <- my.hypTest(x=curRES$npos,k = curRES$highest_rank,n = curRES$npos,N = curRES$totcl)
      
    } else {
      empPval <- NaN
      hypTest_p <- NaN
    }
    
    # attach hypergeometric test p
    curRES <- cbind(curRES,hypTest_p)
    
    # attach empPval
    curRES <- cbind(curRES,empPval)
    
    resTOT <- rbind(resTOT,curRES)
  }

  # sort by RR
  if(nrow(resTOT)>0) {resTOT <- resTOT[order(resTOT$rank_ratio),]}  
  
  # FDR correction (hypergeometric)
  hypTest_FDR <- rep(NaN,ngenes)
  hypTest_FDR[which(!is.na(resTOT$hypTest_p))] <-
    p.adjust(resTOT$hypTest_p[which(!is.na(resTOT$hypTest_p))],'fdr')
  resTOT <- cbind(resTOT,hypTest_FDR)
  
  # FDR correction (empirical p-values)
  empP_FDR <- rep(NaN,ngenes)
  empP_FDR[which(!is.na(resTOT$empPval))] <- p.adjust(resTOT$empPval[which(!is.na(resTOT$empPval))],'fdr')
  resTOT <- cbind(resTOT,empP_FDR)
  
  return(resTOT)
}

# Performs a hypergeometric test to assess whether mutated cell lines are
# enriched among the top-k most dependent cell lines.
my.hypTest <- function(x,k,n,N) {  # mutated in the top k, top k most dependent, number mutated cell lines, tot cell lines
  
  PVALS <- phyper(x-1,n,N-n,k,lower.tail=FALSE)
  
  return(PVALS)
}

# Executes the full tissue-specific analysis pipeline:
# - Subsets CRISPR data and mutations to the selected tissue
# - Filters mutations (coding, non-silent, non-truncating, etc.)
# - Selects candidate genes (excluding highly mutated and core-essential genes)
# - Runs dependency-mutation association testing across genes
# - Annotates results with core fitness gene information
# - Saves result tables, tested variants, plots, and DAM-bearing genes
process_tissue <- function(ctiss, scaled_depFC, bdep, CMP_annot, cl_variants, CRISPRplotsPath, RR_th=1.71, dep_threshold=-0.5, 
n_rantrials=1000, produce_plots=TRUE, ADaM=NULL, Perc_AUC=NULL, t_idx=NULL, ntiss=NULL) {

    # select the tissue and filter the cell lines accordingly
    ts_depFC <- scaled_depFC[,CMP_annot$model_id[CMP_annot$cancer_type==ctiss]] # select cell lines in ctiss
    ts_bdep <- bdep[,CMP_annot$model_id[CMP_annot$cancer_type==ctiss]]          # bin fitness CRISPR tissue specific

    # retains only mutations observed in cell lines of ctiss
    ts_cl_variants <- cl_variants[which(is.element(cl_variants$model_id,colnames(ts_depFC))),]

    # excluding non-coding mutations
    ts_cl_variants <- ts_cl_variants[which(ts_cl_variants$coding),] 
    # excluding start_lost mutations
    ts_cl_variants <- ts_cl_variants[which(ts_cl_variants$effect!='start_lost'),]
    # excluding silent mutations
    ts_cl_variants <- ts_cl_variants[which(ts_cl_variants$effect!='silent'),]
    # excluding nonsense mutations 
    ts_cl_variants <- ts_cl_variants[which(ts_cl_variants$effect!='nonsense'),]
    # excluding stop_lost mutations 
    ts_cl_variants <- ts_cl_variants[which(ts_cl_variants$effect!='stop_lost'),]
    # excluding ess_splice mutations 
    ts_cl_variants <- ts_cl_variants[which(ts_cl_variants$effect!='ess_splice'),]

    # select the genes to analyze (i.e. those with cancer dependency data available)
    ts_cl_variants <- ts_cl_variants[which(is.element(ts_cl_variants$gene_symbol_2023,rownames(ts_depFC))),]

    genesToTest <- unique(ts_cl_variants$gene_symbol_2023)

    # counts the number of observed different mutations per gene
    vs_spec_cardinality <- unlist(lapply(lapply(1:length(genesToTest),
                                            function(x){
                                                dd <- variantSpectrum(cl_var = ts_cl_variants,gene = genesToTest[x])
                                                dd <- dd$protein_mutation}),'length'))

    names(vs_spec_cardinality) <- genesToTest

    # exclude highly mutated genes (e.g. TTN and Tumor suppressor genes) 
    genesToTest <- sort(names(which(vs_spec_cardinality <= 10 & vs_spec_cardinality > 0)))

    # exclude pan-cancer core fitness genes and common-essential genes
    genesToTest <- setdiff(genesToTest,ADaM)
    genesToTest <- setdiff(genesToTest,Perc_AUC)

    # run main function and generate a dataframe with DAMs (tissue specific)
    RESTOT <- multipleGeneClasTests(
      ess_scores = ts_depFC,
      bess_scores = ts_bdep,
      cl_var = ts_cl_variants,
      genes = genesToTest,
      path = file.path(CRISPRplotsPath,ctiss),
      RR_th = RR_th,
      dep_threshold = dep_threshold,
      n_rantrials = n_rantrials,
      produce_plots = produce_plots,
      tissue_idx = t_idx,
      ntiss = ntiss,
      tissue_name = ctiss
    )
    
    RESTOT <- data.frame(ctype=rep(ctiss,nrow(RESTOT)),RESTOT,stringsAsFactors = FALSE)

    # annotate genes with info on core fitness genes
    data(curated_BAGEL_essential)
    cfg <- CoRe.ADaM(ts_bdep,TruePositives = curated_BAGEL_essential,display = FALSE,verbose = FALSE)
    is_cts_cfg <- is.element(RESTOT$GENE,cfg)
    RESTOT <- cbind(RESTOT,is_cts_cfg)

    DAM_bearing_genes <- unique(RESTOT$GENE[RESTOT$rank_ratio < RR_th & RESTOT$medFitEff < dep_threshold])
    DAM_bearing_genes <- unique(DAM_bearing_genes)

    return(list(
      RESTOT = RESTOT,
      ts_cl_variants = ts_cl_variants,
      DAM_bearing_genes = DAM_bearing_genes
    ))

}