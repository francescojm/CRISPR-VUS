# -------- INTOGEN AND COSMIC PATIENTS ANALYSIS -------- #

library(readr)
library(stringr)

# Build an aggregated intOGen count table, assign each entry to a gene-protein variant label
create_intogen_counts_aggr <- function(loc_file, allDAMs_Positions) {

  # Build composite identifiers once
  loc_key <- paste(loc_file$gene_name,
    loc_file$chromosome,
    loc_file$position,
    loc_file$allele_string,
    sep = "\r")

  dam_key <- paste(allDAMs_Positions$gene_name,
    allDAMs_Positions$chromosome,
    allDAMs_Positions$position,
    allDAMs_Positions$allele_string,
    sep = "\r")

  # Match every loc_file row 
  idx <- match(loc_key, dam_key)

  # Allocate output
  var <- rep(NA_character_, nrow(loc_file))

  matched <- !is.na(idx)

  if (any(matched)) {

    gene_name <- allDAMs_Positions$gene_name[idx[matched]]
    protein_mutation <- allDAMs_Positions$protein_mutation[idx[matched]]

    # Remove leading "p."
    protein_mutation <- sub("^p\\.", "", protein_mutation)

    var[matched] <- paste(gene_name, protein_mutation, sep = "-")
  }

  # Count columns and append variant
  intogen_counts_aggr <- cbind(loc_file[, 2:(ncol(loc_file) - 7), drop = FALSE], var = var)

  return(intogen_counts_aggr)
}

# Summarize DAM recurrence across tissues and external cohorts 
# from intOGen and COSMIC
compute_summary_and_tissue_lists <- function(allDAMs, intogen_counts_aggr, COSMIC, COSMICprimarysite, ntiss) {
  
  vars <- paste(allDAMs$GENE, gsub("p.", "", allDAMs$var), sep = "-")
  COSMIC$MUTATION_AA <- gsub("p.", "", COSMIC$MUTATION_AA)

  # to uniform COSMIC nomenclature
  COSMIC$MUTATION_AA[grep("fs\\*", COSMIC$MUTATION_AA)] <-
    gsub("[[:upper:]]fs", "fs", COSMIC$MUTATION_AA[grep("fs\\*", COSMIC$MUTATION_AA)])

  # Harmonize deletion notation in DAMs
  vars[grep("del[[:upper:]]+$", vars)] <-
    gsub("del[[:upper:]]+$", "del", vars[grep("del[[:upper:]]+$", vars)])

  allDAMs$var <- vars

  # Work only on unique variants
  vars_unique <- unique(vars)

  # Summary matrix and tissue-specific occurrence containers
  summary_vars <- matrix(nrow = length(vars_unique), ncol = 10)
  rownames(summary_vars) <- vars_unique

  all_tiss_intogen <- vector("list", length(vars_unique))
  names(all_tiss_intogen) <- vars_unique

  all_tiss_cosmic <- vector("list", length(vars_unique))
  names(all_tiss_cosmic) <- vars_unique

  # Pre-compute indices to avoid repeatedly scanning the full data frames

  # allDAMs lookup 
  allDAMs_ind <- split(seq_len(nrow(allDAMs)), allDAMs$var)

  # intOGen lookup 
  intogen_ind <- split(seq_len(nrow(intogen_counts_aggr)), intogen_counts_aggr$var)

  # COSMIC lookup, build the same gene-mutation identifier used for DAMs 
  COSMIC_var <- paste(COSMIC$GENE_SYMBOL, COSMIC$MUTATION_AA, sep = "-")

  # Only index COSMIC variants that are present among the DAMs 
  keep_cosmic <- COSMIC_var %in% vars_unique 
  COSMIC_ind <- split(which(keep_cosmic), COSMIC_var[keep_cosmic]) 
  
  # Free temporary object 
  rm(COSMIC_var, keep_cosmic) 
  gc()

  # Process each variant independently
  n_vars <- length(vars_unique)

  for(i in seq_along(vars_unique)) {
    dg <- vars_unique[i] 
    cat("Variant", i, "/", n_vars, "\n")

    gene <- unlist(strsplit(dg, "-"))[[1]]
    mut <- unlist(strsplit(dg, "-"))[[2]]
    # Store gene and variant labels
    summary_vars[dg,1] <- gene
    summary_vars[dg,2] <- mut
    # tissues where it is a hit
    ind_dam <- allDAMs_ind[[dg]] 
    if(!is.null(ind_dam)) { 
        tissues_hit <- allDAMs$ctype[ind_dam] 
        summary_vars[dg, 3] <- length(tissues_hit) / ntiss 
        summary_vars[dg, 4] <- paste(tissues_hit, collapse = " | ") }
    
    # intOGen recurrence summary
    ind_dg <- intogen_ind[[dg]] 

    if(!is.null(ind_dg)) {
        summary_vars[dg, 7] <- sum(intogen_counts_aggr[ind_dg, -74])
        summary_vars[dg, 9] <- paste( colnames(intogen_counts_aggr[, -74])[ which(intogen_counts_aggr[ind_dg, -74] > 0) ], collapse = " | " )
        ct_counts <- c()
        for(ct_ind in which(intogen_counts_aggr[ind_dg, -74] > 0)) {
            ct_counts <- c( ct_counts, intogen_counts_aggr[ind_dg, -74][ct_ind] ) 
        }

        if(sum(intogen_counts_aggr[ind_dg, -74])) {
            names(ct_counts) <- colnames(intogen_counts_aggr[, -74])[ which(intogen_counts_aggr[ind_dg, -74] > 0) ]
            all_tiss_intogen[[dg]] <- ct_counts
        }
    }

    # COSMIC recurrence summary
    ind_sel <- COSMIC_ind[[dg]]
    if (is.null(ind_sel)) {
        summary_vars[dg, 10] <- ""
        summary_vars[dg, 8] <- 0
    } else {
        summary_vars[dg, 10] <- paste( unique(COSMICprimarysite[ind_sel]), collapse = " | " )
        summary_vars[dg, 8] <- length( unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]) ) 
        ct_counts <- c() 
        
        for(ct in unique(COSMICprimarysite[ind_sel])) {
            # Search only among rows belonging to the current variant 
            ind_tiss <- ind_sel[which(COSMICprimarysite[ind_sel] == ct)]
            ct_counts <- c( ct_counts, length( unique( COSMIC$COSMIC_SAMPLE_ID[ind_tiss])))
        }

        if(length(unique(COSMICprimarysite[ind_sel])) > 0) {
            names(ct_counts) <- unique(COSMICprimarysite[ind_sel])
            all_tiss_cosmic[[dg]] <- ct_counts
        }
    }
  }

  # Convert the summary matrix into a clean data frame
  summary_vars <- data.frame(
    Gene = summary_vars[, 1],
    Var = summary_vars[, 2],
    Perc_tiss = as.numeric(summary_vars[, 3]),
    Tiss = summary_vars[, 4],
    Num_Intogen = as.numeric(summary_vars[, 7]),
    Num_COSMIC = as.numeric(summary_vars[, 8]),
    Type_Intogen = (summary_vars[, 9]),
    Type_COSMIC = (summary_vars[, 10])
  )
  
  rownames(summary_vars) <- vars_unique
  
  # Replace missing intOGen values with zeros / empty strings
  summary_vars$Num_Intogen[is.na(summary_vars$Num_Intogen)] <- 0
  summary_vars$Type_Intogen[is.na(summary_vars$Type_Intogen)] <- ""

  return(list(
    summary_vars = summary_vars,
    all_tiss_intogen = all_tiss_intogen,
    all_tiss_cosmic = all_tiss_cosmic,
    allDAMs = allDAMs,
    COSMIC = COSMIC
  ))
}

# Reformat and sort the aggregated intOGen count matrix for downstream comparison
finalize_intogen_counts_aggr <- function(intogen_counts_aggr) {

  # Remove entries that did not match any DAM
  intogen_counts_aggr <- intogen_counts_aggr[!is.na(intogen_counts_aggr$var), , drop = FALSE]

  rownames(intogen_counts_aggr) <- str_replace(intogen_counts_aggr$var, "-", " ")
  intogen_counts_aggr <- intogen_counts_aggr[, 1:(ncol(intogen_counts_aggr) - 1)]
  intogen_counts_aggr <- intogen_counts_aggr[order(rowSums(intogen_counts_aggr), decreasing = TRUE), ]
  return(intogen_counts_aggr)
}

# Convert tissue-specific COSMIC counts into a comparable aggregated count matrix
create_cosmic_counts_aggr <- function(all_tiss_cosmic) {
  
  # Collect all unique COSMIC tissue labels present across variants
  cosmic_c_type <- sort(unique(unlist(lapply(all_tiss_cosmic,function(x){names(x)}))))
  
  # Initialize a matrix of zeros: rows = variants, columns = cancer types
  cosmic_counts_aggr <- matrix(
    0,
    length(all_tiss_cosmic),
    length(cosmic_c_type),
    dimnames = list(names(all_tiss_cosmic), cosmic_c_type)
  )
  
  # Fill the matrix with tissue-specific COSMIC counts
  for (i in seq_along(all_tiss_cosmic)) {
    rn <- names(all_tiss_cosmic)[i]
    cnms <- names(all_tiss_cosmic[[i]])
    if (length(cnms) > 0) {
      cosmic_counts_aggr[rn, cnms] <- all_tiss_cosmic[[i]]
    }
  }
  
  # Keep only variants observed at least once in COSMIC
  cosmic_counts_aggr <- cosmic_counts_aggr[rowSums(cosmic_counts_aggr) > 0, ]

  rownames(cosmic_counts_aggr) <- str_replace(rownames(cosmic_counts_aggr), "-", " ")
  cosmic_counts_aggr <- cosmic_counts_aggr[order(rowSums(cosmic_counts_aggr), decreasing = TRUE), ]
  
  tmp <- rownames(cosmic_counts_aggr)
  tmp[which(tmp == "ARID1A Q1334del")] <- "ARID1A Q1334delQ"
  tmp[which(tmp == "EGR2 A309del")] <- "EGR2 A309delA"
  tmp[which(tmp == "GNAI2 K272del")] <- "GNAI2 K272delK"
  tmp[which(tmp == "ATAD5 K481del")] <- "ATAD5 K481delK"
  tmp[which(tmp == "SLX4 I1195del")] <- "SLX4 I1195delI"
  tmp[which(tmp == "E2F3 S73del")] <- "E2F3 S73delS"
  tmp[which(tmp == "PEX12 L123del")] <- "PEX12 L123delL"
  tmp[which(tmp == "EMSY I1112del")] <- "EMSY I1112delI"
  rownames(cosmic_counts_aggr) <- tmp

  return(cosmic_counts_aggr)
}

# Extend the intOGen-COSMIC cancer-type mapping by adding a row for unanalysed categories
create_cmp_mapping_v3 <- function(cmp_mapping, intogen_counts_aggr, cosmic_counts_aggr) {

  # Remove existing "Others" row if present
  cmp_mapping <- cmp_mapping[rownames(cmp_mapping) != "Others (not analysed)", ]
  
  # Collect all cancer types already mapped from intOGen
  allMappedIntogen_types <- sort(unique(unlist(lapply(str_split(cmp_mapping$intOGen.cancer.type,'[|]'),function(x){str_trim(x)}))))
  
  # Identify intOGen cancer types not included in the mapping table
  otherIntoGen_types <- paste(setdiff(colnames(intogen_counts_aggr),allMappedIntogen_types),collapse=' | ')
  
  # Collect all cancer types already mapped from COSMIC
  allMappedCosmic_types <- sort(unique(unlist(lapply(str_split(cmp_mapping$COSMIC.tissue,'[|]'),function(x){str_trim(x)}))))
  
  # Identify COSMIC tissue types not included in the mapping table
  otherCosmic_types <- paste(setdiff(colnames(cosmic_counts_aggr),allMappedCosmic_types),collapse=' | ')
  
  # Append an "Others" row to capture all unmapped categories
  cmp_mapping <- rbind(cmp_mapping, c(otherIntoGen_types, otherCosmic_types))
  rownames(cmp_mapping)[nrow(cmp_mapping)] <- "Others (not analysed)"

  return(cmp_mapping)
}

# Build comparable intOGen and COSMIC matrices and aggregate counts 
# using the shared mapping table
compute_aggregated_counts_aggr_comp <- function(intogen_counts_aggr, cosmic_counts_aggr, cmp_mapping) {
  
  # Identify variants present only in one source
  intogenOnly_vars <- setdiff(rownames(intogen_counts_aggr), rownames(cosmic_counts_aggr))
  cosmicOnly_vars <- setdiff(rownames(cosmic_counts_aggr), rownames(intogen_counts_aggr))
  
  # Add zero-filled rows
  cosmic_counts_aggr_comp <- rbind(
    cosmic_counts_aggr,
    matrix(0, nrow = length(intogenOnly_vars), ncol = ncol(cosmic_counts_aggr),
           dimnames = list(intogenOnly_vars, colnames(cosmic_counts_aggr)))
  )

  intogen_counts_aggr_comp <- rbind(
    intogen_counts_aggr,
    matrix(0, nrow = length(cosmicOnly_vars), ncol = ncol(intogen_counts_aggr),
           dimnames = list(cosmicOnly_vars, colnames(intogen_counts_aggr)))
  )
  
  # Sort both matrices identically by row name
  cosmic_counts_aggr_comp <- cosmic_counts_aggr_comp[order(rownames(cosmic_counts_aggr_comp)), ]
  intogen_counts_aggr_comp <- intogen_counts_aggr_comp[order(rownames(intogen_counts_aggr_comp)), ]
  
  # Aggregate counts across matched intOGen/COSMIC tissue group
  aggregated_counts_aggr_comp <- do.call(cbind, lapply(rownames(cmp_mapping), function(x) {
    intogenTypes <- str_trim(unlist(str_split(cmp_mapping[x, 1], "[|]")))
    cosmicTypes <- str_trim(unlist(str_split(cmp_mapping[x, 2], "[|]")))

    intogenTypes <- intersect(intogenTypes, colnames(intogen_counts_aggr_comp))
    cosmicTypes <- intersect(cosmicTypes, colnames(cosmic_counts_aggr_comp))

    return(rowSums(cbind(intogen_counts_aggr_comp[, intogenTypes], cosmic_counts_aggr_comp[, cosmicTypes])))
  }))

  colnames(aggregated_counts_aggr_comp) <- rownames(cmp_mapping)
  
  # Rank variants by total aggregated recurrence
  aggregated_counts_aggr_comp <- aggregated_counts_aggr_comp[order(rowSums(aggregated_counts_aggr_comp), decreasing = TRUE), ]

  tmp <- rownames(aggregated_counts_aggr_comp)
  tmp[which(tmp == "ARID1A Q1334del")] <- "ARID1A Q1334delQ"
  tmp[which(tmp == "EGR2 A309del")] <- "EGR2 A309delA"
  tmp[which(tmp == "GNAI2 K272del")] <- "GNAI2 K272delK"
  tmp[which(tmp == "ATAD5 K481del")] <- "ATAD5 K481delK"
  tmp[which(tmp == "SLX4 I1195del")] <- "SLX4 I1195delI"
  tmp[which(tmp == "E2F3 S73del")] <- "E2F3 S73delS"
  tmp[which(tmp == "PEX12 L123del")] <- "PEX12 L123delL"
  tmp[which(tmp == "EMSY I1112del")] <- "EMSY I1112delI"
  rownames(aggregated_counts_aggr_comp) <- tmp

  return(aggregated_counts_aggr_comp)
}
