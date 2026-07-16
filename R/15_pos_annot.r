# -------- DAM VARIANT POSITION ANNOTATION -------- # 

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(biomaRt)
})

# Filters variant data to retain only those matching DAMs,
# removes duplicates, orders by gene, and formats a clean table with genomic
# and mutation information.
filter_and_match_DAMs <- function(cl_variants0, allDAMs) {
  
  cl_vars <- cl_variants0[!is.na(cl_variants0$chromosome), ]
  
  # Create unique variant identifiers (gene + mutation)
  variant_signature <- paste(cl_vars$gene_symbol, cl_vars$protein_mutation)
  variant_signature1 <- paste(allDAMs$GENE, allDAMs$var)
  
  # Keep only variants present in DAMs
  ii <- which(is.element(variant_signature,variant_signature1))
  variant_signature <- variant_signature[ii]
  cl_vars <- cl_vars[ii, ]
  
  # Keep unique signatures
  uvariant_signature <- unique(variant_signature)
  cl_vars <- cl_vars[match(uvariant_signature, variant_signature), ]
  
  # Sort by gene
  cl_vars <- cl_vars[order(cl_vars$gene_symbol), ]
  
  # Select/rename columns
  cl_vars_df <- as.data.frame(cl_vars)
  cl_vars_df <- cl_vars_df[,c(1,2,6,7,8,9,10,11,12)]
  colnames(cl_vars_df) <- c('gene_name','gene_id','protein_mutation','rna_mutation', 'cdna_mutation','chromosome','position','ref','alt')

  return(cl_vars_df)
}

# Annotates DAM variants with genomic coordinates (start/end), allele strings, and gene strand information. 
annotate_DAM_positions <- function(cl_vars) {
  
  dt <- as.data.table(cl_vars)
  
  # seq_region_name
  dt[, seq_region_name := str_replace(chromosome, "^chr", "")]
  
  # start/end positions (SNVs)
  dt[, `:=`(start = position, end = position)]
  
  # Function to extract allele string from HGVS notation
  parse_allele <- function(hgvs) {
    if (is.na(hgvs)) return(NA_character_)
    
    # SNV example: c.1799T>A -> T/A
    m <- str_match(hgvs, "c\\.[0-9_]+([ACGT])>([ACGT])")
    if (!is.na(m[1,1])) return(paste0(m[1,2], "/", m[1,3]))
    
    # Deletion: c.123_125delTAA -> TAA/-
    m <- str_match(hgvs, "c\\.([0-9_]+)del([ACGT]+)?")
    if (!is.na(m[1,1])) return(paste0(ifelse(is.na(m[1,3]), "", m[1,3]), "/-"))
    
    # Insertion: c.123_124insT -> -/T
    m <- str_match(hgvs, "c\\.([0-9_]+)ins([ACGT]+)")
    if (!is.na(m[1,1])) return(paste0("-/", m[1,3]))
    
    return(NA_character_)
  }
  
  # Build allele_string using cDNA or RNA mutation notation
  dt[, allele_string := ifelse(!is.na(cdna_mutation),
                               vapply(cdna_mutation, parse_allele, character(1)),
                               vapply(rna_mutation, parse_allele, character(1)))]
  
  # Strand annotation from Ensembl (GRCh38)
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = NULL)
  strand_map <- getBM(
    attributes = c("ensembl_gene_id", "strand"),
    filters = "ensembl_gene_id",
    values = unique(dt$gene_id),
    mart = mart
  )
  dt <- merge(dt, as.data.table(strand_map), by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
  
  return(as.data.frame(dt))
}
