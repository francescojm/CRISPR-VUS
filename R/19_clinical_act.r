# -------- PATIENT-LEVEL DAM ANALYSIS, CLINICAL ACTIONABILITY AND TRACTABILITY -------- #

library(stringr)
library(readr)
library(readxl)
library(data.table)
library(openxlsx)

# Load and preprocess CIViC clinical-evidence data together with the disease-mapping table
load_civic_data <- function(raw_dir, civic_file, civic_mapping_file) {
  
  civicdb_var <- readr::read_tsv(file.path(raw_dir, civic_file))
  
  # select only specific molecular profiles (exclude fusion, deletion, etc)
  civicdb_var_sel <- civicdb_var[-grep(
    "fusion|mutation|expression|deletion|duplication|methylation|Loss-of-function|amplification",
    civicdb_var$molecular_profile, ignore.case = TRUE
  ),]
  
  civicdb_var_sel$molecular_profile <- gsub(" ", "-", civicdb_var_sel$molecular_profile)
  
  civicdb_mapping <- read.csv(file.path(raw_dir, civic_mapping_file), sep=";")
  
  list(
    civicdb_var_sel = civicdb_var_sel,
    civicdb_mapping = civicdb_mapping
  )
}

# Load DAM annotations and derive the list of activating driver genes
load_dams <- function(intogen_drivers, intogen_drivers_s, dams_file, path_results) {
  
  # Extract activating driver genes from the intOGen driver annotation table
  act_drivers <- unique(intogen_drivers$SYMBOL[intogen_drivers$ROLE=="Act"])
  
  # load DAMs with functional impact and clinical relevance
  load(file=file.path(path_results, "_allDAMs_with_funcImpact_and_clinical_rel.RData"))
  allDAMs <- allDAMs_with_funcImpact_and_clinical_rel
  # Build compact variant identifiers: gene-ProteinChange
  vars <- paste(allDAMs$GENE, gsub("p.", "", allDAMs$var), sep="-")
  allDAMs$vars <- vars
  
  list(
    act_drivers = act_drivers,
    allDAMs = allDAMs
  )
}

# Clean and reduce COSMIC data to the subset relevant for DAM analysis
clean_cosmic_data <- function(COSMIC, COSMICprimarysite, allDAMs, act_drivers, civicdb_var_sel, CMP_COSMIC_ctype_mapping) {
  
  # Uniform COSMIC nomenclature
  COSMIC$MUTATION_AA <- gsub("p.", "", COSMIC$MUTATION_AA)
  COSMIC$MUTATION_AA[grep("fs\\*", COSMIC$MUTATION_AA)] <- 
    gsub("[[:upper:]]fs", "fs", COSMIC$MUTATION_AA[grep("fs\\*", COSMIC$MUTATION_AA)])
  
  # Harmonize DAM deletion notation
  vars <- paste(allDAMs$GENE, gsub("p.", "", allDAMs$var), sep="-")
  vars[grep("del[[:upper:]]+$", vars)] <- gsub("del[[:upper:]]+$", "del", vars[grep("del[[:upper:]]+$", vars)])
  allDAMs$var <- vars
  
  # Collect all mapped COSMIC cancer types from the cancer-type mapping table
  allCMP_cosmic_types <- unique(unlist(lapply(strsplit(CMP_COSMIC_ctype_mapping[,2],'\\| '),'str_trim')))

  genes_of_interests <- unique(c(allDAMs$GENE, act_drivers,
                                 unlist(lapply(strsplit(civicdb_var_sel$molecular_profile,'-'), function(x) x[1]))))
  
  # Remove duplicated COSMIC samples
  dupes <- duplicated(COSMIC$COSMIC_SAMPLE_ID)
  COSMIC_tissue_Count <- summary(as.factor(COSMICprimarysite[!dupes]))
  
  # Count how many COSMIC samples fall into each CMP cancer-type grouping
  CMP_ctype_in_COSMIC_Cunt<-unlist(lapply(rownames(CMP_COSMIC_ctype_mapping),function(x) {
    cosmictypes<-setdiff(unlist(strsplit(CMP_COSMIC_ctype_mapping[x,2],'\\ | ')),'|')
    sum(COSMIC_tissue_Count[cosmictypes])
  }))
  names(CMP_ctype_in_COSMIC_Cunt) <- rownames(CMP_COSMIC_ctype_mapping)
  
  # Filter COSMIC for genes of interest and relevant cancer types
  ii <- which(COSMIC$GENE_SYMBOL %in% genes_of_interests & COSMICprimarysite %in% allCMP_cosmic_types)
  COSMIC <- COSMIC[ii,]
  COSMICprimarysite <- COSMICprimarysite[ii]
  
  # Count unique COSMIC samples per tissue after filtering
  ucosmicSamples <- unique(COSMIC$COSMIC_SAMPLE_ID)
  nsamples_per_tissues_in_COSMIC <- summary(as.factor(COSMICprimarysite[match(ucosmicSamples, COSMIC$COSMIC_SAMPLE_ID)]))
  
  list(COSMIC=COSMIC,
       COSMICprimarysite=COSMICprimarysite,
       COSMIC_tissue_Count=COSMIC_tissue_Count,
       CMP_ctype_in_COSMIC_Cunt=CMP_ctype_in_COSMIC_Cunt,
       nsamples_per_tissues_in_COSMIC=nsamples_per_tissues_in_COSMIC)
}

# Identify COSMIC patients carrying CIViC variants.
# For each CIViC molecular profile, the function matches gene and variant in COSMIC,
# optionally restricting to disease-specific tumour types using the CIViC mapping table.
extract_civic_patients <- function(civicdb_var_sel, civicdb_mapping, COSMIC, COSMICprimarysite) {
  
  patients_civic <- c()
  
  for(i in 1:length(civicdb_var_sel$molecular_profile)) {
    gene <- unlist(strsplit(civicdb_var_sel$molecular_profile[i], "-"))[1]
    mut  <- unlist(strsplit(civicdb_var_sel$molecular_profile[i], "-"))[2]
    
    if(!is.na(civicdb_var_sel$disease[i])){
      if(civicdb_var_sel$disease[i]=="Cancer"){
        ind_sel <- which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut)
      } else {
        ct <- civicdb_mapping[match(civicdb_var_sel$disease[i], civicdb_mapping[,1]), 2]
        ind_sel <- which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut & COSMICprimarysite %in% ct)
      }
      patients_civic <- c(patients_civic, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
    }
  }
  
  return(patients_civic)
}

# Select COSMIC patients carrying DAMs in the corresponding matched tumour types.
# For each DAM, the function:
# - extracts gene, variant, and cancer type
# - translates the internal cancer type into COSMIC tumour labels
# - finds COSMIC samples with that gene/variant in the mapped tumour types
# - records patient ID, cancer type, variant, and gene
select_dam_patients <- function(allDAMs, CMP_COSMIC_ctype_mapping, COSMIC, COSMICprimarysite) {
  
  DAMbearigPatient_in_COSMIC <- do.call(rbind,lapply(1:nrow(allDAMs),function(i) {
    
    gene <- allDAMs$GENE[i]
    mut <- unlist(strsplit(allDAMs$vars[i], "-"))[[2]]
    # Retrieve COSMIC tumour types corresponding to the DAM cancer type
    ct <- unlist(strsplit(CMP_COSMIC_ctype_mapping[allDAMs$ctype[i],2], " \\| "))
    
    ind_sel <- which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut & COSMICprimarysite %in% ct)
    
    patients_sel <- unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel])
    ct_sel <- rep(allDAMs$ctype[i], length(patients_sel))
    var_sel <- rep(allDAMs$vars[i], length(patients_sel))
    gene_sel <- rep(allDAMs$GENE[i], length(patients_sel))
    
    cbind(patients_sel, ct_sel, var_sel, gene_sel)
  }))

  if (is.null(DAMbearigPatient_in_COSMIC)) {
    return(data.frame(
      patients_sel = character(0),
      ct_sel = character(0),
      var_sel = character(0),
      gene_sel = character(0)
    ))
  }
  
  return(DAMbearigPatient_in_COSMIC)
}

# Annotate each DAM-bearing patient with:
# - whether the patient already carries a CIViC actionable mutation
# - whether the DAM falls in a known driver-gene 
annotate_patients <- function(DAMbearigPatient_in_COSMIC, patients_civic, drivers) {
  
  DAMbearigPatient_in_COSMIC <- data.frame(DAMbearigPatient_in_COSMIC)

  DAMbearigPatient_in_COSMIC$has_already_an_actionable_mutation <-
    is.element(DAMbearigPatient_in_COSMIC$patients_sel, patients_civic)
  
  DAMbearigPatient_in_COSMIC$isTheDAM_in_an_unreported_DAMbgs <-
    is.element(DAMbearigPatient_in_COSMIC$gene_sel, drivers)
  
  return(DAMbearigPatient_in_COSMIC)
}

# Compute summary statistics 
compute_patient_statistics <- function(DAMbearigPatient_in_COSMIC) {

  if (nrow(DAMbearigPatient_in_COSMIC) == 0) {
    cat("No DAM-bearing patients found in COSMIC\n")
    return(list(
      alreadyCurable = 0,
      coveredWithCancDrivers = 0,
      coveredWithNovel = 0
    ))
  }
  
  # Patients already covered by known CIViC actionable variants
  alreadyCurable <- sum(DAMbearigPatient_in_COSMIC$has_already_an_actionable_mutation)

  # Patients not already actionable but carrying DAMs in known driver genes
  coveredWithCancDrivers <- length(which(!DAMbearigPatient_in_COSMIC$has_already_an_actionable_mutation &
                                          DAMbearigPatient_in_COSMIC$isTheDAM_in_an_unreported_DAMbgs))
  
  # Patients not already actionable and carrying DAMs outside known driver genes
  coveredWithNovel <- length(which(!DAMbearigPatient_in_COSMIC$has_already_an_actionable_mutation &
                                    !DAMbearigPatient_in_COSMIC$isTheDAM_in_an_unreported_DAMbgs))
  
  cat('of ', nrow(DAMbearigPatient_in_COSMIC), ' patients with DAMs in COSMIC, ', nrow(DAMbearigPatient_in_COSMIC)-alreadyCurable,
    ' (', round(100*(nrow(DAMbearigPatient_in_COSMIC)-alreadyCurable)/nrow(DAMbearigPatient_in_COSMIC),2),
    '%) do not have already a clinically actionable mutation\n', sep='')
  
  list(
    alreadyCurable = alreadyCurable,
    coveredWithCancDrivers = coveredWithCancDrivers,
    coveredWithNovel = coveredWithNovel
  )
}

# Plot a pie chart comparing DAM-bearing patients who already have a known clinically
# actionable mutation versus those currently lacking one
plot_actionable_pie <- function(DAMbearigPatient_in_COSMIC, alreadyCurable, figuresPath, produce_plots = TRUE) {
  
  # FIGURE 4B
  if (produce_plots && nrow(DAMbearigPatient_in_COSMIC) > 0) {
    out_file <- file.path(figuresPath, "_patients_with_DAMs_and_actionable_muts.pdf")
    pdf(out_file, 8, 8)
    pie(
      c(alreadyCurable, nrow(DAMbearigPatient_in_COSMIC)-alreadyCurable),
      main=paste(nrow(DAMbearigPatient_in_COSMIC),'Patients with DAMs'),
      labels = c('With co-occurring clinically actionable mutation',
                'Currently lacking clinically actionable mutations'),
      border=NA, col=c('darkgreen','green')
    )
    dev.off()
  }

}

# Plot the distribution of cancer types among patients who:
# - carry a DAM
# - do not already have a CIViC actionable variant
# - carry a DAM outside known driver genes
plot_cancer_type_distribution <- function(DAMbearigPatient_in_COSMIC, CL_colors_v, CMP_ctype_in_COSMIC_Cunt, figuresPath, produce_plots = TRUE) {

  # Keep patients lacking CIViC actionability and carrying DAMs outside known driver genes
  ii <- which(!DAMbearigPatient_in_COSMIC$`has already a CIVIC variant` &
              !DAMbearigPatient_in_COSMIC$`DAM in a known cancer driver`)

  # Patients lacking CIViC actionability, considering all DAM-bearing genes
  ii_all <- which(!DAMbearigPatient_in_COSMIC$`has already a CIVIC variant`)
  
  # Absolute counts by cancer type for unreported
  toPlot <- sort(summary(factor(DAMbearigPatient_in_COSMIC$`cancer type`[ii])), decreasing=TRUE)

  # Absolute counts by cancer type for all DAM-bearing patients
  toPlot_all <- sort(summary(factor(DAMbearigPatient_in_COSMIC$`cancer type`[ii_all])), decreasing = TRUE)
  
  if (produce_plots && length(toPlot) > 0) {
    
    # SUPPLEMENTARY FIGURE 9J
    pdf(file.path(figuresPath, "06_DAMs_in_patients_lacking_act_unreported.pdf"), 9, 6)
    par(mar = c(13, 4, 3, 2))
    barplot(toPlot, las=2, border=FALSE, col=CL_colors_v[names(toPlot)], ylab='n.patients', cex.names = 0.7)
    dev.off()
    
    # Percentage
    stoPlot <- sort(100*toPlot/CMP_ctype_in_COSMIC_Cunt[names(toPlot)])
    stoPlot_all <- sort(100*toPlot_all/CMP_ctype_in_COSMIC_Cunt[names(toPlot_all)], decreasing = TRUE)
    colors_stoPlot_all <- CL_colors_v[match(names(stoPlot_all), names(CL_colors_v))]
    
    # FIGURE 4B
    pdf(file.path(figuresPath, "06_DAMs_in_patients_lacking_act_unreported_perc.pdf"), 7, 5)
    par(mar = c(13, 4, 3, 2))
    barplot(stoPlot, border=FALSE, col=CL_colors_v[names(stoPlot)], las=2, ylab='% of patients', cex.names = 0.7)
    dev.off()
  }

  if (produce_plots && length(toPlot_all) > 0) {

    # SUPPLEMENTARY FIGURE 9I
    pdf(file.path(figuresPath, "06_DAMs_in_patients_lacking_act_perc.pdf"), 7, 5)
    par(mar = c(13, 4, 3, 2))
    barplot(stoPlot_all, border=FALSE, col=colors_stoPlot_all, las=2, log='y', ylim = c(0.01, 50), ylab = "%", cex.names = 0.7, yaxt = "n")
    axis(side = 2, at = c(0.01, 0.05, 0.5, 5, 50), labels = c("0.01", "0.05", "0.50", "5.00", "50.00"), las = 1)
    dev.off()
  }
  
}

# Plotting the percentage of patients carrying a DAM in a known cancer driver gene, where the cancer-type differs
# from those previously associated with that driver (FIGURE 4C)
plot_off_context_driver_DAMs <- function(DAMbearigPatient_in_COSMIC, intogen_drivers, ctype_matching, CMP_ctype_in_COSMIC_Cunt,
    CL_colors_v, figuresPath) {
  
  dat <- as.data.frame(DAMbearigPatient_in_COSMIC, stringsAsFactors = FALSE)
  
  dat$isTheDAM_in_an_unreported_DAMbgs <- as.logical(dat$isTheDAM_in_an_unreported_DAMbgs)
  
  # TRUE means that the DAM occurs in a known IntOGen driver gene
  dat <- dat[dat$isTheDAM_in_an_unreported_DAMbgs,]
  
  # Determine whether the gene is already associated with the cancer type in which the DAM was identified
  known_context <- sapply(seq_len(nrow(dat)), function(i) {
    mapped_intogen_types <- str_trim(str_split(ctype_matching[dat$ct_sel[i], 1], "\\|")[[1]])
    gene_intogen_types <- intogen_drivers$CANCER_TYPE[intogen_drivers$SYMBOL == dat$gene_sel[i]]
    length(intersect(mapped_intogen_types, gene_intogen_types)) > 0
  })
  
  # Keep only off-context DAMs
  dat <- dat[!known_context, ]

  if (nrow(dat) == 0) {
    cat("No off-context DAMs in known driver genes\n")
    return(invisible(NULL))
  }
  
  # Count each patient once per gene and cancer type
  dat <- unique(dat[, c("patients_sel", "ct_sel", "gene_sel")])
  
  counts <- aggregate(patients_sel ~ ct_sel + gene_sel, data = dat, FUN = length)
  colnames(counts) <- c("cancer_type", "gene", "n_patients")
  counts$percentage <- 100 * counts$n_patients / CMP_ctype_in_COSMIC_Cunt[counts$cancer_type]
  counts <- counts[is.finite(counts$percentage) & counts$percentage > 0,]
  
  plot_matrix <- xtabs(percentage ~ cancer_type + gene, data = counts)
  plot_matrix <- plot_matrix[, order(colSums(plot_matrix), decreasing = TRUE), drop = FALSE]
  plot_matrix <- plot_matrix[rowSums(plot_matrix) > 0, , drop = FALSE]
  
  bar_colors <- CL_colors_v[rownames(plot_matrix)]
  bar_colors[is.na(bar_colors)] <- "grey70"
  
  gene_totals <- colSums(plot_matrix)
  max_total <- max(gene_totals)
  
  pdf(file.path(figuresPath, "off_context_DAMs_known_drivers.pdf"), width = 12, height = 7)
  par(mar = c(8, 5, 9, 2), xpd = NA)
  
  bp <- barplot(plot_matrix, col = bar_colors, border = FALSE, space = 0.18, names.arg = rep("", ncol(plot_matrix)),
    ylim = c(-max_total * 0.48, max_total * 1.75), ylab = "% of patients", main = "Off-context DAMs in established\ncancer drivers", font.main = 2)
  
  # Gene names inside bars
  text(x = bp, y = rep(max_total * 0.015, length(bp)), labels = colnames(plot_matrix), srt = 90, adj = c(0, 0.5), font = 2, cex = 0.8)
  
  # Prepare one row per cancer-type segment
  label_data <- do.call(rbind, lapply(seq_len(ncol(plot_matrix)), function(j) {
      present <- which(plot_matrix[, j] > 0)
      if (length(present) == 0) {
        return(NULL)
      }
      cumulative <- cumsum(plot_matrix[, j])
      data.frame(gene_index = j, cancer_index = present, bar_x = bp[j], segment_y = cumulative[present] - plot_matrix[present, j] / 2,
        cancer_type = rownames(plot_matrix)[present], percentage = plot_matrix[present, j], stringsAsFactors = FALSE)}))
  
  # Alternate labels above and below, while spreading nearby labels
  label_data$side <- ifelse(seq_len(nrow(label_data)) %% 2 == 1, "above", "below")
  
  above_index <- which(label_data$side == "above")
  below_index <- which(label_data$side == "below")
  
  # Four vertical lanes above and three below
  label_data$lane <- 1L
  label_data$lane[above_index] <- rep(1:4, length.out = length(above_index))
  label_data$lane[below_index] <- rep(1:3, length.out = length(below_index))
  
  # Horizontal displacement from the bar
  label_data$x_offset <- rep(c(-0.38, 0.38, -0.65, 0.65), length.out = nrow(label_data))
  label_data$label_x <- label_data$bar_x + label_data$x_offset
  label_data$label_y <- ifelse(label_data$side == "above", max_total * (1.10 + 0.16 * label_data$lane), -max_total * (0.12 + 0.11 * label_data$lane))
  
  # Draw labels and angled connector lines
  for (i in seq_len(nrow(label_data))) {
    wrapped_label <- paste(strwrap(label_data$cancer_type[i], width = 20), collapse = "\n")
    elbow_y <- if (label_data$side[i] == "above") {
      label_data$label_y[i] - max_total * 0.04
    } else {
      label_data$label_y[i] + max_total * 0.04
    }
    
    # Vertical line from the bar segment
    segments(
      x0 = label_data$bar_x[i],
      y0 = label_data$segment_y[i],
      x1 = label_data$bar_x[i],
      y1 = elbow_y,
      col = "grey40",
      lwd = 0.7
    )
    
    # Horizontal/diagonal connection to the shifted label
    segments(
      x0 = label_data$bar_x[i],
      y0 = elbow_y,
      x1 = label_data$label_x[i],
      y1 = elbow_y,
      col = "grey40",
      lwd = 0.7
    )
    
    text(x = label_data$label_x[i], y = label_data$label_y[i], labels = wrapped_label, cex = 0.55, adj = if (label_data$x_offset[i] < 0) {
        c(1, 0.5)
      } else {
        c(0, 0.5)
      }
    )
  }
  
  dev.off()
  
  return(list(
    plot_matrix = plot_matrix,
    plot_data = counts,
    label_data = label_data
  ))
}

# Assess DAM tractability / actionability using Open Targets tractability annotations
compute_dam_actionability <- function(raw_dir, tract_file, resultPath, figuresPath, allDAMs_file, totalTestedVariants_file, intOGen_drivers, produce_plots = TRUE){
  
  # Load tractability file
  tract_file <- file.path(raw_dir, tract_file)
  tract <- data.table::fread(tract_file, sep="\t", header=TRUE, quote="", fill=TRUE, showProgress=FALSE)
  
  # Load DAMs and tested variants
  load(file.path(resultPath, allDAMs_file))
  allDAMs <- allDAMs_with_funcImpact_and_clinical_rel
  load(file.path(resultPath, totalTestedVariants_file))
  
  # Get unique Ensembl IDs and corresponding gene symbols
  DAM_unique_ensIds <- unique(totalTestedVariants$ensembl_gene_id[match(allDAMs$GENE, totalTestedVariants$gene_symbol)])
  DAM_unique_genes <- totalTestedVariants$gene_symbol[match(DAM_unique_ensIds, totalTestedVariants$ensembl_gene_id)]
  
  # Extract tractability information for the DAM gene set
  tractability <- tract[match(DAM_unique_ensIds, tract$ensembl_gene_id),
                        c('Top_bucket_sm','Top_bucket_ab',"drug_names_sm","drug_names_ab","drug_names_othercl")]
  
  # Known cancer drivers
  knownDAMbgs <- is.element(DAM_unique_genes, intOGen_drivers)
  
  # Combine into table
  TractabilityTable <- cbind(DAM_unique_ensIds, DAM_unique_genes, knownDAMbgs, tractability)
  
  # Compute min bucket and tiers
  minBucket <- apply(TractabilityTable[,c(4,5)], 1, FUN='min', na.rm=TRUE)
  Tractability_tiers <- cut(minBucket,
                            breaks=c(0,3,6,10),
                            labels=c(1,2,3),
                            right=TRUE,
                            include.lowest=TRUE)
  TractabilityTable <- cbind(TractabilityTable, minBucket, Tractability_tiers)
  
  # Map back to allDAMs
  DAMactionability <- TractabilityTable[match(allDAMs$GENE, TractabilityTable$DAM_unique_genes), ]
  DAMactionability <- data.frame(DAMactionability)
  
  allDAMs <- cbind(allDAMs, DAMactionability[, c(1,4:10)])
  
  # Funnel 2: actionable DAMs summary
  summarize_actionable_dams <- function(allDAMs){
    cat("Actionable DAMs:", length(which(allDAMs$Tractability_tiers==1)), "\n")
    cat("Of which involving unreported cancer drivers:",
        length(which(allDAMs$Tractability_tiers==1 & !allDAMs$in_a_known_driver)), "\n")
    cat("Actionable DAMs with functional impact:",
        length(which(allDAMs$Tractability_tiers==1 & 
                     !is.element(allDAMs$impact_call, c("Unknown impact","Low impact")))), "\n")
    cat("Of which involving unreported cancer drivers:",
        length(which(allDAMs$Tractability_tiers==1 &
                     !is.element(allDAMs$impact_call, c("Unknown impact","Low impact")) &
                     !allDAMs$in_a_known_driver)), "\n")
    cat("Actionable DAMs with functional impact in patients:",
        length(which(allDAMs$Tractability_tiers==1 &
                     !is.element(allDAMs$impact_call, c("Unknown impact","Low impact")) &
                     as.numeric(allDAMs$num_patient_total)>0)), "\n")
    cat("Of which involving unreported cancer drivers:",
        length(which(allDAMs$Tractability_tiers==1 &
                     !is.element(allDAMs$impact_call, c("Unknown impact","Low impact")) &
                     as.numeric(allDAMs$num_patient_total)>0 &
                     !allDAMs$in_a_known_driver)), "\n")
    cat("Actionable DAMs with functional impact in matching ct patients:",
        length(which(allDAMs$Tractability_tiers==1 &
                     !is.element(allDAMs$impact_call, c("Unknown impact","Low impact")) &
                     as.numeric(allDAMs$n_tumour_type_matching_mutant_patients)>0)), "\n")
    cat("Of which involving unreported cancer drivers:",
        length(which(allDAMs$Tractability_tiers==1 &
                     !is.element(allDAMs$impact_call, c("Unknown impact","Low impact")) &
                     as.numeric(allDAMs$n_tumour_type_matching_mutant_patients)>0 &
                     !allDAMs$in_a_known_driver)), "\n")
  }
  
  summarize_actionable_dams(allDAMs)
  
  # Barplot of tractability tiers
  plot_tractability_barplot <- function(TractabilityTable, figuresPath, produce_plots = TRUE){
    
    BAR1 <- table(factor(TractabilityTable$knownDAMbgs[TractabilityTable$Tractability_tiers == 1], levels = c(FALSE, TRUE)))
    BAR2 <- table(factor(TractabilityTable$knownDAMbgs[TractabilityTable$Tractability_tiers == 2], levels = c(FALSE, TRUE)))
    BAR3 <- table(factor(TractabilityTable$knownDAMbgs[TractabilityTable$Tractability_tiers == 3], levels = c(FALSE, TRUE)))
    toPlot <- rbind(BAR3,BAR2,BAR1)
    
    colnames(toPlot) <- c('unreported','known')
    rownames(toPlot) <- c('>6','4-6','1-3')
    
    # FIGURE 4D
    if (produce_plots && sum(toPlot) > 0) {
      pdf(file.path(figuresPath, "DAMbgs_druggability.pdf"), 4, 8)
      barplot(t(toPlot), beside=TRUE, log='y',
              ylab=paste('n. of DAMbearing genes (out of', sum(toPlot), 'with druggability information)'))
      dev.off()
    }

  }
  
  plot_tractability_barplot(TractabilityTable, figuresPath, produce_plots)
  
  # Restrict to DAMs with at least possible functional impact
  ii <- which(!allDAMs$impact_call %in% c("Low impact","Unknown impact"))
  allDAMs_res <- allDAMs[ii, ]
  DAMactionability <- DAMactionability[ii, ]

  DAMactionabilityKnown <- DAMactionability[which(allDAMs_res$impact_call!='Low impact' &
                                           allDAMs_res$impact_call!='Unknown impact' &
                                             is.element(allDAMs_res$GENE,intOGen_drivers)), ]


  DAMactionabilityUnknown <- DAMactionability[which(allDAMs_res$impact_call!='Low impact' &
                                            allDAMs_res$impact_call!='Unknown impact' &
                                             !is.element(allDAMs_res$GENE,intOGen_drivers)), ]
  
  list(allDAMs=allDAMs, DAMactionability=DAMactionability)
}
