# -------- SUMMARY OF PATIENT DATA -------- # 

library(stringr)
library(readr)
library(openxlsx)

# Split the aggregated patient-count matrix into variants occurring in known driver genes
# versus variants occurring in genes not currently listed as known drivers
split_aggregated_by_driver <- function(aggregated_counts_aggr_comp, intogen_drivers_s) {
  gn <- unlist(lapply(str_split(row.names(aggregated_counts_aggr_comp), " "), function(x) { x[1] }))

  count_agg_known <- aggregated_counts_aggr_comp[is.element(gn, intogen_drivers_s), , drop = FALSE]
  count_agg_unreported <- aggregated_counts_aggr_comp[!is.element(gn, intogen_drivers_s), , drop = FALSE]

  return(list(
    gn = gn,
    count_agg_known = count_agg_known,
    count_agg_unreported = count_agg_unreported
  ))
}

# Build a color vector for the cancer types represented in the known-driver matrix.
# Colors are assigned from the clinical color table 'clc'
compute_cols_for_known <- function(count_agg_known, clc) {
  CL_colors_v <- clc$Color_hex_code
  names(CL_colors_v) <- clc$Cancer.Type
  COLS <- CL_colors_v[colnames(count_agg_known)]
  COLS[is.na(COLS)] <- "gray"
  return(COLS)
}

# Combine functional impact predictions with patient recurrence and cancer-type matching
# information from intOGen, COSMIC, and the aggregated comparison matrix.
# The function annotates each DAM with:
# - number of mutant patients in each source
# - total and aggregated recurrence
# - cancer types in which the mutation is observed
# - whether patient tumour types match the analysed DAM cancer type
# - whether the DAM belongs to a known driver gene 
build_funcImpact_and_clinical_rel <- function(allDAMs_with_SIFT_Polyphen, intogen_counts_aggr, cosmic_counts_aggr,
                                              aggregated_counts_aggr_comp, ctype_matching, intogen_drivers_s) {
  
  # Initialize numeric columns for patient counts
  allDAMs_with_SIFT_Polyphen$num_patient_COSMIC <- rep(0, nrow(allDAMs_with_SIFT_Polyphen))
  allDAMs_with_SIFT_Polyphen$num_patient_intOGen <- rep(0, nrow(allDAMs_with_SIFT_Polyphen))
  allDAMs_with_SIFT_Polyphen$num_patient_total <- rep(0, nrow(allDAMs_with_SIFT_Polyphen))
  allDAMs_with_SIFT_Polyphen$num_patient_aggregated <- rep(0, nrow(allDAMs_with_SIFT_Polyphen))

  # Initialize text columns describing the tumour types where the mutation is observed
  allDAMs_with_SIFT_Polyphen$mut_type_COSMIC <- rep("", nrow(allDAMs_with_SIFT_Polyphen))
  allDAMs_with_SIFT_Polyphen$mut_type_Intogen <- rep("", nrow(allDAMs_with_SIFT_Polyphen))
  
  # Build variant signatures
  allVar <- paste(allDAMs_with_SIFT_Polyphen$GENE,
                  unlist(lapply(str_split(allDAMs_with_SIFT_Polyphen$var, "p."), function(x) { x[2] })))
  
  # Recurrence vectors
  NumIntoGen <- rep(0, length(allVar))
  NumCosmic <- rep(0, length(allVar))
  NumAggregated <- rep(0, length(allVar))

  # Retrieve total counts per variant from each matrix
  NumIntoGen <- rowSums(intogen_counts_aggr)[allVar]
  NumCosmic <- rowSums(cosmic_counts_aggr)[allVar]
  NumAggregated <- rowSums(aggregated_counts_aggr_comp)[allVar]

  # Replace missing values for non-observed variants with zero
  NumIntoGen[is.na(NumIntoGen)] <- 0
  NumCosmic[is.na(NumCosmic)] <- 0
  NumAggregated[is.na(NumAggregated)] <- 0
  
  # Vectors storing tumour-type labels
  MutTypesIntoGen <- rep("", length(allVar))
  MutTypesCosmic <- rep("", length(allVar))
  MutTypesAggregated <- rep("", length(allVar))
  MutatedInMatchingPatients <- rep(FALSE, length(allVar)) 

  # For each intOGen variant, record the cancer types with non-zero counts
  tmp <- unlist(lapply(1:nrow(intogen_counts_aggr), function(x) {
    x <- intogen_counts_aggr[x, ]
    return(paste(colnames(x[which(x > 0)]), collapse = " | "))
  }))
  names(tmp) <- rownames(intogen_counts_aggr)
  MutTypesIntoGen <- tmp[allVar]
  MutTypesIntoGen[is.na(MutTypesIntoGen)] <- ""

  # For each COSMIC variant, record the cancer types with non-zero counts
  tmp <- unlist(lapply(1:nrow(cosmic_counts_aggr), function(x) {
    x <- cosmic_counts_aggr[x, ]
    return(paste(names(which(x > 0)), collapse = " | "))
  }))
  names(tmp) <- rownames(cosmic_counts_aggr)
  MutTypesCosmic <- tmp[allVar]
  MutTypesCosmic[is.na(MutTypesCosmic)] <- ""

  # For each aggregated variant, record all combined cancer-type groups with non-zero counts
  tmp <- unlist(lapply(1:nrow(aggregated_counts_aggr_comp), function(x) {
    x <- aggregated_counts_aggr_comp[x, ]
    return(paste(names(which(x > 0)), collapse = " | "))
  }))
  names(tmp) <- rownames(aggregated_counts_aggr_comp)
  MutTypesAggregated <- tmp[allVar]
  MutTypesAggregated[is.na(MutTypesAggregated)] <- ""

  # Flag whether the analysed cancer type overlaps with any patient tumour type
  matching_tumour_patients <-
    unlist(lapply(1:length(allVar), function(x) {
      length(intersect(allDAMs_with_SIFT_Polyphen$ctype[x], str_trim(str_split(MutTypesAggregated[x], "[|]")[[1]]))) > 0
    }))

  # Return the source-specific tumour categories corresponding to the matched cancer type
  matching_tumour_types_in_patients <-
    unlist(lapply(1:length(allVar), function(x) {

      tmp <- intersect(allDAMs_with_SIFT_Polyphen$ctype[x],
                       str_trim(str_split(MutTypesAggregated[x], "[|]")[[1]]))

      if (length(tmp) == 0) {
        return("")
      } else {

        intOGenTypes <- str_trim(unlist(str_split(ctype_matching[tmp, 1], "[|]")))
        cosmicTypes <- str_trim(unlist(str_split(ctype_matching[tmp, 2], "[|]")))

        matchingTypes <- paste(c(intOGenTypes, cosmicTypes), collapse = " | ")

        return(matchingTypes)
      }

    }))

  # Count how many mutant patients are observed specifically in cancer types matching the DAM context
  n_tumour_type_matching_mutant_patients <-
    unlist(lapply(1:length(allVar), function(x) {
      tmp <- intersect(allDAMs_with_SIFT_Polyphen$ctype[x],
                       str_trim(str_split(MutTypesAggregated[x], "[|]")[[1]]))

      if (length(tmp) == 0) {
        return("")
      } else {

        intOGenTypes <- str_trim(unlist(str_split(ctype_matching[tmp, 1], "[|]")))
        cosmicTypes <- str_trim(unlist(str_split(ctype_matching[tmp, 2], "[|]")))

        if (is.element(allVar[x], rownames(cosmic_counts_aggr))) {
          cosmicCounts <- sum(cosmic_counts_aggr[allVar[x],
                                               intersect(cosmicTypes, colnames(cosmic_counts_aggr))],
                              na.rm = TRUE)
        } else { cosmicCounts <- 0 }
        if (is.element(allVar[x], rownames(intogen_counts_aggr))) {
          intogenCounts <- sum(intogen_counts_aggr[allVar[x],
                                                 intersect(intOGenTypes, colnames(intogen_counts_aggr))],
                               na.rm = TRUE)
        } else { intogenCounts <- 0 }

        res <- cosmicCounts + intogenCounts

        return(res)
      }
    }))

  # Build final annotated table by extending the functional-impact table
  allDAMs_with_funcImpact_and_clinical_rel <- allDAMs_with_SIFT_Polyphen
  allDAMs_with_funcImpact_and_clinical_rel$num_patient_COSMIC <- NumCosmic
  allDAMs_with_funcImpact_and_clinical_rel$num_patient_intOGen <- NumIntoGen
  allDAMs_with_funcImpact_and_clinical_rel$num_patient_total <- NumIntoGen + NumCosmic
  allDAMs_with_funcImpact_and_clinical_rel$num_patient_aggregated <- NumAggregated
  allDAMs_with_funcImpact_and_clinical_rel$mut_type_Intogen <- MutTypesIntoGen
  allDAMs_with_funcImpact_and_clinical_rel$mut_type_COSMIC <- MutTypesCosmic
  allDAMs_with_funcImpact_and_clinical_rel$mut_type_aggregated <- MutTypesAggregated
  allDAMs_with_funcImpact_and_clinical_rel$matching_tumour_patients <- matching_tumour_patients
  allDAMs_with_funcImpact_and_clinical_rel$matching_tumour_types_in_patients <- matching_tumour_types_in_patients
  allDAMs_with_funcImpact_and_clinical_rel$n_tumour_type_matching_mutant_patients <- n_tumour_type_matching_mutant_patients

  # Flag DAMs occurring in genes already classified as known drivers
  allDAMs_with_funcImpact_and_clinical_rel$in_a_known_driver <-
    is.element(allDAMs_with_funcImpact_and_clinical_rel$GENE, intogen_drivers_s)

  return(list(
    allDAMs_with_funcImpact_and_clinical_rel = allDAMs_with_funcImpact_and_clinical_rel,
    allVar = allVar,
    MutTypesAggregated = MutTypesAggregated
  ))
}

# Generate summary pie charts showing how many DAMs are observed in cancer patients,
# separated by source. The function produces one chart for all DAMs and one for 
# DAMs outside known driver genes.
plot_and_report_patient_observation_pies <- function(allDAMs_with_funcImpact_and_clinical_rel, figuresPath, intogen_drivers_s, produce_plots = TRUE) {

  # Build variant identifiers
  allVar <- paste(allDAMs_with_funcImpact_and_clinical_rel$GENE,
                  unlist(lapply(str_split(allDAMs_with_funcImpact_and_clinical_rel$var, "p."), function(x) { x[2] })))

  nDAMs <- length(allVar)

  # Count DAMs seen in both datasets
  nDAMs_in_both_intogen_and_cosmic <-
    length(which(allDAMs_with_funcImpact_and_clinical_rel$num_patient_COSMIC > 0 &
                   allDAMs_with_funcImpact_and_clinical_rel$num_patient_intOGen > 0))

  # Count DAMs seen only in intOGen
  nDAMs_in_intogen_only <-
    length(which(allDAMs_with_funcImpact_and_clinical_rel$num_patient_COSMIC == 0 &
                   allDAMs_with_funcImpact_and_clinical_rel$num_patient_intOGen > 0))

  # Count DAMs seen only in COSMIC
  nDAMs_in_cosmic_only <-
    length(which(allDAMs_with_funcImpact_and_clinical_rel$num_patient_COSMIC > 0 &
                   allDAMs_with_funcImpact_and_clinical_rel$num_patient_intOGen == 0))
  
  # Count DAMs not observed in either patient dataset
  nDAMs_not_in_patients <-
    length(which(allDAMs_with_funcImpact_and_clinical_rel$num_patient_COSMIC == 0 &
                   allDAMs_with_funcImpact_and_clinical_rel$num_patient_intOGen == 0))
  
  # Plot pie chart for all DAMs
  # SUPPLEMENTARY FIGURE 20A
  if (produce_plots) {
    pdf(file.path(figuresPath, "04_in_patients_clinical_rap_of_all_DAMs.pdf"), 7, 5)
    pie(c(nDAMs_in_both_intogen_and_cosmic, nDAMs_in_cosmic_only, nDAMs_in_intogen_only, nDAMs_not_in_patients),
        main = paste(nDAMs, "cancer-type-specific DAMs"),
        labels = c("in both", "in COSMIC", "in intOGen", "not observed"),
        col = c("red", "orange", "pink", "gray"), border = FALSE)
    dev.off()
  }

  cat(paste("of the ", nDAMs, " unique DAMs, ",
            nDAMs - nDAMs_not_in_patients, " (",
            round(100 * (nDAMs - nDAMs_not_in_patients) / nDAMs, 2),
            "%) are observed in cancer patients (COSMIC or intOGen)", sep = ""), "\n")

  cat(paste(nDAMs_in_both_intogen_and_cosmic, " (",
            round(100 * nDAMs_in_both_intogen_and_cosmic / nDAMs, 2),
            "%) are observed in both", sep = ""), "\n")

  # Same summary ummary for novel / unreported DAMs only
  novelOnly <- allDAMs_with_funcImpact_and_clinical_rel[!allDAMs_with_funcImpact_and_clinical_rel$in_a_known_driver, ]

  allVar <- paste(novelOnly$GENE,
                  lapply(str_split(novelOnly$var, "p."), function(x) { x[2] }),
                  sep = " ")

  nDAMs <- length(allVar)

  nDAMs_in_both_intogen_and_cosmic <- length(which(novelOnly$num_patient_COSMIC > 0 &
                                                    novelOnly$num_patient_intOGen > 0))

  nDAMs_in_intogen_only <- length(which(novelOnly$num_patient_COSMIC == 0 &
                                          novelOnly$num_patient_intOGen > 0))

  nDAMs_in_cosmic_only <- length(which(novelOnly$num_patient_COSMIC > 0 &
                                         novelOnly$num_patient_intOGen == 0))

  nDAMs_not_in_patients <- length(which(novelOnly$num_patient_COSMIC == 0 &
                                          novelOnly$num_patient_intOGen == 0))
  
  # Plotting
  # SUPPLEMENTARY FIGURE 20B
  if (produce_plots && nDAMs > 0) {
    pdf(file.path(figuresPath, "04_in_patients_clinical_rap_of_novel_DAMs.pdf"), 7, 5)
    pie(c(nDAMs_in_both_intogen_and_cosmic, nDAMs_in_cosmic_only, nDAMs_in_intogen_only, nDAMs_not_in_patients),
        main = paste(nDAMs, "cancer-type-specific DAMs"),
        labels = c("in both", "in COSMIC", "in intOGen", "not observed"),
        col = c("red", "orange", "pink", "gray"), border = FALSE)
    dev.off()
  }

  cat(paste("of the ", nDAMs, " unique DAMs involving unreported driver genes, ",
            nDAMs - nDAMs_not_in_patients, " (",
            round(100 * (nDAMs - nDAMs_not_in_patients) / nDAMs, 2),
            "%) are observed in cancer patients (COSMIC or intOGen)", sep = ""), "\n")

  cat(paste(nDAMs_in_both_intogen_and_cosmic, " (",
            round(100 * nDAMs_in_both_intogen_and_cosmic / nDAMs, 2),
            "%) are observed in both", sep = ""), "\n")

  invisible(list(
    nDAMs_all = nDAMs,
    novelOnly = novelOnly
  ))
}

# Print the main counts used in the Funnel / Figure 1D summary
cat_funnel_figure1D_stats <- function(allDAMs_with_funcImpact_and_clinical_rel, intogen_drivers_s) {

  cat(paste("n cancer type specific DAMs", nrow(allDAMs_with_funcImpact_and_clinical_rel)), "\n")
  cat(paste("of which in unreportd DAMbgs",
            length(which(!is.element(allDAMs_with_funcImpact_and_clinical_rel$GENE, intogen_drivers_s)))), "\n")

  cat(paste("with predicted functional impact",
            length(which(!is.element(allDAMs_with_funcImpact_and_clinical_rel$impact_call,
                                     c("Unknown impact", "Low impact"))))), "\n")

  cat(paste("of which in unreported DAMbgs",
            length(which(!is.element(allDAMs_with_funcImpact_and_clinical_rel$impact_call,
                                     c("Unknown impact", "Low impact")) &
                           !is.element(allDAMs_with_funcImpact_and_clinical_rel$GENE, intogen_drivers_s)))), "\n")

  cat(paste("with predicted functional impact and mutated patients",
            length(which(!is.element(allDAMs_with_funcImpact_and_clinical_rel$impact_call,
                                     c("Unknown impact", "Low impact")) &
                           allDAMs_with_funcImpact_and_clinical_rel$num_patient_total > 0))), "\n")

  cat(paste("of which in unreported DAMbgs",
            length(which(!is.element(allDAMs_with_funcImpact_and_clinical_rel$impact_call,
                                     c("Unknown impact", "Low impact")) &
                           allDAMs_with_funcImpact_and_clinical_rel$num_patient_total > 0 &
                           !is.element(allDAMs_with_funcImpact_and_clinical_rel$GENE, intogen_drivers_s)))), "\n")

  cat(paste("with predicted functional impact and mutated patients (matching cancer type)",
            length(which(!is.element(allDAMs_with_funcImpact_and_clinical_rel$impact_call,
                                     c("Unknown impact", "Low impact")) &
                           allDAMs_with_funcImpact_and_clinical_rel$matching_tumour_patients))), "\n")

  cat(paste("of which in unreported DAMbgs",
            length(which(!is.element(allDAMs_with_funcImpact_and_clinical_rel$impact_call,
                                     c("Unknown impact", "Low impact")) &
                           allDAMs_with_funcImpact_and_clinical_rel$matching_tumour_patients &
                           !is.element(allDAMs_with_funcImpact_and_clinical_rel$GENE, intogen_drivers_s)))), "\n")
}

# Compute per-variant prevalence across patient datasets and generate summary plots
compute_uvar_prevalence_and_outputs <- function(allDAMs_with_SIFT_Polyphen, allDAMs_with_funcImpact_and_clinical_rel, intogen_counts_aggr,
                                               cosmic_counts_aggr, path_results, figuresPath, intogen_drivers_s, produce_plots = TRUE) {
  
  # Build variant signatures 
  allVar <- paste(allDAMs_with_SIFT_Polyphen$GENE,
                  unlist(lapply(str_split(allDAMs_with_SIFT_Polyphen$var, "p."), function(x) { x[2] })))

  uvar <- unique(allVar)

  # For each unique variant, compute total prevalence and prevalence in matching tumour types
  uvarPrevalence <- do.call(rbind, lapply(uvar, function(x) {
    ii <- which(allVar == x)
    
    # Total number of mutant patients across all sources
    TotalMutantPatients <- allDAMs_with_funcImpact_and_clinical_rel$num_patient_total[ii[1]]
    
    # Tumour types matching the analysed cancer-type context
    types <- str_trim(unlist(str_split(allDAMs_with_funcImpact_and_clinical_rel$matching_tumour_types_in_patients[ii], "[|]")))
    
    # Sum COSMIC counts over the matching tumour types
    if (is.element(x, rownames(cosmic_counts_aggr))) {
      cosmicCounts <- sum(cosmic_counts_aggr[x, intersect(types, colnames(cosmic_counts_aggr))], na.rm = TRUE)
    } else { cosmicCounts <- 0 }
    # Sum intOGen counts over the matching tumour types
    if (is.element(x, rownames(intogen_counts_aggr))) {
      intogenCounts <- sum(intogen_counts_aggr[x, intersect(types, colnames(intogen_counts_aggr))], na.rm = TRUE)
    } else { intogenCounts <- 0 }

    TotalMutantPatientsMatchingCtype <- cosmicCounts + intogenCounts

    return(c(TotalMutantPatients,
             TotalMutantPatientsMatchingCtype,
             is.element(unlist(str_split(x, " "))[1], intogen_drivers_s)))
  }))

  colnames(uvarPrevalence) <- c("n mutant patients", "with matching ctype", "is a known driver")
  rownames(uvarPrevalence) <- uvar
  
  # Rank variants by total patient prevalence
  uvarPrevalence <- uvarPrevalence[order(uvarPrevalence[, 1], decreasing = TRUE), ]
  
  # Plot the most prevalent known-driver DAMs
  # SUPPLEMENTARY FIGURE 21A
  howmany <- head(which(uvarPrevalence[, 3] > 0), 38)
  if (produce_plots && length(howmany) > 0) {
    pdf(file.path(figuresPath, "04_Known_DAMs_most_frequently_obser_in_patients.pdf"), 10, 5)
    par(mar = c(8, 4, 2, 2))
    barplot(t(cbind(uvarPrevalence[howmany, 2],
                    uvarPrevalence[howmany, 1] - uvarPrevalence[howmany, 2])) + 1,
            las = 2, ylab = "n. mutant patients + 1", col = c("blue", "gray"),
            border = FALSE, log = "y")
    legend("topright", legend = c("matching cancer types", "others"),
          fill = c("blue", "gray"), border = FALSE)
    dev.off()
  }
  
  # Plot the most prevalent unreported DAMs 
  # SUPPLEMENTARY FIGURE 21B
  howmany <- head(which(uvarPrevalence[, 3] == 0 & uvarPrevalence[, 2] > 0), 50)
  if (produce_plots && length(howmany) > 0) {
    pdf(file.path(figuresPath, "04_unreported_DAMs_most_frequently_obser_in_patients.pdf"), 15, 8)
    par(mar = c(10, 4, 2, 2))
    barplot(t(cbind(uvarPrevalence[howmany, 2],
                    uvarPrevalence[howmany, 1] - uvarPrevalence[howmany, 2])),
            las = 2, ylab = "n. mutant patients + 1", col = c("blue", "gray"),
            border = FALSE)
    legend("topright", legend = c("matching cancer types", "others"),
          fill = c("blue", "gray"), border = FALSE)
    dev.off()
  }

  # Prepare data frame used later for the circular summary plot 
  nCtypSpecAnalysis <- summary(as.factor(allVar), length(unique(allVar)))[uvar]

  df <- data.frame(
    ID = uvar,
    impactCall = allDAMs_with_funcImpact_and_clinical_rel$impact_call[match(uvar, allVar)],
    PrevalenceInTumoursOther = uvarPrevalence[uvar, 1] - uvarPrevalence[uvar, 2],
    PrevalenceInTumoursMatching = uvarPrevalence[uvar, 2],
    nCtypSpecAnalysis[uvar],
    isAknownDriver = (uvarPrevalence[uvar, 3] > 0),
    func_impact = allDAMs_with_funcImpact_and_clinical_rel$impact_call[match(uvar, allVar)]
  )
  
  # Compute total patient prevalence and keep variants with enough matching evidence
  df$totalPatients <- (df$PrevalenceInTumoursOther + df$PrevalenceInTumoursMatching)
  df <- df[which(df$PrevalenceInTumoursMatching > 1), ]
  df <- df[order(df$totalPatients, decreasing = TRUE), ]

  save(df, file = file.path(path_results, "_DF_for_circular_plot.RData"))

  return(list(
    uvarPrevalence = uvarPrevalence,
    df_circular = df
  ))
}

# Helper function for circular plot
polar_rect <- function(theta, w, r0, r1, col, border = NA, nside = 40){
  ang1 <- seq(theta - w/2, theta + w/2, length.out = nside)
  ang2 <- rev(ang1)
  x <- c(r0*cos(ang2), r1*cos(ang1))
  y <- c(r0*sin(ang2), r1*sin(ang1))
  polygon(x, y, col = col, border = border)
}

# Create the main circular summary plot (FIGURE 4A) showing, for each DAM:
# - total number of mutant patients
# - fraction of patients in matching versus non-matching tumour types
# - number of cancer-type-specific analyses uncovering the DAM
# - functional impact class
# - whether the associated gene is a known or unreported driver
# The outward bars encode prevalence in patients, while inward bars encode how many
# cancer-type-specific analyses contributed to that DAM.
plot_circular_dams <- function(df, figuresPath) {

  if (is.null(df) || nrow(df) == 0) {
    cat("No variants available for the circular DAM plot\n")
    return(invisible(NULL))
  }

  # Total number of mutant patients
  Tot <- df$totalPatients
  # Log-scaled prevalence in matching tumour types
  Matching <- log10(df$PrevalenceInTumoursMatching)
  # Log-scaled prevalence in other tumour types
  Other <- log10(df$PrevalenceInTumoursOther)

  # Driver-gene status and functional impact labels
  isAknownDriver <- df$isAknownDriver
  funcImp <- df$impactCall
  IDs <- df$ID
  
  # Number of cancer-type-specific analyses supporting the DAM
  Inward <- as.numeric(df$nCtypSpecAnalysis)
  Inward[!is.finite(Inward) | Inward < 0] <- 0

  functionalSymbol <- rep(NA, length(funcImp))
  functionalSymbol[which(funcImp == "Low impact")] <- 16
  functionalSymbol[which(funcImp == "Possible impact")] <- 17
  functionalSymbol[which(funcImp == "Moderate impact")] <- 15
  functionalSymbol[which(funcImp == "High impact")] <- 8

  Matching[!is.finite(Matching) | Matching < 0] <- 0
  Other[!is.finite(Other) | Other < 0] <- 0

  Total <- Matching + Other

  # Order variants by total prevalence
  ord <- order(Tot, decreasing = TRUE)
  Matching <- Matching[ord]; Other <- Other[ord]; IDs <- IDs[ord]; Total <- Total[ord]
  isAknownDriver <- isAknownDriver[ord]
  Tot <- Tot[ord]
  # not in the original script !!
  # functionalSymbol <- functionalSymbol[ord]
  # Inward <- Inward[ord]

  n <- length(IDs)

  # Geometry / scaling 
  r0 <- 0.35
  r_outer <- 0.60

  scaleFac <- if (max(Total) > 0) (r_outer - r0) / max(Total) else 1

  # Angles and width
  theta <- -seq(0, 2*pi, length.out = n + 1)[-(n + 1)]
  theta <- theta - pi/2   

  barFrac <- 0.90
  bar_w <- (2*pi/n) * barFrac

  # Rotation 
  theta <- theta - pi/2

  pdf(file.path(figuresPath, "circular_plot.pdf"), 14, 14)

  pad <- 1.12
  op <- par(mar = c(1,1,1,1))
  plot(0, 0, type = "n", asp = 1,
       xlim = c(-r_outer*pad, r_outer*pad),
       ylim = c(-r_outer*pad, r_outer*pad),
       axes = FALSE, xlab = "", ylab = "")

  # colors 
  col_total <- adjustcolor("#00a44d", alpha.f = 0.40)
  col_inner <- "#00a44d"

  # total on log scale 
  logTot <- log10(Tot)

  scaleFac <- if (max(logTot, na.rm = TRUE) > 0) (r_outer - r0) / max(logTot, na.rm = TRUE) else 1
  
  # Draw outward prevalence bars
  for (i in seq_len(n)) {
    if (!is.finite(logTot[i]) || logTot[i] <= 0) next

    # full bar = log10(Matching + Other)
    r_tot <- r0 + scaleFac * logTot[i]
    polar_rect(theta[i], bar_w, r0, r_tot, col = col_total, border = NA)

    # inner overlay = log10(Matching) (composition cue)
    if (is.finite(Matching[i]) && Matching[i] > 0) {
      r_split <- r0 + scaleFac * Matching[i]
      polar_rect(theta[i], bar_w, r0, r_split, col = col_inner, border = NA)
    }
  }

  # Add labels, total counts, and functional-impact symbols
  COLS <- rep(adjustcolor("#264292", alpha = 0.70), n)
  COLS[which(!isAknownDriver)] <- "#f79421"

  oldTot <- Tot
  oldTot[duplicated(Tot)] <- NA
  
  for (i in seq_len(n)) {
    r_tip <- r0 + scaleFac * log10(Tot[i])
    r_lab <- r_tip * 1.14
    r_lab1 <- r_tip * 1.015
    r_lab2 <- r_tip * 1.10

    ang <- theta[i]  # rad
    x <- r_lab * cos(ang)
    y <- r_lab * sin(ang)

    x1 <- r_lab1 * cos(ang)
    y1 <- r_lab1 * sin(ang)

    x2 <- r_lab2 * cos(ang)
    y2 <- r_lab2 * sin(ang)

    flip_right <- (cos(ang) < 0)

    srt_deg <- ang * 180 / pi
    if (flip_right) {
      srt_deg <- srt_deg + 180
      adj_val <- 1
    } else {
      adj_val <- 0
    }
    
    # Variant ID label
    text(x, y, labels = IDs[i], srt = srt_deg, adj = adj_val, cex = 1, col = COLS[i], font = 2)
    # Total patient count label
    text(x1, y1, labels = oldTot[i], srt = srt_deg, adj = adj_val, cex = 0.9, col = "black", font = 2)
    # Functional impact symbol
    points(x2, y2, pch = functionalSymbol[i], cex = 1.3, col = COLS[i])
  }

  # Inward bars 
  max_inward_depth <- 0.5 * r0
  scaleFac_in <- if (max(Inward) > 0) max_inward_depth / max(Inward) else 0

  col_inward <- adjustcolor("#A086B7", alpha.f = 0.7)

  for (i in seq_len(n)) {
    t_in <- scaleFac_in * Inward[i]
    if (t_in > 0) {
      r_in1 <- r0
      r_in0 <- max(r0 - t_in, 0)
      polar_rect(theta[i], bar_w, r_in0, r_in1, col = col_inward, border = NA)
    }
  }
  
  # choose which bar defines the axis
  ang0 <- theta[1] + (theta[length(theta)] + theta[1]) / 2
  
  # length of the axis: from baseline r0 inward to the centre
  r_start <- r0
  r_end <- 0.17
  
  # compute coordinates
  x_line <- c(r_start, r_end) * cos(ang0)
  y_line <- c(r_start, r_end) * sin(ang0)

  lines(x_line, y_line, lwd = 1.5, col = "black", lty = 1)

  # Choose tick values in data units 
  ticks_in <- pretty(c(0, max(Inward, na.rm = TRUE)), n = 5)
  ticks_in[1] <- NA

  # Convert to radii along the inward axis
  r_ticks <- r0 - scaleFac_in * ticks_in
  # Clip to [0, r0] so we don't go beyond centre
  r_ticks <- pmax(r_ticks, 0)
  
  # Tick mark length (in radial units)
  tick_len <- 0.02 * r0
  
  # Draw tick marks (small segments perpendicular to the axis)
  for (rt in r_ticks) {
    # along axis
    x0 <- rt * cos(ang0); y0 <- rt * sin(ang0)
    # short perpendicular offset (rotate by +/-90°)
    dx <- tick_len * cos(ang0 + pi/2)
    dy <- tick_len * sin(ang0 + pi/2)
    segments(x0 - dx, y0 - dy, x0 + dx, y0 + dy, col = "black", lwd = 1)
  }
  
  # Labels (data units) slightly inward from the tick
  lab_offset <- 0.04 * r0
  for (k in seq_along(ticks_in)) {
    rt <- r_ticks[k]
    lab <- ticks_in[k]
    xL <- (rt) * cos(ang0)
    yL <- (rt - lab_offset) * sin(ang0)
    text(xL, yL - 0.015, labels = lab, cex = 0.8, adj = c(0.5, 0.5))
  }

  legend(x = -0.13, y = 0.2,
         legend = c(paste("n. of Cancer-type-specific analyses\nunveling the DAM"),
                    paste("n. of patients bearing the DAM\n[matching cancer type(s)]"),
                    paste("n. of patients bearing the DAM\n(other cancer types)")),
         fill = c(adjustcolor("#A086B7", alpha.f = 0.7), col_inner, col_total),
         bty = "n", border = FALSE, y.intersp = 1.8)

  legend(x = -0.05, y = -0.01,
         legend = c("Low", "Possible", "Moderate", "High"),
         border = FALSE, pch = c(16, 17, 15, 8), title = "functional impact")

  legend(x = -0.15, y = -0.18,
         legend = c("Known cancer driver", "Unreported cancer driver"),
         text.col = c(adjustcolor("#264292", alpha = 0.70), "#f79421"),
         pch = NA,
         bty = "n", y.intersp = 1.8, text.font = 2)

  dev.off()

  invisible(TRUE)
}
