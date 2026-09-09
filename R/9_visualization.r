# -------- GLOBAL SUMMARIES AND VISUALIZATION -------- #

library(ggplot2)
library(ggrepel)
library(openxlsx)

# How many cancer types each DAM appears in
# It counts DAM frequency across ctype and generates histogram/barplot
plot_nDAMs_in_nCtypes <- function(allDAMs, figuresPath, produce_plots = TRUE) {

  occurrencies <- summary(as.factor(
    sort(summary(as.factor(paste(allDAMs$GENE,allDAMs$var)),10000))
  ))

  # FIGURE 2E
  if (produce_plots) {
    pdf(file.path(figuresPath,'nDAMs_in_nCtypes.pdf'),6,6)
    barplot(
      rev(-log10(occurrencies+1)),
      border = FALSE
    )
    invisible(dev.off())
  }

  cat(occurrencies['1'], "out of", sum(occurrencies), "unique DAMs (", 100 * occurrencies['1'] / sum(occurrencies),
  "%) are detected in a single ctype specific analysis\n")

  return(occurrencies)
}

# Computes recurrence frequency of DAMbgs across all cancer-type analyses
compute_DAMbgs_across_analyses <- function(allHits) {

  DAMbgsAcrossNanalysis <- sort(summary(as.factor(allHits$GENE), length(unique(allHits$GENE))), decreasing = TRUE)

  occurrencies <- summary(as.factor(DAMbgsAcrossNanalysis))

  cat(occurrencies['1'], 'out of', sum(occurrencies),
      'DAMbgs (', 100*occurrencies['1']/sum(occurrencies),
      '%) are detected in a single ctype specific analysis\n')

  cat(sum(occurrencies[2:length(occurrencies)]), ' DAMbgs (',
      100*sum(occurrencies[2:length(occurrencies)])/sum(occurrencies),
      '%) are detected in multiple ctypes\n')

  return(list(
    DAMbgsAcrossNanalysis = DAMbgsAcrossNanalysis,
    occurrencies = occurrencies
  ))
}

# Plots distribution of DAMbgs recurrence across cancer types
plot_nDAMbgs_in_nCtypes <- function(DAMbgsAcrossNanalysis, figuresPath, produce_plots = TRUE) {
  
  # FIGURE 2E
  if (produce_plots) {
    (file.path(figuresPath,'nDAMbgs_in_nCtypes.'),6,6)
    barplot(log10(summary(as.factor(DAMbgsAcrossNanalysis))+1),
            border=FALSE,col='darkgray')
    invisible(dev.off())
  }
}

# Summarizes frequently recurring DAMbgs and
# compares them to reference IntOGen driver gene list.
# Prints overlap statistics
summarize_DAMbgs_multiplicity <- function(DAMbgsAcrossNanalysis, intogen_drivers) {

  multipleDetectedDAMbgs <- names(which(DAMbgsAcrossNanalysis>1))

  multipleKNOWN <- intersect(multipleDetectedDAMbgs, intogen_drivers)
  multipleUnreported <- setdiff(multipleDetectedDAMbgs, multipleKNOWN)

  cat('Of these ', length(multipleUnreported), ' (',
      100*length(multipleUnreported)/length(multipleDetectedDAMbgs),
      '%) are unreported DAMbgs\n', sep='')

  return(invisible(list(
    multipleDetected = multipleDetectedDAMbgs,
    multipleKNOWN = multipleKNOWN,
    multipleUnreported = multipleUnreported
  )))
}

# Analyzes and visualizes expression patterns
# of the hosting gene in cell lines arbouring the DAMs.
plot_DAM_expression <- function(allHits, figuresPath, produce_plots = TRUE) {
  
  # FIGURE 2B
  if (produce_plots) {
    (file.path(figuresPath,'DAMs_exp_percentiles.'),5,5)
    hist(allHits$percBasalEXP_of_ps_cl,border=FALSE,col='darkcyan',
      main=paste('Basal expression percentile of the hosting gene\nin cell line(s) arbouring the DAMs'),
      xlab='-th')
    dev.off()
  }

  availableExpValues <- length(
    which(!is.na(allHits$percBasalEXP_of_ps_cl) &
            allHits$percBasalEXP_of_ps_cl > -Inf)
  )
  
  # percentage of expressed DAMs
  pp <- round(
    100*length(which(allHits$avgBasalEXP_fpkm_in_ps_cl>=1)) /
      availableExpValues, 2
  )
  
  # FIGURE 2B
  if (produce_plots) {
    (file.path(figuresPath,'percExpressedDAMs.'),5,5)
    pie(c(pp,100-pp),
        col=c('darkcyan','gray'),
        border=FALSE,
        labels = c('Expressed','Not expressed'))
    dev.off()
  }

  cat(pp, "% of DAMs are in genes expressed in the cell line in which the DAM is observed\n")
  cat(round(100*length(which(allHits$percBasalEXP_of_ps_cl > 50)) / availableExpValues, 2),
      "% of DAMs are in genes with a basal expression over the 50th percentile\n")
  cat(round(100*length(which(allHits$percBasalEXP_of_ps_cl > 80)) / availableExpValues, 2),
      "% of DAMs are in genes with a basal expression over the 80th percentile\n")
  cat(round(100*length(which(allHits$percBasalEXP_of_ps_cl == 100)) / availableExpValues, 2),
      "% of DAMs are in genes at the highest percentile of basal expression\n")
  
  return(list(
    pp = pp,
    availableExpValues = availableExpValues
  ))
}

# Computes overlap between significant hits and essential genes.
plot_essentiality_matching <- function(allHits, figuresPath, produce_plots = TRUE) {
  
  # FIGURE 2D
  if (produce_plots) {
    (file.path(figuresPath,'essentialityDAMmatching.'),5,5)
    hist(100*allHits$matching,border=FALSE,col='gray',main='')
    dev.off()
  }

  match_rate <- length(which(allHits$matching==1)) /
    length(allHits$matching)
  
  # FIGURE 2D
  if (produce_plots) {
    (file.path(figuresPath,'essentialityDAMmatchingHGp.'),5,6)
    plot(-log10(allHits$hypTest_p),bg=adjustcolor("blue", alpha.f = 0.3),col=NA,pch=21,cex=2,frame.plot=FALSE,ylab='-log10(HG p)')
    abline(h= -log10(0.05),lty=2)
    dev.off()
  }

  cat("Perfect overlap between DAM-hosting CCLs and gene-dependent CCLs observed in", round(100*match_rate, 2), "% of hits\n")

  return(match_rate)
}

# Creates a ranked summary plot of all significant hits
# by median fitness effect, color-coded by cancer type.
plot_all_hits_summary <- function(allHits, clc, figuresPath, produce_plots = TRUE) {
  
  set.seed(123)

  df <- data.frame(
    x = 1:nrow(allHits),
    y = allHits$medFitEff,
    label = paste(allHits$GENE,allHits$var)
  )

  df$fill <- adjustcolor(clc[allHits$ctype, 1], alpha.f = 0.6)
  
  # top 50 hits
  top_points <- df[order(df$y), ][seq_len(min(50, nrow(df))), ]

  # to reduce overlap between labels
  top_points <- top_points[order(top_points$x), ]

  n_label_rows <- 5
  label_base <- quantile(top_points$y, 0.25, na.rm = TRUE)
  label_step <- 0.08

  top_points$label_x <- top_points$x
  top_points$label_y <- top_points$y - 0.10 - label_step * ((seq_len(nrow(top_points)) - 1) %% n_label_rows)

  # extract unique cancer types
  utype <- unique(allHits$ctype)
  xpos  <- match(utype, allHits$ctype)

  cancer_label_y <- max(df$y, na.rm = TRUE) + 0.18

  manual_labels <- data.frame(
    x = xpos,
    y = rep(cancer_label_y, length(xpos)),
    label = utype,
    colours = clc[unique(allHits$ctype), 1]
  )
  
  # FIGURE 2A
  if (produce_plots) {
    (file.path(figuresPath,'All_Hits_summary.'),18,8)

    print(ggplot(df, aes(x = x, y = y)) +
      geom_point(
        aes(fill = fill),
        shape = 21,
        colour = "black",
        size = 3,
        stroke = 0.5
      ) +
      geom_segment(
        data = top_points,
        aes(x = x, y = y, xend = label_x, yend = label_y),
        color = "grey50",
        linewidth = 0.25
      ) +
      geom_text(
        data = top_points,
        aes(x = label_x, y = label_y, label = label),
        size = 2.4,
        angle = 60,
        hjust = 1,
        vjust = 1
      ) + 
      geom_text(
        data = manual_labels,
        aes(x = x, y = y, label = label, color=colours),
        size = 2.4,
        fontface = "italic",
        hjust = 0,
        vjust = 0.5,
        angle = 70
      ) +
      scale_fill_identity() +
      scale_color_identity() +
      coord_cartesian(
        ylim = c(
          min(top_points$label_y, na.rm = TRUE) - 0.1,
          cancer_label_y + 0.35
        ),
        clip = "off"
      ) +
      theme_minimal() +
      theme(
        plot.margin = margin(t = 40, r = 30, b = 30, l = 30)
      ))

    invisible(dev.off())
  }
}

# Computes correlations between hits, cohort size,
# and mutation burden across cancer types.
# Generates multi-panel plots and returns SNR values.
plot_hits_correlations <- function(allDAMs, totalTestedVariants, incl_cl_annot, clc, figuresPath, produce_plots = TRUE) {
  
  # number of hits per cancer type
  nhits_across_ctypes <- summary(as.factor(allDAMs$ctype))
  # number of tested variants per cancer type
  ntested_variants_across_ctypes <- summary(as.factor(totalTestedVariants$cty))
  # number of analyzed cell lines per cancer type
  ncellLines_across_ctypes <- summary(as.factor(incl_cl_annot$cancer_type))

  commoncl <- intersect(names(nhits_across_ctypes), names(ncellLines_across_ctypes))
  
  # signal-to-noise ratio
  percHitsPerCls <- 100*(nhits_across_ctypes[commoncl]/ncellLines_across_ctypes[commoncl]) /
    (ntested_variants_across_ctypes[commoncl]/ncellLines_across_ctypes[commoncl])

  # SUPPLEMENTARY FIGURE 3BCD
  if (produce_plots && length(commoncl) >= 3) {
    pdf(file.path(figuresPath,'nHits_correlations.pdf'),7,7)
    par(mfrow=c(2,2))

    # n cell lines vs hits
    plot(ncellLines_across_ctypes[commoncl],
        nhits_across_ctypes[commoncl],
        bg=clc[commoncl,1],pch=21,
        xlab='n analysed cell lines',
        ylab='n hits',
        cex=1.5,frame.plot=FALSE)

    tt <- cor.test(ncellLines_across_ctypes[commoncl],nhits_across_ctypes[commoncl])$p.value

    legend('topleft',cex=0.8,c(
      paste('R = ',round(cor(ncellLines_across_ctypes[commoncl],nhits_across_ctypes[commoncl]),2)),
      paste('p = ',round(tt,3))))

    # tested variants vs hits
    plot(ntested_variants_across_ctypes[commoncl],
        nhits_across_ctypes[commoncl],
        bg=clc[commoncl,1],pch=21,
        xlab='n tested variants',
        ylab='n hits',
        cex=1.5,frame.plot=FALSE)

    tt <- cor.test(ntested_variants_across_ctypes[commoncl],nhits_across_ctypes[commoncl])$p.value

    legend('bottomright',cex=0.8,c(
          paste('R = ',round(cor(ntested_variants_across_ctypes[commoncl],nhits_across_ctypes[commoncl]),2)),
          paste('p = ',format(tt,scientific=TRUE,digits=3))
          ))

    tt <- cor.test(ntested_variants_across_ctypes[commoncl],percHitsPerCls[commoncl])$p.value

    plot(ntested_variants_across_ctypes[commoncl],
        percHitsPerCls[commoncl],
        bg=clc[names(ntested_variants_across_ctypes),1],
        pch=21,cex=1.5,
        xlab='n. tested variants',
        ylab='Signal to Noise Ratio',
        frame.plot = FALSE)

    legend('topright',cex=0.8,c(
      paste('R = ',round(cor(ntested_variants_across_ctypes[commoncl],percHitsPerCls[commoncl]),2)),
      paste('p = ',format(tt,scientific=TRUE,digits=3))
    ))

    plot(0,0,type='n',frame.plot=FALSE,xaxt='n',yaxt='n',xlab='',ylab='')
    legend('center',sort(commoncl),
          fill=clc[sort(commoncl),1],cex=0.2)

    dev.off()
  }

  return(sort(percHitsPerCls,decreasing=TRUE))
}

# Stratifies cancer types into genomic instability groups,
# annotates SNR plot, generates boxplot comparison,
# and performs pairwise t-tests between groups.
plot_SNR_genomic_groups <- function(percHitsPerCls, clc, figuresPath, produce_plots = TRUE) {

  percHitsPerCls <- sort(percHitsPerCls, decreasing = TRUE)

  # define groups 
  genomicallyQuite <- c(
    "Acute Myeloid Leukemia",
    "B-Cell Non-Hodgkin's Lymphoma",
    "B-Lymphoblastic Leukemia",
    "Biliary Tract Carcinoma",
    "Burkitt's Lymphoma",
    "Chronic Myelogenous Leukemia",
    "Ewing's Sarcoma",
    "Mesothelioma",
    "Neuroblastoma",
    "Prostate Carcinoma",
    "Rhabdomyosarcoma",
    "T-Cell Non-Hodgkin's Lymphoma",
    "T-Lymphoblastic Leukemia",
    "Thyroid Gland Carcinoma"
  )

  genomicallyIntermediate <- c(
    "Bladder Carcinoma",
    "Breast Carcinoma",
    "Cervical Carcinoma",
    "Esophageal Squamous Cell Carcinoma",
    "Glioma",
    "Head and Neck Carcinoma",
    "Kidney Carcinoma",
    "Oral Cavity Carcinoma",
    "Osteosarcoma",
    "Plasma Cell Myeloma",
    "Squamous Cell Lung Carcinoma"
  )

  genomicallyNoisy <- c(
    "Colorectal Carcinoma",
    "Endometrial Carcinoma",
    "Esophageal Carcinoma",
    "Gastric Carcinoma",
    "Glioblastoma",
    "Hepatocellular Carcinoma",
    "Melanoma",
    "Non-Small Cell Lung Carcinoma",
    "Ovarian Carcinoma",
    "Pancreatic Carcinoma",
    "Small Cell Lung Carcinoma"
  )

  # keep only cancer types included in the current analysis
  genomicallyQuite <- intersect(genomicallyQuite, names(percHitsPerCls))
  genomicallyIntermediate <- intersect(genomicallyIntermediate, names(percHitsPerCls))
  genomicallyNoisy <- intersect(genomicallyNoisy, names(percHitsPerCls))
  
  if (produce_plots) {
    # FIGURE 2F
    pdf(file.path(figuresPath,'SNR_across_ctypes.pdf'),9,7)

    par(mar=c(16,4,2,0))

    dd <- barplot(percHitsPerCls, col = clc[names(percHitsPerCls),1], las = 2, border = FALSE,
                  ylab = 'Avg % of DAMs per cell line', ylim = c(0,0.9))

    # annotated SNR plot 
    par(xpd = TRUE)

    points(dd[sort(match(genomicallyQuite, names(percHitsPerCls)))],
          rep(-0.03, length(genomicallyQuite)),
          pch = 16, col = 'darkgray', cex = 1.8)

    points(dd[sort(match(genomicallyIntermediate, names(percHitsPerCls)))],
          rep(-0.035, length(genomicallyIntermediate)),
          pch = 17, col = 'darkgray', cex = 1.5)

    points(dd[sort(match(genomicallyNoisy, names(percHitsPerCls)))],
          rep(-0.035, length(genomicallyNoisy)),
          pch = 18, col = 'darkgray', cex = 1.8)

    dev.off()

    # boxplot
    # FIGURE 2G
    pdf(file.path(figuresPath, 'SNR_t_tests.pdf'), 7, 8)
    par(mar = c(15, 5, 2, 2))
    cex.axis = 0.8

    boxplot(percHitsPerCls[genomicallyQuite],
            percHitsPerCls[genomicallyIntermediate],
            percHitsPerCls[genomicallyNoisy],
            names = c(paste('genomically quiet cancers\n(few mutations; driven by\nfusions or copy number)'),
                      'genomically intermediate cancers',
                      'genomically noisy cancers'),las=2,frame.plot=FALSE)
    dev.off()
  }

  # statistical tests
  cat("\n--- SNR Group Comparisons ---\n")

  if (length(genomicallyQuite) >= 2 && length(genomicallyNoisy) >= 2) {
    cat("\nGenomically quiet vs Noisy t.test\n")
    print(t.test(
      percHitsPerCls[genomicallyQuite],
      percHitsPerCls[genomicallyNoisy]
    ))
  }

  if (length(c(genomicallyQuite, genomicallyIntermediate)) >= 2 && length(genomicallyNoisy) >= 2) {
    cat("\nGenomically quiet + intermediate vs Noisy t.test\n")
    print(t.test(
      percHitsPerCls[c(genomicallyQuite, genomicallyIntermediate)],
      percHitsPerCls[genomicallyNoisy]
    ))
  }

  if (length(genomicallyIntermediate) >= 2 && length(genomicallyNoisy) >= 2) {
    cat("\nGenomically intermediate vs Noisy t.test\n")
    print(t.test(
      percHitsPerCls[genomicallyIntermediate],
      percHitsPerCls[genomicallyNoisy]
    ))
  }

  if (length(genomicallyQuite) >= 2 && length(genomicallyIntermediate) >= 2) {
    cat("\nGenomically quiet vs Intermediate t.test\n")
    print(t.test(
      percHitsPerCls[genomicallyQuite],
      percHitsPerCls[genomicallyIntermediate]
    ))
  }

  if (length(genomicallyQuite) >= 2 && length(c(genomicallyIntermediate, genomicallyNoisy)) >= 2) {
    cat("\nGenomically quiet vs Intermediate + Noisy t.test\n")
    print(t.test(
      percHitsPerCls[genomicallyQuite],
      percHitsPerCls[c(genomicallyIntermediate, genomicallyNoisy)]
    ))
  }

  return(invisible(list(
    quiet = genomicallyQuite,
    intermediate = genomicallyIntermediate,
    noisy = genomicallyNoisy
  )))
}

# Identifies recurrent DAMbgs not present in the
# reference IntOGen driver list, plots those appearing in >2 analyses,
# and returns their frequency counts.
plot_most_frequent_unreported_DAMbgs <- function(allHits, intogen_drivers, figuresPath, produce_plots = TRUE) {

  mostFreqDAMbgs <- sort(summary(as.factor(allHits$GENE),10000), decreasing=TRUE)

  mostFreqDAMbgs <- mostFreqDAMbgs[
    setdiff(names(mostFreqDAMbgs), intogen_drivers)
  ]

  tmp <- mostFreqDAMbgs[which(mostFreqDAMbgs>2)]
  
  # SUPPLEMENTARY FIGURE 5A
  if (produce_plots && length(tmp) > 0) {
    pdf(file.path(figuresPath,'mostFrequentUnreportedDAMbgs.pdf'),12,5)
    barplot(tmp,las=2,border=FALSE,col='orange',
            ylab='Unreported DAMbgs in n cancer type specific analyses',
            cex.names = 0.7)
    invisible(dev.off())
  }

  return(invisible(tmp))
}


