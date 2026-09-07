# -------- ADD SIFT AND POLYPHEN SCORES -------- #

library(data.table)
library(stringr)
library(biomaRt)
library(httr)
library(jsonlite)

# Filters variant data to retain only those matching DAMs,
# removes duplicates, orders by gene, and formats a clean table with genomic
# and mutation information.
filter_and_match_DAMs_v2 <- function(cl_variants0, allDAMs) {
  
  cl_vars <- cl_variants0[which(!is.na(cl_variants0$position)), ]
  
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

  return(cl_vars)
}

# Prepare variant coordinates and strand information for downstream VEP annotation
annotate_DAM_positions_v2 <- function(cl_vars) {
  
  dt <- as.data.table(cl_vars)
  
  # seq_region_name
  dt[, seq_region_name := str_replace(chromosome, "^chr", "")]
  
  # start/end positions 
  dt[, `:=`(start = as.integer(position), end   = as.integer(position))]
  
  # allele_string from REF/ALT
  dt[, allele_string := ifelse(!is.na(reference) & !is.na(alternative),
                             paste0(reference, "/", alternative), NA_character_)]

  # strand via Ensembl 
  ens_ids <- unique(na.omit(ifelse(!is.na(dt$ensembl_gene_id), dt$ensembl_gene_id, dt$gene_id)))
  
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  strand_map <- getBM(
    attributes = c("ensembl_gene_id","strand"),
    filters    = "ensembl_gene_id",
    values     = ens_ids,
    mart       = mart
  )
  strand_map <- as.data.table(strand_map)

  # Join strand annotation back to the input table
  dt <- strand_map[dt, on = c("ensembl_gene_id" = "ensembl_gene_id")]
  setnames(dt, "strand", "strand_numeric")
  dt[, strand := as.integer(strand_numeric)][, strand_numeric := NULL]

  out <- dt[, .(gene_symbol, model_id, seq_region_name, start, end, strand, allele_string)]

  # Remove intermediate columns
  dt <- dt[,-c(4,5,13,14)]
  
  return(dt)
}

# Query Ensembl VEP in batches and extract SIFT/PolyPhen predictions per variant
vep_sift_polyphen <- function(df, batch_size  = 50, max_retries = 3, sleep_sec = 1, verbose = TRUE) {

  # Check that required columns are present
  req_cols <- c("seq_region_name","start","end","reference","alternative")
  stopifnot(all(req_cols %in% names(df)))

  x <- as.data.table(df)
  # Keep only rows with complete variant information
  x <- x[!is.na(seq_region_name) & !is.na(start) & !is.na(end) &
           !is.na(reference) & !is.na(alternative)]
  x[, `:=`(start = as.integer(start), end = as.integer(end))]
  x <- x[start > 0 & end >= start]

  if (!nrow(x)) {
    warning("No valid variants after basic filtering.")
    return(list(results = data.table(), failed = data.table()))
  }

  # Build VEP input strings
  mk_variant <- function(chr, st, en, ref, alt) sprintf("%s %s %s %s/%s", chr, st, en, ref, alt)
  x[, vep_input := mk_variant(seq_region_name, start, end, reference, alternative)]

  # Split input into batches
  batches <- split(x, ceiling(seq_len(nrow(x))/batch_size))
  n_batches <- length(batches)

  if (verbose) {
    message(sprintf("Submitting %d variants in %d batch(es) to Ensembl VEP REST…",
                    nrow(x), n_batches))
    pb <- utils::txtProgressBar(min = 0, max = n_batches, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  # Convert one VEP response record into a transcript-level table 
  parse_one <- function(rec) {
    if (is.null(rec$transcript_consequences)) return(NULL)
    tc <- rec$transcript_consequences
    dt <- rbindlist(lapply(tc, as.data.frame), fill = TRUE)

    dt[, `:=`(input = rec$input,
              seq_region_name = rec$seq_region_name,
              start = rec$start,
              end = rec$end,
              allele_string = rec$allele_string)]
    dt
  }

  # Submit one batch to VEP
  request_batch <- function(dt_batch, depth = 0, batch_id = NA_integer_) {
    if (nrow(dt_batch) == 0) return(list(ok = NULL, bad = NULL))

    variants <- as.list(dt_batch$vep_input)
    attempt <- 0

    repeat {
      attempt <- attempt + 1
      if (verbose) {
        msg <- sprintf("Batch %s | size=%d | attempt=%d%s",
                       ifelse(is.na(batch_id), "?", batch_id),
                       nrow(dt_batch), attempt,
                       if (depth > 0) sprintf(" | depth=%d (bisect)", depth) else "")
        message(msg)
      }

      res <- try({
        POST(
          url = "https://rest.ensembl.org/vep/human/region",
          add_headers(`Content-Type` = "application/json", `Accept` = "application/json"),
          body = toJSON(list(variants = variants), auto_unbox = TRUE),
          timeout(60)
        )
      }, silent = TRUE)

      # Retry on connection failure
      if (inherits(res, "try-error")) {
        if (verbose) message("  -> HTTP request failed (connection).")
        if (attempt <= max_retries) { Sys.sleep(sleep_sec * attempt); next } else break
      }

      code <- status_code(res)
      
      # Successful response
      if (code == 200) {
        ans <- httr::content(res, as = "parsed", type = "application/json")
        ok  <- rbindlist(lapply(ans, parse_one), fill = TRUE)
        return(list(ok = ok, bad = NULL))
      }

      # Retry on rate-limiting or temporary server unavailability
      if (code %in% c(429, 503)) {
        if (verbose) message(sprintf("  -> Server says %d; backing off…", code))
        if (attempt <= max_retries) { Sys.sleep(sleep_sec * attempt); next } else break
      }

      # Split batch to isolate problematic records
      if (code >= 500) {
        if (verbose) message(sprintf("  -> HTTP %d; bisecting batch to isolate offending records…", code))
        if (nrow(dt_batch) == 1 || depth > 10) {
          if (verbose) message("  -> Reached single record or max depth; marking as failed.")
          return(list(ok = NULL, bad = dt_batch))
        }
        mid   <- nrow(dt_batch) %/% 2
        left  <- request_batch(dt_batch[1:mid], depth + 1, batch_id)
        right <- request_batch(dt_batch[(mid+1):nrow(dt_batch)], depth + 1, batch_id)
        return(list(ok = rbindlist(list(left$ok, right$ok), fill = TRUE),
                    bad = rbindlist(list(left$bad, right$bad), fill = TRUE)))
      }

      # Other client error, mark whole batch as failed
      if (verbose) message(sprintf("  -> HTTP %d (client error); marking batch as failed.", code))
      return(list(ok = NULL, bad = dt_batch))
    }

    # Exit without success
    list(ok = NULL, bad = dt_batch)
  }

  # Process batches 
  ok_list <- list(); bad_list <- list()
  for (i in seq_along(batches)) {
    ans <- request_batch(batches[[i]], batch_id = i)
    if (!is.null(ans$ok))  ok_list[[length(ok_list)+1]]  <- ans$ok
    if (!is.null(ans$bad)) bad_list[[length(bad_list)+1]] <- ans$bad
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) message("\nParsing and summarising transcript-level predictions…")
  resdt <- rbindlist(ok_list, fill = TRUE)

  # Harmonise VEP output fields and summarise scores across transcripts
  if (nrow(resdt)) {
    # Harmonise gene symbol naming
    if (!"gene_symbol" %in% names(resdt) && "hgnc_symbol" %in% names(resdt)) {
      resdt[, gene_symbol := hgnc_symbol]
    }
    # Standardise consequence column naming
    if (!"consequence_terms" %in% names(resdt)) {
      if ("consequence_term" %in% names(resdt)) {
        resdt[, consequence_terms := consequence_term]
      } else if ("consequence" %in% names(resdt)) {
        resdt[, consequence_terms := consequence]
      } else {
        resdt[, consequence_terms := NA_character_]
      }
    }
    # Ensure expected columns exist
    for (cn in c("sift_score","sift_prediction","polyphen_score","polyphen_prediction",
                 "gene_symbol","transcript_id")) {
      if (!cn %in% names(resdt)) resdt[[cn]] <- NA
    }

    # Build join key and merge back to inputs
    resdt[, join_key := sprintf("%s %s %s %s", seq_region_name, start, end, allele_string)]
    x[, join_key := sprintf("%s %s %s %s/%s", seq_region_name, start, end, reference, alternative)]

    keep <- resdt[, .(join_key, gene_symbol, transcript_id, consequence_terms,
                      sift_score, sift_prediction, polyphen_score, polyphen_prediction)]

    suppressWarnings({
      keep[, sift_score := as.numeric(sift_score)]
      keep[, polyphen_score := as.numeric(polyphen_score)]
    })

    # Summarise transcript-level predictions into one row per variant
    summary <- keep[
      , .(sift_min = if (all(is.na(sift_score))) NA_real_ else min(sift_score, na.rm = TRUE),
          sift_pred_any = paste(na.omit(unique(sift_prediction)), collapse = "|"),
          polyphen_max = if (all(is.na(polyphen_score))) NA_real_ else max(polyphen_score, na.rm = TRUE),
          polyphen_pred_any = paste(na.omit(unique(polyphen_prediction)), collapse = "|")),
      by = join_key
    ]

    annotated <- merge(x, summary, by = "join_key", all.x = TRUE)
  } else {
    annotated <- x
    annotated[, c("sift_min","sift_pred_any","polyphen_max","polyphen_pred_any") := NA]
  }

  failed <- rbindlist(bad_list, fill = TRUE)

  if (verbose) {
    message(sprintf("Done. Annotated: %d | Failed: %d", nrow(annotated), nrow(failed)))
    if (nrow(failed)) {
      message("A few examples of failed inputs:")
      print(head(failed[, .(seq_region_name, start, end, reference, alternative)], 5))
      message("Common causes: wrong build, REF/ALT mismatch, large/complex indels.")
    }
  }

  list(
    results = annotated[, !c("join_key","vep_input")],
    failed  = failed
  )
}

# Merge VEP-derived SIFT and PolyPhen annotations back into the DAM table
merge_DAMs_with_VEP <- function(allDAMs, annotated_results) {
  
  # Create matching signatures 
  allDAMs_sigs <- paste(allDAMs$GENE, allDAMs$var)
  eff_sigs <- paste(annotated_results$gene_symbol,
                    annotated_results$protein_mutation)
  ii <- match(allDAMs_sigs, eff_sigs)
  
  # Append summarised prediction columns to the original DAM table
  allDAMs_with_scores <- cbind(allDAMs, annotated_results[ii, c("sift_min", "sift_pred_any", "polyphen_max", "polyphen_pred_any")])
  
  return(allDAMs_with_scores)
}

# Return fallback value
`%||%` <- function(a,b) if (is.null(a) || is.na(a)) b else a
# Check whether a prediction label is present in a text field
has_label <- function(x, lab) grepl(sprintf("\\b%s\\b", lab), x %||% "", ignore.case = TRUE)

# Flag SIFT scores consistent with deleterious effect
is_sift_del <- function(s) !is.na(s) & s < 0.05
# Flag SIFT scores consistent with tolerated effect
is_sift_tol <- function(s) !is.na(s) & s >= 0.05

# Flag PolyPhen scores in the probably damaging range
is_pp_prob <- function(p) !is.na(p) & p >= 0.909
# Flag PolyPhen scores in the possibly damaging range
is_pp_poss <- function(p) !is.na(p) & p >= 0.446 & p < 0.909
# Flag PolyPhen scores in the benign range
is_pp_benign <- function(p) !is.na(p) & p <= 0.445

# Classify variants into broad impact groups using SIFT and PolyPhen evidence
classify_variant <- function(sift_min, sift_pred_any, polyphen_max, polyphen_pred_any) {
  # Numeric evidence
  s_del <- is_sift_del(sift_min)
  s_tol <- is_sift_tol(sift_min)

  p_prob <- is_pp_prob(polyphen_max)
  p_poss <- is_pp_poss(polyphen_max)
  p_ben  <- is_pp_benign(polyphen_max)

  # text evidence 
  s_del_l <- has_label(sift_pred_any, "deleterious")
  s_tol_l <- has_label(sift_pred_any, "tolerated")

  p_prob_l <- has_label(polyphen_pred_any, "probably")
  p_poss_l <- has_label(polyphen_pred_any, "possibly")
  p_ben_l  <- has_label(polyphen_pred_any, "benign")

  # collapse signals
  any_del <- s_del || s_del_l || p_prob || p_prob_l || p_poss || p_poss_l
  any_ben <- s_tol || s_tol_l || p_ben || p_ben_l

  # No usable prediction data available
  if (all(is.na(sift_min), is.na(polyphen_max)) &&
      !any(nzchar(sift_pred_any), nzchar(polyphen_pred_any))) {
    return("Unknown impact")
  }

  # high / moderate / possible impact strata
  if ((s_del || s_del_l) && (p_prob || p_prob_l)) return("High impact")
  if ((s_del || s_del_l) && (p_poss || p_poss_l)) return("Moderate impact")
  if ((p_prob || p_prob_l) && !(s_del || s_del_l)) return("Moderate impact")
  if ((s_del || s_del_l) && !(p_prob || p_prob_l)) return("Possible impact")
  if ((p_poss || p_poss_l) && !(s_del || s_del_l)) return("Possible impact")

  # low: at least one benign/tolerated and no damaging evidence
  if (any_ben && !any_del) return("Low impact")

  # Mixed or unclear evidence
  "Unknown impact"
}

# Plot the distribution of VEP impact classes across DAMs and return class counts
plot_DAM_VEP_summary <- function(allDAMs_with_scores, figuresPath, driver_genes = NULL, prefix = "DAM", produce_plots = TRUE) {
  
  data_use <- allDAMs_with_scores
  
  # Optionally remove known driver genes from summary
  if (!is.null(driver_genes)) {
    data_use <- data_use[
      !data_use$GENE %in% driver_genes, ]
  }
  
  # Count variants per impact category 
  VEPs <- summary(as.factor(data_use$impact_call))
  VEPs <- VEPs[c("High impact",
                 "Moderate impact",
                 "Possible impact",
                 "Low impact",
                 "Unknown impact")]
  
  VEPs[is.na(VEPs)] <- 0
  
  # SUPPLEMENTARY FIGURE 19 A-B
  if (produce_plots && sum(VEPs) > 0) {
    pdf(file.path(figuresPath, paste0(prefix, "_VEP_prediction.pdf")), 11, 6)
    pie(VEPs,
        col = c("#800026","#fc4e2a","#feb24c","#ffeda0","grey85"))
    dev.off()
  }

  cat(prefix, ":", sum(VEPs[1:2]), "DAMs have a high/moderate functional impact\n")
  cat(prefix, ":", sum(VEPs[1:3]), "DAMs have a high/moderate/possible functional impact\n")
  
  return(invisible(VEPs))
}

# Transfer DAM-based impact calls to SAMs by matching gene and variant identifiers
annotate_SAMs_with_VEP <- function(allSAMs, allDAMs_with_scores) {
  
  # Match SAM gene/variant pairs to annotated DAM gene/variant pairs
  ii <- match(
    paste(allSAMs$genes, allSAMs$vars),
    paste(allDAMs_with_scores$GENE,
          allDAMs_with_scores$var)
  )
  
  allSAMs_with_scores <- cbind(
    allSAMs,
    allDAMs_with_scores[ii, "impact_call"]
  )
  
  colnames(allSAMs_with_scores)[ncol(allSAMs_with_scores)] <- "impact_call"
  
  return(allSAMs_with_scores)
}

# Plot the distribution of VEP impact classes across unique SAMs
plot_SAM_VEP_summary <- function(allSAMs_with_scores, figuresPath, driver_genes = NULL, prefix = "SAM", produce_plots = TRUE) {
  
  # Create unique SAM signatures by cancer type, gene and variant
  samSigs <- paste(allSAMs_with_scores$cancer_type, allSAMs_with_scores$genes, allSAMs_with_scores$vars)
  
  usigs <- unique(samSigs)
  udata <- allSAMs_with_scores[match(usigs, samSigs), ]
  
  # Optionally exclude driver genes
  if (!is.null(driver_genes)) {
    udata <- udata[!udata$genes %in% driver_genes, ]
  }
  
  # Count unique SAMs per impact class
  VEPs <- summary(as.factor(udata$impact_call))
  VEPs <- VEPs[c('High impact','Moderate impact','Possible impact','Low impact','Unknown impact')]
  VEPs[is.na(VEPs)] <- 0
  
  # SUPPLEMENTARY FIGURE 19 C-D
  if (produce_plots && sum(VEPs) > 0) {
    pdf(file.path(figuresPath, paste0(prefix, "_VEP_prediction.pdf")), 11, 6)
    pie(VEPs,
        col = c("#800026","#fc4e2a","#feb24c","#ffeda0","grey85"))
    dev.off()
  }
  
  return(invisible(VEPs))
}

