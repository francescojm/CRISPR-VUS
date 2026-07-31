## -------- TISSUE SELECTION -------- ##

# Filters CMP cell line annotations to retain well-represented and informative
# cancer types, and saves the filtered annotation table.
# Steps:
# - Counts the number of cell lines per cancer_type and removes tissue types represented by fewer than min_n (default = 5) cell lines
# - Excludes generic or uninformative cancer categories
#   (e.g., "Other Solid Cancers", "Non-Cancerous")
select_and_save_tissues <- function(CMP_annot, min_n = 5, sel_tissue_idx = NULL) {

    tissues <- CMP_annot$cancer_type

    # remove tissues that appear fewer than min_n times
    st <- summary(as.factor(tissues))
    tissues <- sort(setdiff(tissues, names(which(st < min_n))))

    # remove generic/uninformative categories
    final_tissues <- setdiff(tissues, c('Other Solid Carcinomas','Other Solid Cancers','Other Sarcomas','Other Blood Cancers', 'Non-Cancerous'))

    # restrict to selected tissue indices, if provided
    if (!is.null(sel_tissue_idx)) {
      final_tissues <- final_tissues[sel_tissue_idx]
    }

    CMP_annot <- CMP_annot[which(CMP_annot$cancer_type %in% final_tissues), ]

    return(list(
      CMP_annot = CMP_annot,
      tissues = final_tissues
    ))
}
