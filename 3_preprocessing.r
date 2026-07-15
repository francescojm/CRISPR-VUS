## -------- PREPROCESSING -------- ##

# Preprocesses the gene dependency matrix for downstream analysis.
# Steps:
# - Aligns cell lines to CMP annotation using BROAD_ID and replaces row names
#   with model_id, removing unmatched models
# - Transposes the matrix and removes genes with missing values
# - Binarizes dependency scores using the specified dependency threshold
#   (default: <= -0.5 classified as dependent = 1, otherwise 0)
preprocess_crispr_matrix <- function(scaled_depFC, CMP_annot0, cl_variants, threshold = -0.5) {

    # summary of available data
    cat(nrow(CMP_annot0), "annotated models in the Cell Models Passports\n")
    cat("of which", length(which(is.element(CMP_annot0$model_id, cl_variants$model_id))), "with mutation data\n")

    # align to CMP annotation
    toremove <- which(is.na(CMP_annot0$model_id[match(rownames(scaled_depFC),CMP_annot0$BROAD_ID)]))
     
    scaled_depFC <- scaled_depFC[-toremove,]
    rownames(scaled_depFC) <- CMP_annot0$model_id[match(rownames(scaled_depFC),CMP_annot0$BROAD_ID)]

    # transpose 
    scaled_depFC <- t(scaled_depFC)

    # remove genes with missing values
    scaled_depFC <- scaled_depFC[-which(rowSums(is.na(scaled_depFC))>0),]

    # select only cell lines with data in both CRISPR screens and CMP annotation file
    CMP_annot <- CMP_annot0[which(is.element(CMP_annot0$model_id, colnames(scaled_depFC))), ]

    cat("of which", nrow(CMP_annot), "with high quality CRISPR data\n")

    # binarize
    bdep <- scaled_depFC
    bdep[scaled_depFC <= threshold] <- 1
    bdep[scaled_depFC > threshold] <- 0  

    return(list(
        scaled_depFC = scaled_depFC,
        bdep = bdep,
        CMP_annot = CMP_annot
    ))
}
