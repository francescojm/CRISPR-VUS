# -------- INTOGEN DRIVER ANNOTATION -------- #

# Annotates gene-level dependency results with IntOGen driver information.
# For each gene in RESTOT, retrieves associated cancer types and driver roles
# (Activating, Loss-of-Function, ambiguous) from the IntOGen database,
# computes role percentages, and appends these annotations to the results table.
annotate_IntOGen <- function(RESTOT, intOGen_drivers) {
  
  intOGen_annot <- do.call(rbind, lapply(RESTOT$GENE, function(x) {
    ids <- which(intOGen_drivers$SYMBOL == x)
    
    # extract cancer types
    intOGen_driver_for <- paste(sort(unique(intOGen_drivers$CANCER_TYPE[ids])), collapse = ' | ')
    if (intOGen_driver_for == "") {
        intOGen_driver_for <- 'none'
    }
    
    # extract roles 
    intOGen_Role_Perc <- c(Act=0, LoF=0, ambiguous=0)
    rolSums <- summary(as.factor(intOGen_drivers$ROLE[ids]))
    intOGen_Role_Perc[names(rolSums)] <- round(100 * rolSums / sum(rolSums), 2)
    
    data.frame(intOGen_driver_for=intOGen_driver_for,Act=intOGen_Role_Perc[1],LoF=intOGen_Role_Perc[2],ambigous=intOGen_Role_Perc[3],row.names = x)
  }))
  
  RESTOT <- cbind(RESTOT, intOGen_annot)
  
  return(RESTOT)
}