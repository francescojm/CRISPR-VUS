# -------- STRING ANALYSIS -------- #

suppressPackageStartupMessages({
  library(STRINGdb)
  library(tidyverse)
  library(ggplot2)
})

# compute STRING interaction score between non-driver hit genes and known driver genes
compute_STRING_score <- function(hits_nodriver, driver_genes, string_db) {
  
  df <- data.frame(gene = unique(c(hits_nodriver, driver_genes)))
  Genesmapped <- string_db$map(df, "gene", removeUnmappedRows = TRUE)

  # Retrieve all protein-protein interactions among the mapped STRING IDs
  int <- string_db$get_interactions(Genesmapped$STRING_id)
  
  # Compute the interaction score between hit genes and driver genes
  score <- sum(
    int[int[,1] %in% Genesmapped[Genesmapped[,1] %in% hits_nodriver, 2] &
          int[,2] %in% Genesmapped[Genesmapped[,1] %in% driver_genes, 2], 3]
  ) +
    sum(
      int[int[,2] %in% Genesmapped[Genesmapped[,1] %in% hits_nodriver, 2] &
            int[,1] %in% Genesmapped[Genesmapped[,1] %in% driver_genes, 2], 3]
    )
  
  return(list(score = score, Genesmapped = Genesmapped, interactions = int))
}

# compute empirical p-value by random sampling
compute_empirical_p <- function(hits_nodriver, driver_genes, all_genes_nodriver,
                                string_db, score, n_random = 1000, figuresPath = NULL, produce_plots = TRUE) {
  
  set.seed(84905750)
  scorerand <- c()
  
  for(j in 1:n_random) {
    # Randomly sample genes from the background gene set, build df, and map IDs
    rand_genes <- sample(all_genes_nodriver, length(hits_nodriver), replace = FALSE)
    dfrand <- data.frame(gene = unique(c(rand_genes, driver_genes)))
    Randmapped <- string_db$map(dfrand, "gene", removeUnmappedRows = TRUE)

    # Retrieve interactions among the mapped STRING IDs
    intrand <- string_db$get_interactions(Randmapped$STRING_id)
    
    # Compute the interaction score between random genes and driver genes
    scorerand <- c(scorerand, sum(
      intrand[intrand[,1] %in% Randmapped[Randmapped[,1] %in% rand_genes,2] &
               intrand[,2] %in% Randmapped[Randmapped[,1] %in% driver_genes,2],3]
    ) +
      sum(
        intrand[intrand[,2] %in% Randmapped[Randmapped[,1] %in% rand_genes,2] &
                 intrand[,1] %in% Randmapped[Randmapped[,1] %in% driver_genes,2],3]
      ))
  }
  
  # Compute empirical p-value
  empP <- length(which(scorerand >= score)) / n_random
  
  # FIGURE 3I
  if (produce_plots) {
    # Plot density
    df_plot <- data.frame(First_neighbours = scorerand)
    pdf(file.path(figuresPath, "STRING_PPI_connections_emp.pdf"),5,4)
    print(ggplot(df_plot, aes(x = First_neighbours)) +
      geom_density() +
      geom_vline(xintercept = score, linetype = "dashed", color = "red", size = 1) +
      theme_classic())
    dev.off()
  }
  
  return(list(empP = empP, scorerand = scorerand))
}

# run the full STRING interaction analysis workflow
run_STRING_analysis <- function(allHits, driver_genes, totalTestedVariants, string_db, figuresPath = NULL, n_random = 1000, produce_plots = TRUE) {
  
  hits <- unique(allHits$GENE)

  # Remove known driver genes from the hits
  hits_nodriver <- setdiff(hits, driver_genes)

  # Define background genes excluding driver genes
  all_genes_nodriver <- setdiff(unique(totalTestedVariants$gene_symbol), driver_genes)
  
  # compute observed STRING score
  res <- compute_STRING_score(hits_nodriver, driver_genes, string_db)
  score <- res$score
  
  # compute empirical p-value
  emp_res <- compute_empirical_p(hits_nodriver, driver_genes, all_genes_nodriver,
                                 string_db, score = score, n_random = n_random, figuresPath = figuresPath, produce_plots = produce_plots)
  empP <- emp_res$empP

  #cat("Empirical p-value:", empP, "\n")
  
  return(list(score = score, empP = empP, scorerand = emp_res$scorerand))
}