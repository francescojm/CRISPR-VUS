## -------- DATA LOADING -------- ##

library(httr)
library(readr)
library(readxl)
suppressPackageStartupMessages(library(rvest))

# Downloads and prepares the latest versions of all required external datasets
# Specifically:
# - Downloads and decompresses Sanger gene annotation, model list, and mutation files
# - Detects and downloads the latest GDSC1 and GDSC2 fitted dose–response data
# - Identifies and extracts the most recent IntOGen Drivers release, saving the
#   Compendium_Cancer_Genes.tsv file
# - Retrieves the latest DepMap release via API and downloads the
#   CRISPRGeneEffect dataset
download_latest_inputs <- function(raw_dir) {

    options(timeout = 600)
  
    # -------- Sanger files (gzipped) -------- #
    gene_annot_url <- "https://cog.sanger.ac.uk/cmp/download/gene_identifiers_latest.csv.gz"
    model_list_url <- "https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz"
    cl_variants_url <- "https://cog.sanger.ac.uk/cmp/download/mutations_all_latest.csv.gz" 
  
    gene_annot_path <- file.path(raw_dir, "gene_identifiers_latest.csv")
    model_list_path <- file.path(raw_dir, "model_list_latest.csv")
    cl_variants_path <- file.path(raw_dir, "mutations_all_latest.csv") 
  
    # temporary paths for gzipped files
    gene_annot_gz <- file.path(raw_dir, "gene_identifiers_latest.csv.gz")
    model_list_gz <- file.path(raw_dir, "model_list_latest.csv.gz")
    cl_variants_gz <- file.path(raw_dir, "mutations_all_latest.csv.gz")

    # download & decompress Sanger files
    download.file(gene_annot_url, destfile = gene_annot_gz, mode = "wb")
    message("Gene identifiers downloaded (gzipped)")
    download.file(model_list_url, destfile = model_list_gz, mode = "wb")
    message("Model list downloaded (gzipped)")
    download.file(cl_variants_url, destfile = cl_variants_gz, mode = "wb")  
    message("Cell line variants downloaded (gzipped)")
  
    R.utils::gunzip(gene_annot_gz, destname = gene_annot_path, remove = TRUE, overwrite = TRUE)
    message("Gene identifiers decompressed")
    R.utils::gunzip(model_list_gz, destname = model_list_path, remove = TRUE, overwrite = TRUE)
    message("Model list decompressed")
    R.utils::gunzip(cl_variants_gz, destname = cl_variants_path, remove = TRUE, overwrite = TRUE) 
    message("Cell line variants decompressed")

    # -------- GDSC dose-response files from bulk download page -------- #
    bulk_url <- "https://cellmodelpassports.sanger.ac.uk/downloads"
    
    response <- GET(bulk_url)
    
    if (status_code(response) != 200) {
        stop("Failed to access bulk download page. Status: ", status_code(response))
    }
    
    page <- read_html(content(response, as = "text", encoding = "UTF-8"))
    
    # extract all links
    all_links <- page %>% 
        html_nodes("a") %>% 
        html_attr("href")
    
    # find GDSC fitted dose response links
    gdsc1_link <- all_links[grepl("GDSC1_fitted_dose_response.*\\.xlsx$", all_links)]
    gdsc2_link <- all_links[grepl("GDSC2_fitted_dose_response.*\\.xlsx$", all_links)]
    
    if (length(gdsc1_link) == 0 || length(gdsc2_link) == 0) {
        stop("Could not find GDSC fitted dose response files on the bulk download page")
    }
    
    # take the first match (the latest)
    gdsc1_url <- gdsc1_link[1]
    gdsc2_url <- gdsc2_link[1]
    
    # extract release version from URL
    release_version <- gsub(".*GDSC_release([0-9.]+)/.*", "\\1", gdsc1_url)
    
    message("Latest GDSC release detected: ", release_version)
    
    # define file paths
    gdsc1_xlsx_path <- file.path(raw_dir, paste0("GDSC1_fitted_dose_response_release", release_version, ".xlsx"))
    gdsc2_xlsx_path <- file.path(raw_dir, paste0("GDSC2_fitted_dose_response_release", release_version, ".xlsx"))
    gdsc1_csv_path <- file.path(raw_dir, "GDSC1_fitted_dose_response_latest.csv")
    gdsc2_csv_path <- file.path(raw_dir, "GDSC2_fitted_dose_response_latest.csv")
    
    # download GDSC files
    download.file(gdsc1_url, destfile = gdsc1_xlsx_path, mode = "wb")
    message("GDSC1 fitted dose-response downloaded")
    download.file(gdsc2_url, destfile = gdsc2_xlsx_path, mode = "wb")
    message("GDSC2 fitted dose-response downloaded")
    
    # convert to CSV
    gdsc1_data <- read_excel(gdsc1_xlsx_path)
    gdsc2_data <- read_excel(gdsc2_xlsx_path)
    
    write.csv(gdsc1_data, gdsc1_csv_path, row.names = FALSE)
    write.csv(gdsc2_data, gdsc2_csv_path, row.names = FALSE)
    message("GDSC files converted to CSV")

    # -------- IntOGen Drivers -------- #
    intogen_url <- "https://www.intogen.org/download"
    
    response_intogen <- GET(intogen_url)
    
    if (status_code(response_intogen) != 200) {
        stop("Failed to access IntOGen download page. Status: ", status_code(response_intogen))
    }
    
    page_intogen <- read_html(content(response_intogen, as = "text", encoding = "UTF-8"))
    
    # extract all links
    all_links_intogen <- page_intogen %>% 
        html_nodes("a") %>% 
        html_attr("href")
    
    # find IntOGen drivers links
    driver_links <- all_links_intogen[grepl("IntOGen-Drivers-[0-9]{8}\\.zip$", all_links_intogen)]
    
    if (length(driver_links) == 0) {
        stop("Could not find IntOGen-Drivers files on the download page")
    }
    
    # extract dates and find the latest
    dates <- gsub(".*IntOGen-Drivers-([0-9]{8})\\.zip", "\\1", driver_links)
    latest_idx <- which.max(as.numeric(dates))
    latest_file <- driver_links[latest_idx]
    intogen_date <- dates[latest_idx]
    
    message("Latest IntOGen-Drivers release detected: ", intogen_date)
    
    # construct full download URL
    intogen_download_url <- paste0("https://www.intogen.org/download", latest_file)
    
    # define file paths
    zip_path <- file.path(raw_dir, paste0("IntOGen-Drivers-", intogen_date, ".zip"))
    extract_dir <- file.path(raw_dir, "intogen_temp")
    final_tsv_path <- file.path(raw_dir, "Compendium_Cancer_Genes_latest.tsv")
    
    # download the zip file
    download.file(intogen_download_url, destfile = zip_path, mode = "wb")
    message("IntOGen-Drivers downloaded")
    
    # extract the zip file
    message("Extracting IntOGen zip file...")
    unzip(zip_path, exdir = extract_dir)
    
    # find the Compendium_Cancer_Genes.tsv file
    tsv_files <- list.files(extract_dir, pattern = "Compendium_Cancer_Genes\\.tsv$", 
                            recursive = TRUE, full.names = TRUE)
    
    if (length(tsv_files) == 0) {
        stop("Could not find Compendium_Cancer_Genes.tsv in the extracted files")
    }
    
    # copy to the raw directory with standardized name
    file.copy(tsv_files[1], final_tsv_path, overwrite = TRUE)
    message("Compendium_Cancer_Genes.tsv extracted and saved")
    
    # clean up temporary extraction directory
    unlink(extract_dir, recursive = TRUE)
  
    # -------- DepMap CRISPR file using API -------- # 
    api_url <- "https://depmap.org/portal/api/download/files"
    response <- GET(api_url)
    if (status_code(response) != 200) {stop("Failed to access DepMap API. Status: ", status_code(response))}
    files_df <- read_csv(content(response, as = "text", encoding = "UTF-8"), show_col_types = FALSE)
  
    # find the latest release
    latest_release <- files_df$release[1]  # files are sorted by release
    message("Latest DepMap release: ", latest_release)
  
    # filter for CRISPRGeneEffect.csv in the latest release
    crispr_file <- files_df[files_df$filename == "CRISPRGeneEffect.csv" & files_df$release == latest_release, ]
  
    if (nrow(crispr_file) == 0) {
        stop("Could not find CRISPRGeneEffect.csv in latest release")
    }
  
    crispr_url <- crispr_file$url[1]
  
    crispr_path <- file.path(raw_dir, "CRISPRGeneEffect_latest.csv")
  
    httr::GET(crispr_url, httr::write_disk(crispr_path, overwrite = TRUE), timeout(600))
    message("CRISPR data downloaded")
  
    return(list(
        gene_annot_path = gene_annot_path,
        model_list_path = model_list_path,
        crispr_path = crispr_path,
        cl_variants_path = cl_variants_path,
        gdsc1_path = gdsc1_csv_path,
        gdsc2_path = gdsc2_csv_path,
        intogen_path = final_tsv_path
    ))
}
