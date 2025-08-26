library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)
library(hdf5r)
library(dittoSeq)
library(SplineDV)
library(presto)
library(ggpubr)
library(tidyr)
library(stringr)
library(tibble)
library(cowplot)
library(openxlsx)
library(enrichR)

#base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
base_dir <- "/mnt/SCDC/Optimus/selim_working_dir/2023_nr4a1_colon/results"
setwd(base_dir)
data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))

# Ensure `cell_types` is properly extracted
cell_types <- unique(data$CellType)

# Define output directory for DV results
output_dir <- "DV_results_genes"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (cell in cell_types) {
  # Print the current cell type being processed
  cat("Processing cell type:", cell, "\n")
  
  # Subset the Seurat object for the current cell type
  cell_subset <- subset(data, CellType == cell)
  
  # Further subset for KO and WT conditions
  data_ko <- subset(cell_subset, Condition == 'Nr4a1 KO')
  data_wt <- subset(cell_subset, Condition == 'WT')
  
  # Extract counts for KO and WT
  data_ko <- GetAssayData(data_ko, layer = 'counts')  # slot = 'counts' instead of layers
  data_wt <- GetAssayData(data_wt, layer = 'counts')
  
  if (data_ko@Dim[2] > 100 && data_wt@Dim[2] > 100) {
    # Compute Differential Variability (DV)
    DV_res <- data.frame(splineDV(X = data_ko, Y = data_wt))
    }
  else{
    next
  }
  # Filter for significant results
  DV_res <- DV_res[DV_res$Pval < 0.05,]
  
  # Create a new workbook for DV results
  wb <- createWorkbook()
  
  # Add a worksheet to the workbook for the overall DV results
  addWorksheet(wb, "DV_Results")
  writeData(wb, sheet = "DV_Results", DV_res)
  
  # Ensure consistent file naming and save the workbook
  save_cell_name <- gsub(" ", "_", cell)  # Replace spaces with underscores *once*
  output_file <- file.path(output_dir, paste0(save_cell_name, "_DV_results.xlsx"))
  saveWorkbook(wb, file = output_file, overwrite = TRUE)
  
  cat("Saved DV results for:", cell, "\n")
}


# Define paths
input_dir <- "DV_results_genes"
output_dir <- "Enrichr_results/DV_Enrichr"

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set Enrichr database and parameters
setEnrichrSite('Enrichr')
dbs <- "GO_Biological_Process_2025"
pval_cutoff <- 0.05
top_n <- 500
# Function to perform Enrichr and process results with waiting time
run_enrichr <- function(genes, direction, dbs, pval_cutoff, wait_time = 5) {
  results_df <- data.frame()
  
  for (db in dbs) {
    enrichr_res <- enrichr(genes, db)[[1]]
    
    if (!is.null(enrichr_res) && nrow(enrichr_res) > 0) {
      enrichr_res <- enrichr_res %>% filter(Adjusted.P.value < pval_cutoff)
      
      if (nrow(enrichr_res) > 0) {
        results_df <- bind_rows(results_df, enrichr_res %>%
                                  transmute(database = db,
                                            Term, Overlap, `Odds.Ratio`,
                                            P.value, Adjusted.P.value,
                                            Genes)) # Use transmute
      }
    }
  }
  Sys.sleep(wait_time)  # Prevent API overload
  return(results_df)
}

# Get unique cell types (assuming 'data' is defined elsewhere)
cell_types <- unique(data$CellType)
for (cell in cell_types) {
  print(paste("Working with", cell))
  file_path <- file.path(input_dir, paste0(gsub(" ", "_", cell), "_DV_results.xlsx"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  
  wb <- loadWorkbook(file_path)
  dv_genes <- read.xlsx(wb, sheet = "DV_Results")

  # Subsample top_n genes for upregulated (Direction == 1)
  top_dv_genes <- dv_genes %>% head(top_n) %>% pull(genes)
  
  # Run Enrichr for upregulated and downregulated genes
  enrichr_results <- run_enrichr(top_dv_genes, paste0("Variable in KO ", cell), dbs, pval_cutoff)

  # Write results to Excel
  output_file <- file.path(output_dir, paste0(gsub("[^[:alnum:]_]", "", gsub(" ", "_", cell)), "_enrichr_DV_KO_var_results.xlsx"))
  wb_out <- createWorkbook()
  
  # Add upregulated results to the workbook if available
  if (nrow(enrichr_results) > 0) {
    addWorksheet(wb_out, "DV_Pathways")
    writeData(wb_out, "DV_Pathways", enrichr_results)
  } else {
    print(paste("No DV pathways for", cell))
  }
  
  # Save the workbook after both sheets have been added
  saveWorkbook(wb_out, output_file, overwrite = TRUE)
  print(paste("Enrichment results saved in", output_file))
}

