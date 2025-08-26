# Load necessary libraries
library(dplyr)
library(Seurat) # Assuming 'data' object might be a Seurat object or related structure
library(ggplot2)
library(Matrix)
library(dittoSeq)
library(presto)
library(ggpubr)
library(tidyr)
library(stringr)
library(tibble)
library(cowplot)
library(openxlsx)
library(enrichR)

# --- Configuration ---
# !!! IMPORTANT: Please verify this base directory path is correct for your system !!!
base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
setwd(base_dir)

# Load main data object (assuming it contains a 'CellType' column)
data_file_path <- file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds")
if (!file.exists(data_file_path)) {
  stop(paste("Main data file not found:", data_file_path))
}
data <- readRDS(data_file_path)

input_dv_dir <- "DV_results_genes"
input_de_dir <- "DE_results_genes"
output_dir <- "Enrichr_results/DV_plus_DE_Enrichr"

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Enrichr configuration
setEnrichrSite('Enrichr') # Ensure this is the correct site ('Enrichr', 'FlyEnrichr', 'HumanEnrichr', etc.)
dbs <- "GO_Biological_Process_2025" # Can be a single string or a vector like c("GO_Biological_Process_2023", "KEGG_2021_Human")
pval_cutoff <- 0.05
api_wait_time <- 2 # Seconds to wait between Enrichr API calls per database

# --- Function to perform Enrichr and process results ---
run_enrichr <- function(genes, direction_description, dbs_list, pval_cutoff_value, wait_time = 2) {
  results_df <- data.frame() # Initialize an empty data frame
  
  if (!is.character(genes) || length(genes) == 0) {
    message(paste("No genes provided or invalid format for", direction_description, "- skipping Enrichr for this set."))
    return(results_df)
  }
  
  for (db_item in dbs_list) {
    message(paste("Querying Enrichr with database:", db_item, "for", direction_description, "(", length(genes), "genes )"))
    enrichr_output <- NULL 
    
    tryCatch({
      enrichr_output <- enrichR::enrichr(genes, db_item) 
    }, error = function(e) {
      message(paste("Error during Enrichr API call for database '", db_item, "' in '", direction_description, "': ", e$message, sep=""))
    })
    
    if (!is.null(enrichr_output) && length(enrichr_output) > 0 && inherits(enrichr_output[[1]], "data.frame")) {
      enrichr_res_raw <- enrichr_output[[1]] 
      
      if (nrow(enrichr_res_raw) > 0 && "Adjusted.P.value" %in% colnames(enrichr_res_raw)) {
        enrichr_res_filtered <- enrichr_res_raw %>% 
          filter(Adjusted.P.value < pval_cutoff_value)
        
        if (nrow(enrichr_res_filtered) > 0) {
          required_cols <- c("Term", "Overlap", "Odds.Ratio", "P.value", "Adjusted.P.value", "Genes")
          if (all(required_cols %in% names(enrichr_res_filtered))) {
            results_df <- bind_rows(results_df, enrichr_res_filtered %>%
                                      transmute(database = db_item,
                                                Term, Overlap, Odds.Ratio, 
                                                P.value, Adjusted.P.value,
                                                Genes))
          } else {
            message(paste("One or more required columns missing in filtered results for database '", db_item, "' for '", direction_description, "'. Missing: ", paste(setdiff(required_cols, names(enrichr_res_filtered)), collapse=", "), sep=""))
          }
        } else {
          message(paste("No terms passed the p-value cutoff for database '", db_item, "' in '", direction_description, "'.", sep=""))
        }
      } else {
        message(paste("No results or 'Adjusted.P.value' column missing after querying database '", db_item, "' for '", direction_description, "'. Raw rows: ", nrow(enrichr_res_raw), sep=""))
      }
    } else {
      message(paste("No valid data frame returned from Enrichr for database '", db_item, "' for '", direction_description, "'.", sep=""))
    }
    Sys.sleep(wait_time) 
  }
  return(results_df)
}




# --- Main processing loop for each cell type ---
cell_types_vector <- if (inherits(data, "Seurat")) unique(data@meta.data$CellType) else unique(data$CellType)

for (cell in cell_types_vector) {
  cell_label <- gsub(" ", "_", cell) 
  cell_label <- gsub("[^[:alnum:]_]", "", cell_label) # Further sanitize for file names
  print(paste("--- Working with cell type:", cell, "(Label:", cell_label, ") ---"))
  
  file_path_dv <- file.path(input_dv_dir, paste0(cell_label, "_DV_results.xlsx"))
  file_path_de <- file.path(input_de_dir, paste0(cell_label, "_markers_Wilcoxon.xlsx"))
  
  if (!file.exists(file_path_dv)) {
    warning(paste("DV results file not found for", cell, ":", file_path_dv, "- Skipping this cell type."))
    next
  }
  if (!file.exists(file_path_de)) {
    warning(paste("DE results file not found for", cell, ":", file_path_de, "- Skipping this cell type."))
    next
  }
  
  dv_genes_df <- data.frame()
  tryCatch({
    wb_dv <- loadWorkbook(file_path_dv)
    if ("DV_Results" %in% names(wb_dv)) {
      dv_genes_df <- read.xlsx(wb_dv, sheet = "DV_Results")
    } else {
      warning(paste("'DV_Results' sheet not found in", file_path_dv))
    }
  }, error = function(e) {
    warning(paste("Error loading DV workbook for", cell, ":", e$message))
  })
  if (nrow(dv_genes_df) == 0 || !"genes" %in% names(dv_genes_df)) {
    warning(paste("No data or 'genes' column in DV_Results for", cell, "- Skipping enrichment for this cell type."))
    next
  }
  top_dv_genes <- dv_genes_df %>% pull(genes)
  
  de_genes_df <- data.frame()
  tryCatch({
    wb_de <- loadWorkbook(file_path_de)
    if ("All_Pval_0.05" %in% names(wb_de)) {
      de_genes_df <- read.xlsx(wb_de, sheet = "All_Pval_0.05")
    } else {
      warning(paste("'All_Pval_0.05' sheet not found in", file_path_de))
    }
  }, error = function(e) {
    warning(paste("Error loading DE workbook for", cell, ":", e$message))
  })
  if (nrow(de_genes_df) == 0 || !"Gene" %in% names(de_genes_df)) {
    warning(paste("No data or 'Gene' column in DE results for", cell, "- Skipping enrichment for this cell type."))
    next
  }
  
  # Overlap of DV and DE genes
  select_genes <- intersect(de_genes_df$Gene, top_dv_genes)
  if (length(select_genes) == 0) {
    print(paste("No common genes found between DV and DE for", cell, "- Skipping enrichment for this cell type."))
    next
  }
  
  de_genes_filtered <- de_genes_df %>% filter(Gene %in% select_genes)
  
  top_upregulated <- character(0) 
  top_downregulated <- character(0) 
  
  if (!"avg_log2FC" %in% names(de_genes_filtered)) {
    warning(paste("'avg_log2FC' column not found in filtered DE results for", cell, "- Cannot separate up/down regulated genes."))
  } else {
    upregulated_df <- de_genes_filtered %>% filter(avg_log2FC > 0)
    downregulated_df <- de_genes_filtered %>% filter(avg_log2FC < 0)
    
    if (nrow(upregulated_df) > 0) top_upregulated <- upregulated_df %>% pull(Gene)
    if (nrow(downregulated_df) > 0) top_downregulated <- downregulated_df %>% pull(Gene)
  }
  
  enrichr_results_up <- data.frame()
  enrichr_results_dn <- data.frame()
  
  description_base <- paste0("DV_plus_DE_in_KO_", cell_label) 
  
  if (length(top_upregulated) > 0) {
    enrichr_results_up <- run_enrichr(genes = top_upregulated, 
                                      direction_description = paste(description_base, "UP"), 
                                      dbs_list = dbs, 
                                      pval_cutoff_value = pval_cutoff,
                                      wait_time = api_wait_time)
  } else {
    print(paste("Skipping enrichment for UPregulated genes in", cell, "- No genes in the list."))
  }
  
  if (length(top_downregulated) > 0) {
    enrichr_results_dn <- run_enrichr(genes = top_downregulated, 
                                      direction_description = paste(description_base, "DN"), 
                                      dbs_list = dbs, 
                                      pval_cutoff_value = pval_cutoff,
                                      wait_time = api_wait_time)
  } else {
    print(paste("Skipping enrichment for DOWNregulated genes in", cell, "- No genes in the list."))
  }
  
  # --- Write results to Excel ---
  # Create workbook if there are input gene lists OR enrichment results
  should_create_workbook <- (length(top_upregulated) > 0 || length(top_downregulated) > 0 || 
                               nrow(enrichr_results_up) > 0 || nrow(enrichr_results_dn) > 0)
  
  if (should_create_workbook) {
    output_file_name <- paste0(cell_label, "_enrichr_DV_DE_KO_results.xlsx")
    output_file_path <- file.path(output_dir, output_file_name)
    wb_out <- createWorkbook()
    
    # Add Input Gene Lists sheet
    if (length(top_upregulated) > 0 || length(top_downregulated) > 0) {
      len_up <- length(top_upregulated)
      len_down <- length(top_downregulated)
      max_len_lists <- max(len_up, len_down)
      
      df_gene_lists <- data.frame(
        Upregulated_Genes = if(len_up > 0) c(top_upregulated, rep(NA_character_, max_len_lists - len_up)) else rep(NA_character_, max_len_lists),
        Downregulated_Genes = if(len_down > 0) c(top_downregulated, rep(NA_character_, max_len_lists - len_down)) else rep(NA_character_, max_len_lists)
      )
      # Ensure data frame has columns even if max_len_lists is 0 (both lists empty)
      if (max_len_lists == 0) {
        df_gene_lists <- data.frame(Upregulated_Genes = character(0), Downregulated_Genes = character(0))
      }
      
      addWorksheet(wb_out, "Input_Gene_Lists")
      writeData(wb_out, "Input_Gene_Lists", df_gene_lists)
    }
    
    # Add UP pathways sheet
    if (nrow(enrichr_results_up) > 0) {
      addWorksheet(wb_out, "DV_DE_UP_Pathways")
      writeData(wb_out, "DV_DE_UP_Pathways", enrichr_results_up)
    } else {
      addWorksheet(wb_out, "DV_DE_UP_Pathways_EMPTY")
      writeData(wb_out, "DV_DE_UP_Pathways_EMPTY", data.frame(Message = paste("No significant UP pathways found or no UP genes for", cell)))
      print(paste("No significant DV_DE_UP pathways to write for", cell))
    }
    
    # Add DN pathways sheet
    if (nrow(enrichr_results_dn) > 0) {
      addWorksheet(wb_out, "DV_DE_DN_Pathways")
      writeData(wb_out, "DV_DE_DN_Pathways", enrichr_results_dn)
    } else {
      addWorksheet(wb_out, "DV_DE_DN_Pathways_EMPTY")
      writeData(wb_out, "DV_DE_DN_Pathways_EMPTY", data.frame(Message = paste("No significant DN pathways found or no DN genes for", cell)))
      print(paste("No significant DV_DE_DN pathways to write for", cell))
    }
    
    tryCatch({
      saveWorkbook(wb_out, output_file_path, overwrite = TRUE)
      print(paste("Results and input gene lists saved in", output_file_path))
    }, error = function(e){
      warning(paste("Error saving workbook for", cell, ":", e$message))
    })
    
  } else {
    print(paste("No input genes to list and no significant enrichment results for UP or DN pathways in", cell, "- Excel file not created."))
  }
}

print("--- Script finished ---")

