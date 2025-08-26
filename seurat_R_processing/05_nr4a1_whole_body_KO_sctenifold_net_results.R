library(scTenifoldNet)
library(scSGS)
library(openxlsx)
library(Seurat)
library(dplyr)
library(enrichR)
library(ggplot2)
library(pbapply)


base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
#base_dir <- "/mnt/SCDC/Optimus/selim_working_dir/2023_nr4a1_colon/results"
setwd(base_dir)
so <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))

um_tcells <- DimPlot(object = so, reduction = "umap", group.by = "CellType", 
                     repel = TRUE, alpha = 1, pt.size = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) +
  coord_cartesian(clip = "off")
um_tcells

#----------
# Load Data (All cell types)
#----------
setEnrichrSite('enrichr')
table(so$Condition, so$CellType)

# Drop Genes
rownames(so)[grepl('orf',rownames(so))]

soNoMito <- so[!startsWith(rownames(so), 'MT-'), ]
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'Mt-'), ]
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'mt-'), ]
# Remove unnamed genes (ENSG0)
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'ENSG0'), ]
# Remove ORF genes (case-insensitive)
soNoMito <- soNoMito[!grepl('orf', rownames(soNoMito), ignore.case = TRUE), ]
# Remove antisense genes (case-insensitive)
soNoMito <- soNoMito[!grepl('-AS', rownames(soNoMito), ignore.case = TRUE), ]
soNoMito
gc()
rm(so)
#----------
# scTenifoldNet Analysis in loop
#----------
nCores <- 10

for (CT in unique(soNoMito$CellType)){
  so1 <- subset(soNoMito, subset = CellType == CT)
  print(CT)
  
  soKO <- GetAssayData(subset(so1, subset = Condition == 'Nr4a1 KO'), layer = 'count')
  soWT <- GetAssayData(subset(so1, subset = Condition == 'WT'), layer = 'count')
  rm(so1)
  gc()
  
  # Filter by Dropout. Preserve genes that express in at least 10% cells
  soKO <- soKO[Matrix::rowSums(soKO > 0) / ncol(soKO) > 0.1, ]
  soWT <- soWT[Matrix::rowSums(soWT > 0) / ncol(soWT) > 0.1, ]
  
  # Create Excel Workbook
  wb <- createWorkbook()
  
  # scTenifoldNet
  tryCatch({
    # Comparing gene ids.
    sharedGenes <- intersect(rownames(soKO), rownames(soWT))
    
    # Filtering out non-shared genes
    numSharedGenes <- length(sharedGenes)
    print(paste("Number of shared genes:", numSharedGenes))      
    
    X <- soKO[sharedGenes,]
    Y <- soWT[sharedGenes,]
    
    # Construction of gene-regulatory networks based on principal component regression (pcNet).
    set.seed(1)
    old <- Sys.time()
    xNet <- pcNet(X, nCores = nCores)
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    print('Net1 constructed...')
    set.seed(1)
    
    old <- Sys.time()
    yNet <- pcNet(Y, nCores = nCores)
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    print('Net2 constructed...')
    
    set.seed(1)
    # Manifold alignment
    mA <- manifoldAlignment(xNet, yNet, d = 30, nCores = nCores)
    rownames(mA) <- c(paste0('X_', sharedGenes),paste0('Y_', sharedGenes))
    
    # Differential regulation testing
    dR <- dRegulation(manifoldOutput = mA)
    
    resDF <- dR
    
    addWorksheet(wb, sheetName = 'scTenifoldNet', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'scTenifoldNet', x = resDF)
    
    # Enrichr
    #dbs <- c('GO_Biological_Process_2023','GO_Molecular_Function_2023')
    dbs <- c('GO_Biological_Process_2023')
    glist <- resDF %>% arrange(p.value) %>% head(250) %>% filter(p.value < 0.05) %>%
      pull(gene)
    addWorksheet(wb, sheetName = 'DRGenesForEnrichr', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'DRGenesForEnrichr',
                   x = resDF[resDF$gene %in% glist,])
    
    enrichrRes <- enrichr(glist, dbs)
    
    addWorksheet(wb, sheetName = 'enrichr_GOBP', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'enrichr_GOBP',
                   x = enrichrRes[[1]] %>% filter(Adjusted.P.value < 0.05))
    addWorksheet(wb, sheetName = 'enrichr_GOMF', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'enrichr_GOMF',
                   x = enrichrRes[[2]] %>% filter(Adjusted.P.value < 0.05))
    
    # Save workbook
    saveWorkbook(wb, paste0('scTenifoldNet_results/',CT,'_sctenifold_results.xlsx'),
                 overwrite = TRUE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



#----------
# scTenifoldNet Analysis in loop for SubPopulations
#----------
nCores <- 10
for (CT in unique(soNoMito$SubCellType)){
  so1 <- subset(soNoMito, subset = SubCellType == CT)
  print(CT)
  
  soKO <- GetAssayData(subset(so1, subset = Condition == 'Nr4a1 KO'), assay = 'RNA', layer = 'count')
  soWT <- GetAssayData(subset(so1, subset = Condition == 'WT'), assay = 'RNA', layer = 'count')

  rm(so1)
  gc()
  
  # Filter by Dropout. Preserve genes that express in at least 10% cells
  soKO <- soKO[Matrix::rowSums(soKO > 0) / ncol(soKO) > 0.1, ]
  soWT <- soWT[Matrix::rowSums(soWT > 0) / ncol(soWT) > 0.1, ]
  
  # Create Excel Workbook
  wb <- createWorkbook()
  
  # scTenifoldNet
  tryCatch({
    # Comparing gene ids.
    sharedGenes <- intersect(rownames(soKO), rownames(soWT))
    
    # Filtering out non-shared genes
    numSharedGenes <- length(sharedGenes)
    print(paste("Number of shared genes:", numSharedGenes))      
    
    X <- soKO[sharedGenes,]
    Y <- soWT[sharedGenes,]
    
    # Construction of gene-regulatory networks based on principal component regression (pcNet).
    set.seed(1)
    old <- Sys.time()
    xNet <- pcNet(X, nCores = nCores)
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    print('Net1 constructed...')
    set.seed(1)
    
    old <- Sys.time()
    yNet <- pcNet(Y, nCores = nCores)
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    print('Net2 constructed...')
    
    set.seed(1)
    # Manifold alignment
    mA <- manifoldAlignment(xNet, yNet, d = 30, nCores = nCores)
    rownames(mA) <- c(paste0('X_', sharedGenes),paste0('Y_', sharedGenes))
    
    # Differential regulation testing
    dR <- dRegulation(manifoldOutput = mA)
    
    resDF <- dR
    
    addWorksheet(wb, sheetName = 'scTenifoldNet', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'scTenifoldNet', x = resDF)
    
    # Enrichr
    #dbs <- c('GO_Biological_Process_2023','GO_Molecular_Function_2023')
    dbs <- c('GO_Biological_Process_2025')
    
    glist <- resDF %>% arrange(p.value) %>% head(250) %>% filter(p.value < 0.05) %>%
      pull(gene)
    addWorksheet(wb, sheetName = 'DRGenesForEnrichr', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'DRGenesForEnrichr',
                   x = resDF[resDF$gene %in% glist,])
    
    enrichrRes <- enrichr(glist, dbs)
    
    addWorksheet(wb, sheetName = 'enrichr_GOBP', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'enrichr_GOBP',
                   x = enrichrRes[[1]] %>% filter(Adjusted.P.value < 0.05))
    # addWorksheet(wb, sheetName = 'enrichr_GOMF', gridLines = FALSE)
    # writeDataTable(wb = wb, sheet = 'enrichr_GOMF',
    #                x = enrichrRes[[2]] %>% filter(Adjusted.P.value < 0.05))
    
    # Save workbook
    saveWorkbook(wb, paste0('scTenifoldNet_results/subcolonocytes/',CT,'_sctenifold_results.xlsx'),
                 overwrite = TRUE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}




## Stem cells only cell type analysis wit Wnt pathway genes

base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
#base_dir <- "/mnt/SCDC/Optimus/selim_working_dir/2023_nr4a1_colon/results"
setwd(base_dir)
so <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))

wnt_pathway_genes <- c('Ncoa3', 'Csnk1a1', 'Ryk', 'Jup', 'Wnt5a', 'Cav1', 'Wnt5b', 'Axin2', 
                       'Wnt9a', 'Wnt9b', 'Wnt16', 'Sfrp4', 'Sfrp5', 'Bcl9', 'Sfrp2', 'Cpe', 
                       'Cd24', 'Lef1', 'Ddx3x', 'Gsk3b', 'Usp34', 'Lrrk2', 'Lrp6', 'Wnt8b', 
                       'Lrp5', 'Fzd10', 'Wnt8a', 'Celsr2', 'Sfrp1', 'Wnt6', 'Frat1', 'Frat2', 
                       'Frzb', 'Wnt11', 'Tmem67', 'Dvl1', 'Dvl2', 'Dvl3', 'Gpc4', 'Myoc', 
                       'Wnt1', 'Wnt2', 'Wnt3', 'Wnt4', 'Fzd2', 'Siah1', 'Fzd1', 'Fzd4', 
                       'Tcf7l2', 'Wnt10b', 'Tcf7l1', 'Wnt10a', 'Fzd6', 'Wnt3a', 'Fzd5', 
                       'Fzd8', 'Fzd7', 'Siah2', 'Wnt7a', 'Fzd9', 'Wnt7b', 'Csnk1e', 
                       'Ccdc88c', 'Fermt2', 'Tcf7', 'Prkd2', 'Snai1', 'Strn', 
                       'Rab5a', 'Nxn', 'Pkd1', 'Wnt2b', 'Arl6', 'Rarg', 'Ctnnb1', 
                       'Dixdc1', 'Klhl12', 'Gata3', 'Porcn', 'Plcg2', 'Vps35', 
                       'Bcl9l', 'Nr4a2', 'Nr4a1')


soNoMito <- subset(so, subset = SubCellType == "Stem")
table(soNoMito$Condition, soNoMito$SubCellType)

#----------
# Load Data (Sub cell types in colonocytes)
#----------
setEnrichrSite('enrichr')

# Drop Genes
rownames(soNoMito)[grepl('orf',rownames(soNoMito))]

soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'MT-'), ]
rm(so)
gc()
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'Mt-'), ]
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'mt-'), ]
gc()
# Remove unnamed genes (ENSG0)
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'ENSG0'), ]
# Remove ORF genes (case-insensitive)
soNoMito <- soNoMito[!grepl('orf', rownames(soNoMito), ignore.case = TRUE), ]
gc()
# Remove antisense genes (case-insensitive)
soNoMito <- soNoMito[!grepl('-AS', rownames(soNoMito), ignore.case = TRUE), ]
gc()

soNoMito <- soNoMito[wnt_pathway_genes, ]
DefaultAssay(soNoMito)<- "RNA"
#----------
# scTenifoldNet Analysis in loop
#----------
nCores <- 2
for (CT in unique(soNoMito$SubCellType)){
  so1 <- subset(soNoMito, subset = SubCellType == CT)
  print(CT)
  
  soKO <- GetAssayData(subset(so1, subset = Condition == 'Nr4a1 KO'), assay = 'RNA', layer = 'count')
  soWT <- GetAssayData(subset(so1, subset = Condition == 'WT'), assay = 'RNA', layer = 'count')
  
  rm(so1)
  gc()
  
  # Filter by Dropout. Preserve genes that express in at least 10% cells
  soKO <- soKO[Matrix::rowSums(soKO > 0) / ncol(soKO) > 0.1, ]
  soWT <- soWT[Matrix::rowSums(soWT > 0) / ncol(soWT) > 0.1, ]
  
  # Create Excel Workbook
  wb <- createWorkbook()
  
  # scTenifoldNet
  tryCatch({
    # Comparing gene ids.
    sharedGenes <- intersect(rownames(soKO), rownames(soWT))
    
    # Filtering out non-shared genes
    numSharedGenes <- length(sharedGenes)
    print(paste("Number of shared genes:", numSharedGenes))      
    
    X <- soKO[sharedGenes,]
    Y <- soWT[sharedGenes,]
    
    # Construction of gene-regulatory networks based on principal component regression (pcNet).
    set.seed(1)
    old <- Sys.time()
    xNet <- pcNet(X, nCores = nCores)
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    print('Net1 constructed...')
    set.seed(1)
    
    old <- Sys.time()
    yNet <- pcNet(Y, nCores = nCores)
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    print('Net2 constructed...')
    
    set.seed(1)
    # Manifold alignment
    mA <- manifoldAlignment(xNet, yNet, d = 30, nCores = nCores)
    rownames(mA) <- c(paste0('X_', sharedGenes),paste0('Y_', sharedGenes))
    
    # Differential regulation testing
    dR <- dRegulation(manifoldOutput = mA)
    
    resDF <- dR
    
    addWorksheet(wb, sheetName = 'scTenifoldNet', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'scTenifoldNet', x = resDF)
    
    # Enrichr
    #dbs <- c('GO_Biological_Process_2023','GO_Molecular_Function_2023')
    dbs <- c('GO_Biological_Process_2023')
    
    glist <- resDF %>% arrange(p.value) %>% head(250) %>% filter(p.value < 0.05) %>%
      pull(gene)
    addWorksheet(wb, sheetName = 'DRGenesForEnrichr', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'DRGenesForEnrichr',
                   x = resDF[resDF$gene %in% glist,])
    
    enrichrRes <- enrichr(glist, dbs)
    
    addWorksheet(wb, sheetName = 'enrichr_GOBP', gridLines = FALSE)
    writeDataTable(wb = wb, sheet = 'enrichr_GOBP',
                   x = enrichrRes[[1]] %>% filter(Adjusted.P.value < 0.05))
    # addWorksheet(wb, sheetName = 'enrichr_GOMF', gridLines = FALSE)
    # writeDataTable(wb = wb, sheet = 'enrichr_GOMF',
    #                x = enrichrRes[[2]] %>% filter(Adjusted.P.value < 0.05))
    
    # Save workbook
    saveWorkbook(wb, paste0('scTenifoldNet_results/subcolonocytes/',CT,'_sctenifold_results_wnt.xlsx'),
                 overwrite = TRUE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}





## Colonocytes sub cell type analysis

base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
#base_dir <- "/mnt/SCDC/Optimus/selim_working_dir/2023_nr4a1_colon/results"
setwd(base_dir)
so <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))

#----------
# Load Data (Sub cell types in colonocytes)
#----------
setEnrichrSite('enrichr')
table(so$Condition, so$CellType)

# Drop Genes
rownames(so)[grepl('orf',rownames(so))]

soNoMito <- so[!startsWith(rownames(so), 'MT-'), ]
rm(so)
gc()
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'Mt-'), ]
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'mt-'), ]
gc()
# Remove unnamed genes (ENSG0)
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'ENSG0'), ]
# Remove ORF genes (case-insensitive)
soNoMito <- soNoMito[!grepl('orf', rownames(soNoMito), ignore.case = TRUE), ]
gc()
# Remove antisense genes (case-insensitive)
soNoMito <- soNoMito[!grepl('-AS', rownames(soNoMito), ignore.case = TRUE), ]
soNoMito <- subset(soNoMito, subset = CellType == "Colonocytes")
table(soNoMito$Condition, soNoMito$SubCellType)
gc()


# --- Re-run Enrichr GOBP Analysis with 2025 database ---
# The scTenifoldNet computation (pcNet, manifoldAlignment, dRegulation) is entirely skipped.
# We will load the previously computed 'resDF' from the Excel files.

# Loop through each unique SubCellType identified in the data
for (CT in unique(soNoMito$SubCellType)){
  print(paste("Attempting to update GOBP analysis for SubCellType:", CT))
  
  # Construct the full path to the previously generated Excel file for the current SubCellType
  excel_path <- file.path(base_dir, 'scTenifoldNet_results', 'subcolonocytes', paste0(CT,'_sctenifold_results.xlsx'))
  
  # Check if the Excel file exists before attempting to load it.
  # If the file doesn't exist, we cannot perform the Enrichr analysis.
  if (!file.exists(excel_path)) {
    warning(paste("Skipping", CT, ": Previous scTenifoldNet results file not found at", excel_path, ". Please ensure it exists."))
    next # Skip to the next SubCellType in the loop
  }
  
  # Use a tryCatch block to gracefully handle any errors during processing a specific file.
  tryCatch({
    # Load the existing workbook from the specified path
    wb <- loadWorkbook(excel_path)
    print(paste("Workbook successfully loaded for", CT))
    
    # Read the 'scTenifoldNet' sheet to retrieve the differential regulation results (resDF).
    # This sheet should contain the 'gene' and 'p.value' columns needed for Enrichr.
    if (!("scTenifoldNet" %in% names(wb))) {
      warning(paste("Skipping", CT, ": Sheet 'scTenifoldNet' not found in the workbook at", excel_path, ". It's needed to extract genes for Enrichr."))
      next
    }
    resDF <- read.xlsx(wb, sheet = 'scTenifoldNet')
    print(paste("Differential regulation data ('scTenifoldNet' sheet) read for", CT))
    
    # Extract the gene list ('glist') for Enrichr analysis.
    # This logic remains the same: top 250 genes with p.value < 0.05.
    glist <- resDF %>% 
      arrange(p.value) %>% 
      head(250) %>% 
      filter(p.value < 0.05) %>%
      pull(gene)
    
    # If no significant genes are found, skip Enrichr for this SubCellType.
    if (length(glist) == 0) {
      message(paste("No significant genes (p.value < 0.05 among top 250) found for Enrichr analysis in", CT, ". Skipping Enrichr for this SubCellType."))
      next
    }
    print(paste("Gene list prepared for", CT, "with", length(glist), "genes."))
    
    # Specify the Enrichr database to use: 'GO_Biological_Process_2025'.
    dbs <- c('GO_Biological_Process_2025')
    
    # Perform the Enrichr analysis with the extracted gene list and specified database.
    enrichrRes <- enrichr(glist, dbs)
    print(paste("Enrichr analysis completed for", CT, "using GO_Biological_Process_2025."))
    
    # Before adding the new sheet, remove the old 'enrichr_GOBP' sheet if it exists.
    # This ensures that the sheet is cleanly overwritten with the new 2025 data.
    if ("enrichr_GOBP" %in% names(wb)) {
      removeWorksheet(wb, sheet = 'enrichr_GOBP')
      print(paste("Existing 'enrichr_GOBP' sheet removed for", CT))
    }
    
    # Add a new worksheet specifically for the updated GOBP results.
    addWorksheet(wb, sheetName = 'enrichr_GOBP', gridLines = FALSE)
    
    # Filter the Enrichr results for terms with an Adjusted.P.value < 0.05.
    gobp_results <- enrichrRes[[1]] %>% filter(Adjusted.P.value < 0.05)
    
    # Write the filtered GOBP results to the new sheet.
    if (nrow(gobp_results) == 0) {
      message(paste("No significant GOBP terms (Adjusted.P.value < 0.05) found for", CT, ". An empty sheet will be added to indicate processing."))
      # Write an empty data frame if no significant terms, to show processing occurred.
      writeDataTable(wb = wb, sheet = 'enrichr_GOBP', x = data.frame()) 
    } else {
      writeDataTable(wb = wb, sheet = 'enrichr_GOBP', x = gobp_results)
      print(paste("Updated 'enrichr_GOBP' sheet added for", CT, "with", nrow(gobp_results), "significant terms."))
    }
    
    # Save the workbook, overwriting the original file. This updates the Excel file in place.
    saveWorkbook(wb, excel_path, overwrite = TRUE)
    print(paste("Workbook successfully saved (updated) for", CT, "at", excel_path))
    
  }, error=function(e){
    # Catch and print any errors that occur for a specific SubCellType,
    # allowing the loop to continue processing other files.
    cat("ERROR processing SubCellType", CT, ":", conditionMessage(e), "\n")
  })
}

print("GOBP analysis update process complete for all colonocyte sub-cell types.")



base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
#base_dir <- "/mnt/SCDC/Optimus/selim_working_dir/2023_nr4a1_colon/results"
setwd(base_dir)
so <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))
DefaultAssay(so) <- "RNA"
#----------
# Load Data (Sub cell types in colonocytes) analysis scSGS
#----------
table(so$Condition, so$CellType)

so1 <- subset(so, subset = SubCellType == "Stem")
# Drop Genes
rownames(so1)[grepl('orf', rownames(so1))]

soNoMito <- so1[!startsWith(rownames(so1), 'MT-'), ]
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'Mt-'), ]
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'mt-'), ]
# Remove unnamed genes (ENSG0)
soNoMito <- soNoMito[!startsWith(rownames(soNoMito), 'ENSG0'), ]
# Remove ORF genes (case-insensitive)
soNoMito <- soNoMito[!grepl('orf', rownames(soNoMito), ignore.case = TRUE), ]
# Remove antisense genes (case-insensitive)
soNoMito <- soNoMito[!grepl('-AS', rownames(soNoMito), ignore.case = TRUE), ]

table(soNoMito$Condition, soNoMito$SubCellType)
gc()
rm(so, so1)

#soWT <- GetAssayData(subset(soNoMito, subset = Condition == 'WT'), assay = 'RNA', layer = 'count')
soKO <- subset(soNoMito, subset = Condition == 'Nr4a1 KO')
soWT <- subset(soNoMito, subset = Condition == 'WT')
gc()

HVG_res = HVG_splinefit(GetAssayData(soWT, assay = 'RNA', layer = 'counts'), diptest=TRUE,
                         dropout.filter=TRUE, nHVGs=5000)
# HVG results
HVG_only <- HVG_res[HVG_res$HVG == TRUE, ]
# Get the row number for "Nr4a1"
row_number_Nr4a1 <- which(rownames(HVG_res) == "Nr4a1")
print(row_number_Nr4a1)
# Get the row number for "Nr4a1"
row_number_Nr4a2 <- which(rownames(HVG_res) == "Nr4a2")
print(row_number_Nr4a2)


SGS_res = scSGS(GetAssayData(soWT, layer='counts'), 'Nr4a1', nHVG=5000,
                rm.mt=T, rm.rp=T, calcHVG=TRUE, verbose=TRUE,)

# scSGS DE results
head(SGS_res$DE)
View(SGS_res$HVG_df)

clipr::write_clip(rownames(SGS_res$DE %>% filter(abs(avg_log2FC) > 1)))

wb <- createWorkbook()
addWorksheet(wb, sheetName = 'scSGS_HVGs', gridLines = FALSE)
writeDataTable(wb = wb, sheet = 'scSGS_HVGs', x = SGS_res$HVG_df)
addWorksheet(wb, sheetName = 'scSGS_DE', gridLines = FALSE)
writeDataTable(wb = wb, sheet = 'scSGS_DE', x = SGS_res$DE)
saveWorkbook(wb, paste0('stem_cells_scSGS.xlsx'),
             overwrite = TRUE)

