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
library(patchwork)
library(openxlsx)
library(AUCell)
library(ggsignif)

base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
setwd(base_dir)
data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))

process_and_extract_cell_type <- function(data, cell_type_name, assay_name = "RNA", num_hvg = 2000, dims_pca = 50, resolution = 1.0) {
  DefaultAssay(data) <- assay_name
  data_cell_type <- subset(data, subset = CellType == cell_type_name)
  data_cell_type <- FindVariableFeatures(data_cell_type, selection.method = "vst", nfeatures = num_hvg)
  data_cell_type <- ScaleData(data_cell_type)
  data_cell_type <- RunPCA(data_cell_type)
  data_cell_type <- RunUMAP(data_cell_type, dims = 1:dims_pca, n.epochs = 500)
  data_cell_type <- FindNeighbors(data_cell_type, dims = 1:dims_pca)
  data_cell_type <- FindClusters(data_cell_type, resolution = resolution)
  return(data_cell_type)
}

DefaultAssay(data) <- "RNA"
# Substract colonocytes cells
data_colo <- process_and_extract_cell_type(data, "Colonocytes", resolution = 1.0)
# Re-factor for plotting 
# data_colo$SubCellType = factor(data_colo$SubCellType,
#                                  levels = c('Stem','TA','DCS','Imm. Goblet',
#                                             'Goblet','Aqp8+ Abs.',
#                                             'Car1+ Abs.',
#                                             'Tuft','EC',
#                                             'Cck+ EECs'))
data_colo$SubCellType = factor(data_colo$SubCellType,
                               levels = c('Stem','TA','DCS','Goblet',
                                          'Abs. Col', 'EECs', 'Tuft'))

# Save t cells rds
saveRDS(data_colo, file.path(base_dir, "colonocytes_only.rds"))

# Load this for faster processing
data_colo <- readRDS(file.path(base_dir, "colonocytes_only.rds"))


DimPlot(data_colo, reduction = "umap", group.by = 'SubCellType',
        label = T, label.box = T)+ NoLegend()


## Canonical marker dotplots ####
p1 = DotPlot(data_colo, features = c('Lgr5','Smoc2',
                                     'Mki67','Kcnq1','Cftr',
                                     'Reg4','Ccl6',
                                     'Muc2','Spdef','Tff3',
                                     'Aqp8','Ceacam20',
                                     'Tph1','Chga',
                                     'Trpm5','Dclk1'),dot.min = 0.15,
             dot.scale = 8, scale = T,
             group.by = 'SubCellType') +
  scale_color_gradientn(colours = c('dodgerblue3','grey95','firebrick3')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = 10, color = "black"),  # Keep axis titles
    legend.text = element_text(size = 16),  # Increase legend label size
    legend.title = element_text(size = 14)  # Increase legend title size
  ) 

p1
#png('plots/Colonocytes_Dot.png', height = 4, width = 10, res = 600, units = 'in')
ggsave("plots/Colonocytes_Dot.svg", plot = p1, height = 4, width = 10, device = "svg")
p1
dev.off()


# UMAP for SubCellType
um_colo <- DimPlot(object = data_colo, reduction = "umap", group.by = "SubCellType", 
                   label = TRUE, repel = TRUE, label.size = 7, 
                   label.box = TRUE, alpha = 1, raster = TRUE) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 10, color = "black"),  # Keep axis titles
    legend.text = element_text(size = 16),  # Increase legend label size
    legend.title = element_text(size = 14)  # Increase legend title size
  ) + NoLegend() +
  coord_cartesian(clip = "off")

um_colo
# Save BatchID UMAP
#png('plots/UMAP_colonocytes.png', height = 6, width = 8, res = 300, units = 'in')
svg('plots/UMAP_colonocytes.svg', height = 6, width = 8)
print(um_colo)
dev.off()

# UMAP for Condition
um_colo <- DimPlot(object = data_colo, reduction = "umap", group.by = "Condition", 
                   label = TRUE, repel = TRUE, label.size = 7, 
                   label.box = TRUE, alpha = 1, raster = TRUE) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    legend.text = element_text(size = 16),  # Increase legend label size
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) + NoLegend() +
  coord_cartesian(clip = "off")

um_colo
# Save BatchID UMAP
#png('plots/UMAP_colonocytes_condition.png', height = 6, width = 8, res = 300, units = 'in')
svg('plots/UMAP_colonocytes_condition.svg', height = 5, width = 7)
um_colo
dev.off()

# Apply custom theme
cust_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 20), # Center and enlarge title
  axis.text.x = element_text(angle = 60, hjust = 1, size = 18),
  axis.text.y = element_text(size = 14),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 18),
  strip.text = element_text(size = 16)
)


#data_colo$Condition <- factor(ifelse(grepl("WT", data_colo$BatchID), "WT", "Nr4a1 KO"))
# Create DittoSeq plot
cc <- dittoBarPlot(data_colo, var = 'SubCellType', group.by = 'Condition', scale = 'count') +
  cust_theme + ggtitle("Colonocyte cell counts")  # Set title correctly

cc
# Save plot
#png('plots/colonocytes_counts_subct.png', height = 6, width = 8, res = 300, units = 'in')
svg('plots/colonocytes_counts_subct.svg', height = 6, width = 8)
print(cc)  # Use 'cc' instead of 'um_batch'
dev.off()

#data_colo$Condition <- factor(ifelse(grepl("WT", data_colo$BatchID), "WT", "Nr4a1 KO"))
# Create DittoSeq plot
cc <- dittoBarPlot(data_colo, var = 'SubCellType', group.by = 'Condition', scale = 'percent') +
  cust_theme + ggtitle("Colonocyte cell counts (%)")  # Set title correctly

cc
# Save plot
svg('plots/colonocytes_counts_pct_subct.svg', height = 6, width = 8)
print(cc)  # Use 'cc' instead of 'um_batch'
dev.off()



##########
# Cell Scores for Colonocytes (Apoptosis, cell potency, etc)
##########
# Load the Excel file
library(readxl)
genelists <- read_excel('cellscores.xlsx')

# Define the score lists you are interested in
# score_lists <- c("Yuan et al (2019) Cell Cycle", "Yuan et al (2019) Differentiation",
#                  "Yuan et al (2019) Inflammation", "Yuan et al (2019) Proliferation",
#                  "Yuan et al (2019) Stemness", "Yuan et al (2019) Angiogenesis",
#                  "Yuan et al (2019) Apoptosis", "Stemness")

# Define the score lists you are interested in
score_lists <- c("Yuan et al (2019) Cell Cycle", 
                 "Yuan et al (2019) Inflammation",
                 "Yuan et al (2019) Stemness", "Stemness")

# Initialize an empty list to store the extracted data
extracted_data <- list()
# Initialize a list to store the UCell scores
all_scores <- list()

# Loop through each score list
for (score_type in score_lists) {
  # Create a logical index for the current score type
  idx <- genelists$ScoreType == score_type
  # Extract the rows corresponding to the current score type
  subset_df <- genelists[idx, ]
  # Initialize a vector to store the individual genes for this score type
  all_genes <- character(0)
  # Loop through the PositiveMarkers column of the subset
  for (genes_string in subset_df$PositiveMarkers) {
    # Split the comma-separated string into individual genes
    genes <- unlist(strsplit(genes_string, ","))
    # Trim any leading/trailing whitespace from the gene names
    genes <- trimws(genes)
    # Append the genes to the all_genes vector
    all_genes <- c(all_genes, genes)
  }
  # Convert all_genes to the format in data_colo (first letter uppercase, rest lowercase)
  formatted_all_genes <- sapply(tolower(all_genes), function(x) {
    if (nchar(x) > 0) {paste0(toupper(substring(x, 1, 1)), substring(x, 2))} else {""}
  })
  # Find the intersection with the rownames of data_colo (which are already in the correct format)
  intersected_genes <- intersect(formatted_all_genes, rownames(data_colo))
  # Create the gene signature list using the correctly cased intersected genes
  gene_signature <- list(signature = intersected_genes)
  # Check if there are any overlapping genes
  if (length(gene_signature$signature) > 0) {
    # Score the cells using UCell
    data_colo[[paste0("UCell_", score_type)]] <- AUCell_run(GetAssayData(data_colo),
                                                            gene_signature,
                                                            normAUC = TRUE)@assays@data@listData$AUC[1,]
    # Store the scores
    score_column_name <- paste0(score_type, "_UCell")
    all_scores[[score_column_name]] <- data_colo[[paste0("UCell_", score_type)]]
    cat(paste0("UCell scores calculated for: ", score_type, "\n"))
  } else {
    cat(paste0("No overlapping genes found for: ", score_type, ". Skipping UCell scoring.\n"))
    all_scores[[paste0(score_type, "_UCell")]] <- NA # Store NA if no overlap
  }
}


# Stemness Custom for gut
genes <- c("Lgr5", "Smoc2", "Olfm4", "Ascl2", "Sox9", "Cd44", "Ephb2", "Rnf43", "Notch1")
stemness_cust <- list(signature = intersect(genes, rownames(data_colo)))
data_colo[[paste0("UCell_", "Stemness_custom")]] <- AUCell_run(GetAssayData(data_colo),
                                                                                 stemness_cust,
                                                                                 normAUC = TRUE)@assays@data@listData$AUC[1,]
# Immune activation + stress response
genes <- c("Isg15", "Gbp2", "Saa1", "Cd74", "Atf3", "Hspa1b", "Gpx2", "Duox2", "Klf6", "Foxo3")
immune_act_stress_resp <- list(signature = intersect(genes, rownames(data_colo)))

data_colo[[paste0("UCell_", "Immune Activation|Stress Response")]] <- AUCell_run(GetAssayData(data_colo),
                                                                                   immune_act_stress_resp,
                                                                       normAUC = TRUE)@assays@data@listData$AUC[1,]

# Wnt beta catinent
genes <- c("FZD1","JAG1","JAG2","FZD8","WNT5B","AXIN2","AXIN1","PTCH1","CSNK1E","DKK1","DKK4","ADAM17","KAT2A",
                        "NCOR2","TCF7","PPARD","TP53","NUMB","HDAC5","LEF1","NOTCH4","NOTCH1","HDAC11","HDAC2","CTNNB1","MAML1","CUL1",
                        "PSEN2","DLL1","GNAI1","NKD1","RBPJ","WNT6","FRAT1","CCND2","HEY2","HEY1","DVL2","MYC","WNT1","NCSTN","SKP2")
# Format the custom inflammation gene list to match data_colo rownames
genes <- sapply(tolower(genes), function(x) {
  if (nchar(x) > 0) {
    paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  } else {
    ""
  }
})
wnt_beta_catin_sig <- list(signature = intersect(genes, rownames(data_colo)))

data_colo[[paste0("UCell_", "Wnt-beta Catenin Signaling - Hallmark")]] <- AUCell_run(GetAssayData(data_colo),
                                                                          wnt_beta_catin_sig,
                                                                          normAUC = TRUE)@assays@data@listData$AUC[1,]


# Establishment of Skin Barrier (GO:0061436)
genes <- c("Cldn4", "Stmn1")
skin_bar <- list(signature = intersect(genes, rownames(data_colo)))
data_colo[[paste0("UCell_", "Establishment of Skin Barrier")]] <- AUCell_run(GetAssayData(data_colo),
                                                                         skin_bar,
                                                                     normAUC = TRUE)@assays@data@listData$AUC[1,]



# Kegg Signaling pathways regulating pluripotency of stem cells
# The original human gene list:
# genes <- c("ACVR1", "PCGF5", "PCGF6", "WNT5A", "ESX1", "PCGF3", "WNT5B", "LIFR", "ISL1", "MAPK13", "MAPK14", "MAPK11", "MAPK12", "APC", "KAT6A", "IL6ST", "RAF1", "HOXB1", "FZD10", "NRAS", "AKT3", "MYC", "AKT1", "AKT2", "SMARCAD1", "FZD2", "ESRRB", "SMAD3", "SMAD2", "FZD1", "SMAD5", "MAP2K2", "FZD4", "WNT10B", "FZD3", "SMAD4", "WNT10A", "FZD6", "WNT3A", "FZD5", "MAP2K1", "SMAD1", "FZD8", "FZD7", "INHBA", "FZD9", "INHBC", "IGF1", "INHBB", "INHBE", "SMAD9", "HNF1A", "PIK3R1", "FGFR1", "NANOG", "GRB2", "LHX5", "TCF3", "FGFR4", "NODAL", "MYF5", "FGFR3", "FGFR2", "WNT2B", "RIF1", "FGF2", "IGF1R", "SOX2", "HESX1", "ZFHX3", "SETDB1", "AXIN2", "LIF", "AXIN1", "WNT9A", "WNT9B", "PAX6", "DUSP9", "ACVR2B", "ACVR2A", "WNT16", "BMP4", "MEIS1", "PIK3CB", "ID1", "PIK3CA", "ID3", "ID2", "ID4", "BMPR1A", "BMPR1B", "BMPR2", "GSK3B", "DLX5", "ONECUT1", "WNT8B", "PIK3R3", "PIK3R2", "WNT8A", "BMI1", "WNT6", "WNT11", "DVL1", "DVL2", "OTX1", "DVL3", "WNT1", "JARID2", "SKIL", "WNT2", "WNT3", "WNT4", "LEFTY2", "APC2", "STAT3", "WNT7A", "WNT7B", "POU5F1", "KLF4", "TBX3", "HAND1", "REST", "KRAS", "NEUROG1", "CTNNB1", "PIK3CD", "ACVR1C", "ACVR1B", "ZIC3", "MAPK3", "PCGF1", "PCGF2", "MAPK1", "HRAS", "JAK3", "JAK1", "JAK2")

# Converted to mouse gene names. Note: Some human genes may not have a direct mouse ortholog or may be represented differently.
genes <- c("Acvr1", "Pcgf5", "Pcgf6", "Wnt5a", "Esx1", "Pcgf3", "Wnt5b", "Lifr", 
           "Isl1", "Mapk13", "Mapk14", "Mapk11", "Mapk12", "Apc", "Kat6a", "Il6st", 
           "Raf1", "Hoxb1", "Fzd10", "Nras", "Akt3", "Myc", "Akt1", "Akt2", "Smarcad1", 
           "Fzd2", "Esrrb", "Smad3", "Smad2", "Fzd1", "Smad5", "Map2k2", "Fzd4", 
           "Wnt10b", "Fzd3", "Smad4", "Wnt10a", "Fzd6", "Wnt3a", "Fzd5", "Map2k1", 
           "Smad1", "Fzd8", "Fzd7", "Inhba", "Fzd9", "Inhbc", "Igf1", "Inhbb", 
           "Inhbe", "Smad9", "Hnf1a", "Pik3r1", "Fgfr1", "Nanog", "Grb2", "Lhx5", 
           "Tcf3", "Fgfr4", "Nodal", "Myf5", "Fgfr3", "Fgfr2", "Wnt2b", "Rif1", 
           "Fgf2", "Igf1r", "Sox2", "Hesx1", "Zfhx3", "Setdb1", "Axin2", "Lif", 
           "Axin1", "Wnt9a", "Wnt9b", "Pax6", "Dusp9", "Acvr2b", "Acvr2a", "Wnt16", 
           "Bmp4", "Meis1", "Pik3cb", "Id1", "Pik3ca", "Id3", "Id2", "Id4", 
           "Bmpr1a", "Bmpr1b", "Bmpr2", "Gsk3b", "Dlx5", "Onecut1", "Wnt8b", 
           "Pik3r3", "Pik3r2", "Wnt8a", "Bmi1", "Wnt6", "Wnt11", "Dvl1", "Dvl2", 
           "Otx1", "Dvl3", "Wnt1", "Jarid2", "Skil", "Wnt2", "Wnt3", "Wnt4", 
           "Lefty2", "Apc2", "Stat3", "Wnt7a", "Wnt7b", "Pou5f1", "Klf4", "Tbx3", 
           "Hand1", "Rest", "Kras", "Neurog1", "Ctnnb1", "Pik3cd", "Acvr1c", "Acvr1b", 
           "Zic3", "Mapk3", "Pcgf1", "Pcgf2", "Mapk1", "Hras", "Jak3", "Jak1", "Jak2")

kegg_pluripotency <- list(signature = intersect(genes, rownames(data_colo)))
data_colo[[paste0("UCell_", "Kegg_pluripotency")]] <- AUCell_run(GetAssayData(data_colo),
                                                                 kegg_pluripotency,
                                                                 normAUC = TRUE)@assays@data@listData$AUC[1,]

library(msigdbr)
msigdbr_df <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
pathways_gobp <- split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
available_pathways <- unique(msigdbr_df$gs_name)
wnt_pathways <- available_pathways[grepl("WNT_SIGNALING", available_pathways, ignore.case = TRUE)]
can_wnt_s <- pathways_gobp["GOBP_CANONICAL_WNT_SIGNALING_PATHWAY"]
ncan_wnt_s <- pathways_gobp["GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY"]
wnt_s <- pathways_gobp["GOBP_WNT_SIGNALING_PATHWAY"]
#pathways_gobp["GOBP_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY"]
#pathways_gobp["GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY"]
#pathways_gobp["GOBP_REGULATION_OF_NON_CANONICAL_WNT_SIGNALING_PATHWAY"]

genes <- list(signature = intersect(ncan_wnt_s, rownames(data_colo)))
data_colo[[paste0("UCell_", "Non Canonical WNT Signaling Pathway - GOBP")]] <- AUCell_run(GetAssayData(data_colo),
                                                                         ncan_wnt_s,
                                                                         normAUC = TRUE)@assays@data@listData$AUC[1,]
genes <- list(signature = intersect(can_wnt_s, rownames(data_colo)))
data_colo[[paste0("UCell_", "Canonical WNT Signaling Pathway - GOBP")]] <- AUCell_run(GetAssayData(data_colo),
                                                                                ncan_wnt_s,
                                                                                normAUC = TRUE)@assays@data@listData$AUC[1,]
genes <- list(signature = intersect(wnt_s, rownames(data_colo)))
data_colo[[paste0("UCell_", "WNT Signaling Pathway - GOBP")]] <- AUCell_run(GetAssayData(data_colo),
                                                                                ncan_wnt_s,
                                                                                normAUC = TRUE)@assays@data@listData$AUC[1,]



# Apply custom theme
cust_theme <- theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 16), # Center and enlarge title
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
  axis.title.x = element_blank(),
  axis.text.y = element_text(size = 15),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 13),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 14)
)

# Get the names of the UCell score columns
score_columns <- grep("^UCell_", colnames(data_colo@meta.data), value = TRUE)

# Loop through each score column to generate and save a plot
for (score_col in score_columns) {
  # Extract the base score name for the title and filename
  score_name <- gsub("^UCell_", "", score_col)
  clean_score_name <- gsub("Yuan et al \\(2019\\) ", "", score_name)
  
  # Create the plot
  p <- data_colo[[c(score_col, 'Condition', 'SubCellType')]] %>%
    rename(score = !!sym(score_col)) %>% # Rename the score column for easier use
    ggplot(aes(x = Condition, y = score, fill = Condition)) +
    geom_violin() +
    geom_boxplot(color = 'black', width = 0.5, outlier.size = 0.3) +
    geom_signif(comparisons = list(c("WT", "Nr4a1 KO")),
                map_signif_level = TRUE, step_increase = 0.12, vjust = 2) +
    ggtitle(clean_score_name) + ylab('AUCell score') +
    facet_wrap(~SubCellType,  ncol = 4,  scales = 'free') + cust_theme + NoLegend()
  
  # Define the filename for the SVG
  filename <- paste0('plots/', gsub("[^[:alnum:]]", "_", score_name), '_Condition_subct_colonocytes.svg')
  
  # Save the plot as an SVG
  svg(filename, height = 4, width = 10) #5 cols
  #svg(filename, height = 6, width = 10) #4 cols
  
  print(p) # Important: Use print() to draw the ggplot object to the device
  dev.off()
  
  cat(paste0("Saved plot: ", filename, "\n"))
}

####################################
# Plotting differentiation potency
# Extract the base score name for the title and filename
score_name <- "CellPotency"
clean_score_name <- gsub("Yuan et al \\(2019\\) ", "", score_name)
score_col <- score_name

# Create the plot
p <- data_colo[[c(score_col, 'Condition', 'SubCellType')]] %>%
  rename(score = !!sym(score_col)) %>% # Rename the score column for easier use
  ggplot(aes(x = Condition, y = score, fill = Condition)) +
  geom_violin() +
  geom_boxplot(color = 'black', width = 0.5, outlier.size = 0.3) +
  geom_signif(comparisons = list(c("WT", "Nr4a1 KO")),
              map_signif_level = TRUE, step_increase = 0.12, vjust = 2) +
  ggtitle("Differentiation Potency") + ylab('Score') +
  facet_wrap(~SubCellType,  ncol = 4,  scales = 'free') + cust_theme + NoLegend()

# Define the filename for the SVG
filename <- paste0('plots/', gsub("[^[:alnum:]]", "_", score_name), '_Condition_subct_colonocytes.svg')

# Save the plot as an SVG
svg(filename, height = 4, width = 10) #5 cols
#svg(filename, height = 6, width = 10) #4 cols

print(p) # Important: Use print() to draw the ggplot object to the device
dev.off()

cat("Finished saving SVG plots for each score.\n")



# Load necessary libraries (if not already loaded)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer, pivot_wider
library(ggpubr) # For geom_signif (though not used directly in this specific bar plot style, good to have if needed for error bars)
library(Hmisc) # For mean_cl_normal (if error bars were to be added to this specific bar plot)
library(writexl)

# --- Define the score name ---
score_name <- "CellPotency"
clean_score_name <- gsub("Yuan et al \\(2019\\) ", "", score_name) # Adjust if your score name doesn't have this prefix

# --- 1. Extract and prepare data for CellPotency score ---
# Ensure these columns exist in your data_colo@meta.data
if (!all(c(score_name, 'Condition', 'SubCellType') %in% colnames(data_colo@meta.data))) {
  stop("One or more required metadata columns (CellPotency, Condition, SubCellType) not found in data_colo.")
}

# Create a data frame from Seurat object metadata
df_cell_potency <- data_colo[[c(score_name, 'Condition', 'SubCellType')]] %>%
  rename(score = !!sym(score_name)) # Rename the score column for easier use

# Ensure Condition is a factor with desired levels
df_cell_potency$Condition <- factor(df_cell_potency$Condition, levels = c("Nr4a1 KO", "WT"))
# Ensure SubCellType is a factor (order can be set later if desired)
df_cell_potency$SubCellType <- factor(df_cell_potency$SubCellType)


# --- 2. Perform Wilcoxon test and summarize data for CellPotency ---
results_table <- df_cell_potency %>%
  group_by(SubCellType) %>% # Group by SubCellType, equivalent to CellType in your example
  summarise(
    # Get statistics for KO
    N_KO = sum(Condition == "Nr4a1 KO"),
    Mean_KO = mean(score[Condition == "Nr4a1 KO"], na.rm = TRUE),
    # Get statistics for WT
    N_WT = sum(Condition == "WT"),
    Mean_WT = mean(score[Condition == "WT"], na.rm = TRUE),
    # Perform Wilcoxon test
    p_value = {
      ko_values <- score[Condition == "Nr4a1 KO"]
      wt_values <- score[Condition == "WT"]
      
      # Check if both KO and WT conditions exist and have enough data (at least 2 observations for wilcox.test)
      if (length(ko_values) >= 2 && length(wt_values) >= 2) {
        wilcox_result <- tryCatch(
          wilcox.test(ko_values, wt_values, exact = FALSE)$p.value,
          error = function(e) NA # Handle cases where test still can't be performed (e.g., all values identical)
        )
        wilcox_result
      } else {
        NA # Not enough data for the test (e.g., <2 cells in a group)
      }
    },
    .groups = 'drop' # Drop grouping after summarise
  ) %>%
  ungroup() # Ungroup the data frame

# Display the initial results table
cat("\n--- CellPotency Score Wilcoxon Test Results by SubCellType ---\n")
print(results_table)


# --- 3. Format the p-value for better readability and align with geom_signif ---
results_table_formatted <- results_table %>%
  mutate(p_value_formatted = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    is.na(p_value) ~ "N/A", # For cases where test couldn't be run
    TRUE ~ "NS"
  ))

cat("\n--- Formatted CellPotency Score Wilcoxon Test Results ---\n")
print(results_table_formatted)
write_xlsx(results_table_formatted, "cell_potency_score_colonocytes_stats.xlsx")

# --- If you want to recreate final_summary_table with the correct p-values (as in your example): ---
final_summary_table <- results_table_formatted %>%
  select(SubCellType, N_KO, N_WT, Mean_KO, Mean_WT, p_value, p_value_formatted) %>%
  pivot_longer(
    cols = starts_with("N_") | starts_with("Mean_"),
    names_to = c(".value", "Condition"),
    names_pattern = "([A-Za-z]+)_(KO|WT)"
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(N, Mean),
    names_glue = "{.value}_{Condition}"
  ) %>%
  select(SubCellType, N_KO, N_WT, Mean_KO, Mean_WT, p_value, p_value_formatted)

cat("\n--- Final Summary Table for CellPotency Score ---\n")
print(final_summary_table)


#################
### Cell Cycle score
################
s_genes <- c("Atad2","Brip1","Casp8ap2","Ccne2","Cdc45","Cdc6","Cdca7","Cenpu","Chaf1b","Clspn","Dscc1","Dtl","E2f8","Exo1","Fen1","Gins2","Gmnn","Hells","Mchm2","Mcm4","Mcm5","Mcm6","Msh2","Nasp","Pcna","Pola1","Pold3","Polr1b","Prim1","Rad51ap1","Rfc2","Rrm1","Rrm2","Slbp","Tipin","Tyms","Ubr7","Uhrf1","Ung","Usp1","Vps50","Wdr76","Xrcc2")
g2m_genes <- c("Anln","Anp32e","Aurka","Aurkb","Birc5","Bub1","Cbx5","Ccnb2","Cdc20","Cdc25c","Cdca2","Cdca3","Cdca8","Cdk1","Cenpa","Cenpe","Cenpf","Ckap2","Ckap2l","Ckap5","Cks1b","Cks2","Ctcf","Dazl","Dlgap5","Ect2","G2e3","Gas2l3","Gtse1","Hjurp","Hmgb2","Hmmr","Jpt1","Kif11","Kif20b","Kif23","Kif2c","Lbr","Mki67","Ncapd2","Ndc80","Nek2","Nuf2")
data_colo <- CellCycleScoring(object = data_colo, 
                               s.features = s_genes, 
                               g2m.features = g2m_genes)
DimPlot(data_colo, group.by = 'Phase')

library(dplyr)
library(ggplot2)
# Explicitly convert the metadata to a data frame (if needed)
metadata_df <- as.data.frame(data_colo@meta.data)

# Use base R bracket notation to select the columns
cell_data <- metadata_df[, c("Phase", "Condition", "CellType", "SubCellType", "S.Score", "G2M.Score")]

cell_data_filtered <- cell_data
# Calculate counts and percentages
plot_data <- cell_data_filtered %>%
  group_by(Condition, SubCellType, Phase) %>%
  summarise(CellCount = n(), .groups = 'drop') %>%
  group_by(Condition, SubCellType) %>%
  mutate(Percentage = CellCount / sum(CellCount) * 100)

p <- ggplot(plot_data, aes(x = Condition, y = Percentage, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 4) +
  facet_wrap(~SubCellType, scales = "free_x", ncol = 4) +
  labs(title = "Cell Cycle Phase Percentage by Condition and SubPopulations",
       x = "Condition",
       y = "Percentage of Cells") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 10, color = "black"),  # Keep axis titles
        legend.text = element_text(size = 12),  # Increase legend label size
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12)
        )
p
ggsave("plots/CellCycle_proportions_colonocytes.svg", plot = p, height = 5, width = 7, device = "svg")


library(tidyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)

# Pivot the data to a long format
cell_data_filtered <- cell_data
plot_data_long <- cell_data_filtered %>%
  pivot_longer(
    cols = c("S.Score", "G2M.Score"),
    names_to = "Score_Type",
    values_to = "Score_Value"
  )

# Create the combined violin plot
p <- ggplot(plot_data_long, aes(x = Condition, y = Score_Value, fill = Condition)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7, outlier.size = 0.3) +
  facet_wrap(Phase ~ SubCellType, scales = "free_y", ncol = 7) +
  labs(
    title = "Cell Cycle Score Distribution by Subcell Type and Condition",
    x = "SubCellType",
    y = "Score Value"
  ) +
  geom_signif(comparisons = list(c('Nr4a1 KO','WT')), map_signif_level = T,
              vjust = 2, margin_top = 0.2)+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

p
ggsave("plots/CellCycle_scores_colonocytes.svg", plot = p, height = 7, width = 10, device = "svg")


# Filter only for stem cells 
# Filter the data to include only SubCellType == "Stem"
cell_data_filtered <- cell_data %>%
  filter(SubCellType == "Stem")

# Calculate counts and percentages
plot_data <- cell_data_filtered %>%
  group_by(Condition, SubCellType, Phase) %>%
  summarise(CellCount = n(), .groups = 'drop') %>%
  group_by(Condition, SubCellType) %>%
  mutate(Percentage = CellCount / sum(CellCount) * 100)

# Generate the plot
ggplot(plot_data, aes(x = Condition, y = Percentage, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  facet_grid(SubCellType ~ Condition, scales = "free_x", space = "free_x") +
  labs(title = "Cell Cycle Phase Percentage for Stem Cells by Condition, and Cell Type",
       x = "Condition",
       y = "Percentage of Cells") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(tidyr)
library(ggplot2)
library(ggsignif)
# Pivot the data to a long format
plot_data_long <- cell_data_filtered %>%
  pivot_longer(
    cols = c("S.Score", "G2M.Score"),
    names_to = "Score_Type",
    values_to = "Score_Value"
  )

# Create the combined violin plot
ggplot(plot_data_long, aes(x = Condition, y = Score_Value, fill = Condition)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7) +
  facet_wrap(~Score_Type, scales = "free_y") +
  labs(
    title = "Cell Cycle Score Distribution in Stem Cells",
    x = "Condition",
    y = "Score Value"
  ) +
  geom_signif(comparisons = list(c("WT", "Nr4a1 KO")), test = "wilcox.test",
              map_signif_level = TRUE, step_increase = 0.12, vjust = 2) +
  cust_theme + NoLegend()



####################################
# MANAscore 2025 stemness score table | or Stemness_custom
####################################
# Load necessary libraries (if not already loaded)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer, pivot_wider
library(ggpubr) # For geom_signif (though not used directly in this specific bar plot style, good to have if needed for error bars)
library(Hmisc) # For mean_cl_normal (if error bars were to be added to this specific bar plot)
library(writexl)

# --- Define the score name ---
score_name <- "UCell_Stemness_custom"
clean_score_name <- gsub("Yuan et al \\(2019\\) ", "", score_name) # Adjust if your score name doesn't have this prefix

# --- 1. Extract and prepare data for CellPotency score ---
# Ensure these columns exist in your data_colo@meta.data
if (!all(c(score_name, 'Condition', 'SubCellType') %in% colnames(data_colo@meta.data))) {
  stop("One or more required metadata columns (UCell_Stemness, Condition, SubCellType) not found in data_colo.")
}

# Create a data frame from Seurat object metadata
df_cell_potency <- data_colo[[c(score_name, 'Condition', 'SubCellType')]] %>%
  rename(score = !!sym(score_name)) # Rename the score column for easier use

# Ensure Condition is a factor with desired levels
df_cell_potency$Condition <- factor(df_cell_potency$Condition, levels = c("Nr4a1 KO", "WT"))
# Ensure SubCellType is a factor (order can be set later if desired)
df_cell_potency$SubCellType <- factor(df_cell_potency$SubCellType)


# --- 2. Perform Wilcoxon test and summarize data for CellPotency ---
results_table <- df_cell_potency %>%
  group_by(SubCellType) %>% # Group by SubCellType, equivalent to CellType in your example
  summarise(
    # Get statistics for KO
    N_KO = sum(Condition == "Nr4a1 KO"),
    Mean_KO = mean(score[Condition == "Nr4a1 KO"], na.rm = TRUE),
    # Get statistics for WT
    N_WT = sum(Condition == "WT"),
    Mean_WT = mean(score[Condition == "WT"], na.rm = TRUE),
    # Perform Wilcoxon test
    p_value = {
      ko_values <- score[Condition == "Nr4a1 KO"]
      wt_values <- score[Condition == "WT"]
      
      # Check if both KO and WT conditions exist and have enough data (at least 2 observations for wilcox.test)
      if (length(ko_values) >= 2 && length(wt_values) >= 2) {
        wilcox_result <- tryCatch(
          wilcox.test(ko_values, wt_values, exact = FALSE)$p.value,
          error = function(e) NA # Handle cases where test still can't be performed (e.g., all values identical)
        )
        wilcox_result
      } else {
        NA # Not enough data for the test (e.g., <2 cells in a group)
      }
    },
    .groups = 'drop' # Drop grouping after summarise
  ) %>%
  ungroup() # Ungroup the data frame

# Display the initial results table
cat("\n--- UCell_Stemness Score Wilcoxon Test Results by SubCellType ---\n")
print(results_table)

# --- 3. Format the p-value for better readability and align with geom_signif ---
results_table_formatted <- results_table %>%
  mutate(p_value_formatted = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    is.na(p_value) ~ "N/A", # For cases where test couldn't be run
    TRUE ~ "NS"
  ))

cat("\n--- Formatted Stemness Score Wilcoxon Test Results ---\n")
print(results_table_formatted)
write_xlsx(results_table_formatted, "stemness_score_colonocytes_stats.xlsx")

# --- If you want to recreate final_summary_table with the correct p-values (as in your example): ---
final_summary_table <- results_table_formatted %>%
  select(SubCellType, N_KO, N_WT, Mean_KO, Mean_WT, p_value, p_value_formatted) %>%
  pivot_longer(
    cols = starts_with("N_") | starts_with("Mean_"),
    names_to = c(".value", "Condition"),
    names_pattern = "([A-Za-z]+)_(KO|WT)"
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(N, Mean),
    names_glue = "{.value}_{Condition}"
  ) %>%
  select(SubCellType, N_KO, N_WT, Mean_KO, Mean_WT, p_value, p_value_formatted)

cat("\n--- Final Summary Table for CellPotency Score ---\n")
print(final_summary_table)



###############################
# --- PLOTTING SECTION ---
###############################
# Define custom theme for consistent plotting style (as provided by you)
cust_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title.x = element_blank(), # Remove x-axis title
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1), # Rotate x-axis labels
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray80"),
    panel.grid.minor.y = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # Add some margin around the plot
  )

NoLegend <- function() {
  theme(legend.position = "none")
}


# --- Prepare data for combined facet labels AND order them ---
# Get the desired order of SubCellTypes from your final_summary_table
desired_subcelltype_order <- final_summary_table$SubCellType

facet_annotation_data <- final_summary_table %>%
  # Ensure SubCellType is a factor with the desired order *before* creating Facet_Label
  mutate(
    SubCellType = factor(SubCellType, levels = desired_subcelltype_order),
    Facet_Label = paste0(SubCellType, "\nW. p-val: ", p_value_formatted)
  ) %>%
  # Convert Facet_Label to a factor with levels ordered by SubCellType
  mutate(
    Facet_Label = factor(Facet_Label, levels = unique(Facet_Label[order(SubCellType)]))
  ) %>%
  select(SubCellType, Facet_Label, p_value_formatted)


# 1. Transform the data into a long format for plotting, carrying N values
plot_data <- final_summary_table %>%
  # Ensure SubCellType is a factor with the desired order *before* pivoting
  mutate(
    SubCellType = factor(SubCellType, levels = desired_subcelltype_order)
  ) %>%
  pivot_longer(
    cols = c(Mean_KO, Mean_WT, N_KO, N_WT),
    names_to = c(".value", "Condition"),
    names_pattern = "([A-Za-z]+)_(KO|WT)",
    values_drop_na = FALSE
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("KO", "WT")) # Ensure correct factor order for Condition
  ) %>%
  rename(Mean_Score = Mean, N_Cells = N) %>% # Renamed to Mean_Score for clarity with CellPotency
  select(SubCellType, Condition, Mean_Score, N_Cells)

# Now, crucial for geom_text for p-value stars:
# Make sure plot_data also has the Facet_Label with the correct factor order,
# so that geom_text can inherit it correctly for faceting.
plot_data <- plot_data %>%
  left_join(select(facet_annotation_data, SubCellType, Facet_Label), by = "SubCellType") %>%
  # Ensure Facet_Label in plot_data is also a factor with the correct order
  mutate(Facet_Label = factor(Facet_Label, levels = levels(facet_annotation_data$Facet_Label)))


# 2. Create the bar plots with mean values and significance stars
p_cellpotency_bar_annotated <- plot_data %>%
  ggplot(aes(x = Condition, y = Mean_Score, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~ Facet_Label, ncol = 4, scales = 'free_y') + # Facet by Facet_Label
  labs(
    y = paste0(clean_score_name, " Mean Score"), # Dynamic y-axis label
    x = "Condition"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", margin = margin(b = 10, t = 5), size = 14),
    strip.background = element_rect(fill = "white", color = NA),
    
    panel.spacing.y = unit(0.5, "lines"),
    panel.spacing.x = unit(1, "lines"),
    
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
    
    axis.title.y = element_text(margin = margin(r = 2), size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_blank(), # This is the correct and desired line for x-axis text
    
    legend.position = "top",
    legend.box.margin = margin(b = 2),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  scale_fill_manual(values = c("KO" = scales::alpha("red", 0.6),
                               "WT" = scales::alpha("blue", 0.6))) +
  # Increase y range using expansion to ensure labels are visible
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  # Add Mean Value above each bar
  geom_text(
    aes(
      label = sprintf("%.3f", Mean_Score) # Use Mean_Score
    ),
    vjust = -0.3,
    position = position_dodge(width = 0.8),
    size = 5,
    color = "black"
  ) +
  # ADDED BACK: geom_text FOR CELL COUNTS (N_Cells)
  geom_text(
    aes(
      label = paste0("N=", N_Cells), # Format as "N=X"
      y = Mean_Score # Start N at the top of the bar where Mean_Score is
    ),
    vjust = 1.5, # Adjust vertical justification relative to Mean_Score: 1.5 moves it below Mean_Score
    position = position_dodge(width = 0.8),
    size = 4, # Slightly smaller size for count
    color = "black"
  )

print(p_cellpotency_bar_annotated)

# Ensure 'plots' directory exists
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Save plot
filename_svg <- paste0('plots/', gsub("[^[:alnum:]]", "_", clean_score_name), '_MeanBarPlot_Annotated_Condition_subct.svg')
svg(filename_svg, height = 5, width = 10)
print(p_cellpotency_bar_annotated)
dev.off()

cat(paste0("Finished saving SVG plot for ", clean_score_name, " score with annotations.\n"))




########################################
## Differential expression analysis 
########################################

# Ensure `cell_types` is properly extracted
cell_types <- unique(data_colo$SubCellType)
# Set the grouping variable for Seurat
Idents(data_colo) <- "Condition"
# Define output directory
output_dir <- "DE_results_genes/colonocytes"
# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (cell in cell_types) {
  # Print the current cell type being processed
  cat("Processing cell type:", cell, "\n")
  
  # Subset the Seurat object for the current cell type
  cell_subset <- subset(data_colo, SubCellType == cell)
  
  # Perform FindMarkers for this cell type comparing WT vs Nr4a1 KO
  markers <- FindMarkers(cell_subset,
                         ident.1 = "Nr4a1 KO",
                         ident.2 = "WT",
                         test.use = "wilcox",  # Wilcoxon rank-sum test
                         logfc.threshold = 0.1,
                         min.pct = 0.01,
                         min.diff.pct = -Inf,
                         verbose = TRUE)
  
  # Filter significant markers
  #all_markers <- markers[markers$p_val < 0.05, ]
  all_markers <- markers[markers$p_val_adj < 0.05, ]
  
  # Sort by absolute value of log2FC
  all_markers <- all_markers[order(abs(all_markers$avg_log2FC), decreasing = TRUE), ]
  
  upregulated <- all_markers[all_markers$avg_log2FC > 0, ]
  downregulated <- all_markers[all_markers$avg_log2FC < 0, ]
  
  # Add gene names as a column
  all_markers$Gene <- rownames(all_markers)
  upregulated$Gene <- rownames(upregulated)
  downregulated$Gene <- rownames(downregulated)
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add sheets to the workbook
  addWorksheet(wb, "All_Pval_0.05")
  addWorksheet(wb, "Up_regulated")
  addWorksheet(wb, "Down_regulated")
  
  # Write data to respective sheets
  writeData(wb, sheet = "All_Pval_0.05", all_markers)
  writeData(wb, sheet = "Up_regulated", upregulated)
  writeData(wb, sheet = "Down_regulated", downregulated)
  
  # Ensure consistent file naming
  save_cell_name <- gsub(" ", "_", cell)  # Replace spaces with underscores *once*
  output_file <- file.path(output_dir, paste0(save_cell_name, "_markers_Wilcoxon.xlsx"))
  
  # Save workbook
  saveWorkbook(wb, file = output_file, overwrite = TRUE)
  
  cat("Saved:", output_file, "\n")
}





library(enrichR)
library(openxlsx)
library(dplyr)

# Define paths
input_dir <- "DE_results_genes/colonocytes"
output_dir <- "Enrichr_results/colonocytes"

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set Enrichr database and parameters
setEnrichrSite('Enrichr')
#hallmark_dbs_stringr <- available_dbs[str_detect(available_dbs$libraryName, regex("hallmark", ignore_case = TRUE)), ]
#dbs <- "MSigDB_Hallmark_2020"
dbs <- "GO_Biological_Process_2025"

pval_cutoff <- 0.05
top_n <- 200

# Function to perform Enrichr and process results with waiting time
run_enrichr <- function(genes, direction, dbs, pval_cutoff, wait_time = 2) {
  results_df <- data.frame()
  
  for (db in dbs) {
    if (length(genes) > 0) { # Check if genes vector is not empty
      enrichr_res_list <- enrichr(genes, db)
      if (!is.null(enrichr_res_list) && length(enrichr_res_list) > 0) {
        enrichr_res <- enrichr_res_list[[1]]
        
        if (!is.null(enrichr_res) && nrow(enrichr_res) > 0) {
          #enrichr_res <- enrichr_res %>% filter(P.value < pval_cutoff)
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
    } else {
      print(paste("No genes provided for", direction))
    }
    Sys.sleep(wait_time)  # Prevent API overload
  }
  return(results_df)
}

# Assuming 'data_colo' is a dataframe containing T cell data
# and has a column named 'SubCellType'

# Get unique cell types
cell_types <- unique(data_colo$SubCellType)

for (cell in cell_types) {
  print(paste("Working with", cell))
  file_path <- file.path(input_dir, paste0(gsub(" ", "_", cell), "_markers_Wilcoxon.xlsx"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  
  wb <- loadWorkbook(file_path)
  upregulated <- read.xlsx(wb, sheet = "Up_regulated")
  downregulated <- read.xlsx(wb, sheet = "Down_regulated")
  
  top_upregulated <- upregulated %>% head(top_n) %>% pull(Gene)
  top_downregulated <- downregulated %>% head(top_n) %>% pull(Gene)
  
  enrichr_results_up <- run_enrichr(top_upregulated, paste0("Up-regulated ", cell), dbs, pval_cutoff)
  enrichr_results_dn <- run_enrichr(top_downregulated, paste0("Down-regulated ", cell), dbs, pval_cutoff)
  
  # Skip file writing if no significant pathways are found
  if (nrow(enrichr_results_up) == 0 && nrow(enrichr_results_dn) == 0) {
    print(paste("Skipping file for", cell, "- No significant pathways found"))
    next
  }
  
  # Write results to Excel
  output_file <- file.path(output_dir, paste0(gsub("[^[:alnum:]_]", "", gsub(" ", "_", cell)), "_enrichr_Wilcoxon_results.xlsx"))
  wb_out <- createWorkbook()
  
  if (nrow(enrichr_results_up) > 0) {
    addWorksheet(wb_out, "Upregulated_Pathways")
    writeData(wb_out, "Upregulated_Pathways", enrichr_results_up)
  }
  
  if (nrow(enrichr_results_dn) > 0) {
    addWorksheet(wb_out, "Downregulated_Pathways")
    writeData(wb_out, "Downregulated_Pathways", enrichr_results_dn)
  }
  
  # Only save workbook if at least one pathway exists
  saveWorkbook(wb_out, output_file, overwrite = TRUE)
  print(paste("Enrichment results saved in", output_file))
}






# Define custom theme for consistent plotting style
cust_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_blank(), # Remove x-axis title
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1), # Rotate x-axis labels
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray80"),
    panel.grid.minor.y = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # Add some margin around the plot
  )

# Function to remove legend (for consistent use with your previous code)
NoLegend <- function() {
  theme(legend.position = "none")
}

# List of genes to plot D GRN
genes_to_plot <- c("Ascl2", "Wnt2b", "Id4", "Ephb3", "Melk", "Birc5",
                   "Bhlha15", "Gcg", "Insl5", "Clca1", "Fcgbp", "Smoc2", "Gas6",
                   "Atf3", "Clu", "Sgk1")

# Define the target cell type subset (e.g., "Stem")
ct_sub <- "Stem"

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each gene and create a violin plot
for (gene in genes_to_plot) {
  # Check if the gene exists in the data_colo object
  if (gene %in% rownames(data_colo)) {
    # Extract expression data for the current gene, specifically for cells in 'ct_sub'
    # This pulls from the 'data' slot (normalized and/or scaled data)
    gexp <- data_colo@assays$RNA$data[gene, data_colo$SubCellType == ct_sub]
    # Extract corresponding 'Condition' metadata for the same subset of cells
    cond <- data_colo$Condition[data_colo$SubCellType == ct_sub]
    
    # Create a data frame for ggplot
    df_plot <- data.frame(
      Expression = as.numeric(gexp),
      Condition = cond # 'Condition' metadata for the subsetted cells
    )
    
    # Ensure Condition is a factor with desired levels for consistent plotting order
    df_plot$Condition <- factor(df_plot$Condition, levels = c("Nr4a1 KO", "WT"))
    
    # --- Dynamic Y-limit Calculation ---
    # Calculate the maximum expression for the current gene to set a dynamic y-limit.
    # We add a small buffer (e.g., 10%) to the max value for better visualization.
    max_expression_for_gene <- max(df_plot$Expression) * 1.3
    
    # Determine the minimum y-limit:
    # If all expression values are non-negative, start the axis at 0.
    # If there are negative values (e.g., from certain scaling methods), extend the min limit.
    min_expression_for_gene <- min(df_plot$Expression)
    if (min_expression_for_gene < 0) {
      min_expression_for_gene <- min_expression_for_gene * 1.1 # Extend min if negative
    } else {
      min_expression_for_gene <- 0 # Start from 0 if all values are non-negative
    }
    
    # Add a small buffer to the max if the range is very narrow (e.g., all values are similar)
    # This prevents the plot from looking squashed if expression is very flat.
    if (abs(max_expression_for_gene - min_expression_for_gene) < 0.1) {
      max_expression_for_gene <- max_expression_for_gene + 0.5
    }
    # Ensure min is not greater than max, and set a default if max is still problematic
    if (max_expression_for_gene <= min_expression_for_gene) {
      max_expression_for_gene <- min_expression_for_gene + 1 # Ensure at least a range of 1
    }
    
    # Create the violin plot
    p <- df_plot %>%
      ggplot(aes(y = Expression, x = Condition, fill = Condition)) +
      geom_violin(trim = TRUE, color = "black", alpha = 0.7) + # Add black border and transparency
      geom_boxplot(width = 0.4, outlier.shape = NA, fill = "white", alpha = 0.7) + # Boxplot on top, no outliers
      # geom_point(aes(color = Condition), size = 0.5, position = position_jitter(width = 0.1), alpha = 0.3) + # Add jittered points (commented out as per user request)
      # Add significance bars
      geom_signif(map_signif_level = TRUE, comparisons = list(c("Nr4a1 KO", "WT")),
                  vjust = -0.3, tip_length = 0.02, # Adjust vjust and tip_length for better appearance
                  test = "wilcox.test", test.args = list(exact = FALSE), # exact = FALSE for larger datasets
                  textsize = 3) + # Adjust text size for p-value
      labs(title = paste0(gene, " - ", ct_sub, " Cells"), # Dynamically include cell type in title
           y = "Expression Level") +
      scale_fill_manual(values = c("Nr4a1 KO" = "#E41A1C", "WT" = "#377EB8")) + # Custom colors
      scale_color_manual(values = c("Nr4a1 KO" = "#E41A1C", "WT" = "#377EB8")) + # Custom colors for points
      cust_theme +
      NoLegend() +
      # Dynamically set y-axis limits for each plot using coord_cartesian
      coord_cartesian(ylim = c(min_expression_for_gene, max_expression_for_gene))
    
    # Add the plot to the list
    plot_list[[gene]] <- p
  } else {
    message(paste0("Gene '", gene, "' not found in the Seurat object. Skipping plot."))
  }
}

# --- Displaying the plots ---
# Combine all plots into a single figure
p1 <- wrap_plots(plot_list, ncol = 4)
print(p1)

# Save plot
svg('plots/diff_grn_exp_cond_stem.svg', height = 10, width = 12)
print(p1)
dev.off()
# To save all plots to a PDF file:
# pdf("Stem_Cell_Gene_Violin_Plots.pdf", width = 10, height = 8 * ceiling(length(plot_list) / 4)) # Adjust width/height for 4 columns
# for (p in plot_list) {
#   print(p)
# }
# dev.off()




# Define custom theme for consistent plotting style
cust_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title.x = element_blank(), # Remove x-axis title
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1), # Rotate x-axis labels
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray80"),
    panel.grid.minor.y = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # Add some margin around the plot
  )

# Function to remove legend (for consistent use with your previous code)
NoLegend <- function() {
  theme(legend.position = "none")
}


# List of genes to plot
genes_to_plot <- c("Ascl2", "Wnt2b", "Id4", "Ephb3", "Melk", "Birc5",
                   "Bhlha15", "Gcg", "Insl5", "Clca1", "Fcgbp", "Smoc2", "Gas6",
                   "Atf3", "Clu", "Sgk1")

# Define the target cell type subset (e.g., "Stem")
ct_sub <- "Stem"

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each gene and create a bar plot of mean expression with confidence intervals
for (gene in genes_to_plot) {
  # Check if the gene exists in the data_colo object
  if (gene %in% rownames(data_colo)) {
    # Extract expression data for the current gene, specifically for cells in 'ct_sub'
    # This pulls from the 'data' slot (normalized and/or scaled data)
    gexp <- data_colo@assays$RNA$data[gene, data_colo$SubCellType == ct_sub]
    # Extract corresponding 'Condition' metadata for the same subset of cells
    cond <- data_colo$Condition[data_colo$SubCellType == ct_sub]
    
    # Create a data frame for ggplot
    df_plot <- data.frame(
      Expression = as.numeric(gexp),
      Condition = cond # 'Condition' metadata for the subsetted cells
    )
    
    # Ensure Condition is a factor with desired levels for consistent plotting order
    df_plot$Condition <- factor(df_plot$Condition, levels = c("Nr4a1 KO", "WT"))
    
    # Calculate the maximum mean expression + a buffer for y-axis scaling
    # We use mean_cl_normal to get the upper bound of the confidence interval
    summary_data <- df_plot %>%
      group_by(Condition) %>%
      summarise(
        mean_expr = mean(Expression),
        upper_ci = mean_cl_normal(Expression)$ymax # Get upper confidence interval
      )
    
    y_max_current_gene <- max(summary_data$upper_ci, na.rm = TRUE) * 1.2 # Add 20% buffer
    
    # Create the bar plot with mean expression and confidence intervals
    p <- df_plot %>%
      ggplot(aes(x = Condition, y = Expression, fill = Condition)) +
      # Bar for mean expression
      stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), alpha = 0.7, color = "black") +
      # Error bars for 95% confidence interval
      stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
                   position = position_dodge(width = 0.8), width = 0.2, size = 0.5) +
      # Add significance bars
      geom_signif(map_signif_level = TRUE, comparisons = list(c("Nr4a1 KO", "WT")),
                  vjust = 0.5, tip_length = 0.02, # Adjust vjust and tip_length for better appearance
                  test = "wilcox.test", test.args = list(exact = FALSE), # exact = FALSE for larger datasets
                  textsize = 3) + # Adjust text size for p-value
      labs(title = paste0(gene, " Mean Expression in ", ct_sub, " Cells"), # Dynamically include cell type in title
           y = "Mean Expression Level") +
      scale_fill_manual(values = c("Nr4a1 KO" = "#E41A1C", "WT" = "#377EB8")) + # Custom colors for fill
      # scale_color_manual(values = c("Nr4a1 KO" = "#E41A1C", "WT" = "#377EB8")) + # Custom colors for error bars (if needed, though fill handles color)
      # Set y-axis limits dynamically
      coord_cartesian(ylim = c(0, y_max_current_gene)) + # Ensure y-axis starts at 0
      cust_theme +
      NoLegend()
    
    # Add the plot to the list
    plot_list[[gene]] <- p
  } else {
    message(paste0("Gene '", gene, "' not found in the Seurat object. Skipping plot."))
  }
}


# --- Displaying the plots ---
# You can display them one by one:
# print(plot_list[[1]]) # To print the first plot
# print(plot_list$Ascl2) # To print the plot for 'Ascl2'

# Or combine and display them using patchwork (recommended for multiple plots)
# You can arrange them in a grid, e.g., 4 columns:
wrap_plots(plot_list, ncol = 4)

# To save all plots to a PDF file:
# pdf("Stem_Cell_Gene_Mean_Bar_Plots.pdf", width = 10, height = 8 * ceiling(length(plot_list) / 4)) # Adjust width/height for 4 columns
# for (p in plot_list) {
#   print(p)
# }
# dev.off()



