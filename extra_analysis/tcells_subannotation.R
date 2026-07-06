library(dplyr)
library(Seurat)
library(ggplot2)
#library(celda)
library(Matrix)
library(hdf5r)
#library(SCEVAN)
library(dittoSeq)
library(SplineDV)
library(presto)
library(ggpubr)
library(tidyr)
library(stringr)
library(tibble)
library(cowplot)
library(scSGS)
# Annotate T cells
data <- readRDS(file.path(base_dir, "AllSamplesMerged_OldGenesRemoved_StringentFiltered_DoubletsRemoved_AmbientRNARemoved_AllGenesUMAP_FixedKOvsWT_Annotated_SampleSplit.rds"))

data_tcells <- process_and_extract_cell_type(data, "T cells", resolution = 1.0)

Idents(data_tcells) <- "CellType"

DimPlot(object = data_tcells, reduction = "umap", group.by = c("seurat_clusters"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1)

DotPlot(data_tcells, features = c('Cd8a','Cd8b1',
                                  'Cd4', 
                                  'Cd3e', 'Cd3g',
                                  'Ncr1','Klrb1b',
                                  'Batf3','',
                                  'Trdc', 
                                  'Gata3',
                                  'Ms4a2','Kit',
                                  'Tpm2', "Myh11", "Tagln",      
                                  'Pecam1',"Lyve1", "Flt4", 
                                  'Plvap', "Flt1",
                                  "Vip", "Tubb3","Vat1l", 
                                  "Mki67"),
        group.by = 'seurat_clusters', dot.min = 0.05, cols = 'RdBu' ) + # Set color range
  coord_flip()

#> Mast cells <- Kit, Ms4a2, Cd3-
#> Innate Lymphoid cells <- Gata3, Cd3-
#> Gamma-delta T cells <- Trdc, Cd3, Cd4-, Cd8-
#> NK cells <- Ncr1, Klrb1b, Cd3-
#> SMCs <- Tpm2
DimPlot(object = data_tcells, reduction = "umap", group.by = c("seurat_clusters"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1,pt.size = 2)

# Find markers for cluster 9 ONLY
# cluster_markers <- FindMarkers(data_tcells,
#                                 ident.1 = 9, only.pos = TRUE,
#                                 group.by = "seurat_clusters",
#                                 logfc.threshold = 0.5, # Adjust as needed
#                                 min.pct = 0.1) # Adjust as needed
# top10_genes <- rownames(head(cluster_markers[order(cluster_markers$p_val_adj), ], 20))
# print(top10_genes)

data_tcells$new_celltypes = recode_factor(data_tcells$seurat_clusters,
                                          '0' = 'Innate Lymphoid',
                                          '1' = 'T cells',
                                          '2' = 'Cd4+ T cells',
                                          '3' = 'Cd4+ T cells',
                                          '4' = 'NK cells',
                                          '5' = 'Gamma-delta T cells',
                                          '6' = 'Cd8+ T cells',
                                          '7' = 'Cd4+ T cells',
                                          '8' = 'Innate Lymphoid',
                                          '9' = 'Cycling T cells',
                                          '10' = 'NK cells',
                                          '11' = 'SMCs',
                                          '12' = 'Mast cells',
                                          '13' = 'NK cells',
                                          '14' = 'Vascular endothelial',
)

DimPlot(object = data_tcells, reduction = "umap", group.by = 'new_celltypes')

data_tcells$new_celltypes = factor(data_tcells$new_celltypes,
                                   levels = c( 'T cells','Cd4+ T cells',
                                               'Cd8+ T cells','NK cells',
                                               'Gamma-delta T cells',
                                               'Innate Lymphoid',
                                               'Mast cells','SMCs')
)

DotPlot(data_tcells, features = c('Cd3e', 'Cd3g', 'Cd4', 'Cd8a','Cd8b1','Ncr1','Klrb1b','Trdc', 'Trgc1', 'Gata3','Ms4a2','Kit','Tpm2', 'Myh11'),
        group.by = 'new_celltypes', dot.min = 0.05, cols = 'RdBu' ) 


# Import cell types that do not belong to T cells
str_ntc <- c('SMCs', 'NK cells','Innate Lymphoid', 'Mast cells', 'Vascular endothelial')
idx <- data_tcells$new_celltypes %in% str_ntc
barcodes_ntc <- colnames(data_tcells)[idx]  # Works only if colnames(data_tcells) are barcodes
idx2 <- which(colnames(data) %in% barcodes_ntc)
data$CellType[idx2] <- as.character(data_tcells$new_celltypes[idx])


DimPlot(object = data, reduction = "umap", group.by = c("CellType"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1)
DimPlot(object = data, reduction = "umap", split.by = c("CellType"), group.by = c("CellType"),
        alpha = 1, ncol = 4)

# Create sub cell type annotation layer
data$SubCellType <- data$CellType  # Preserve broad annotation
common_cells_tcell <- intersect(rownames(data@meta.data), rownames(data_tcells@meta.data))
new_cell_types_subset <- as.character(data_tcells@meta.data[common_cells_tcell, "new_celltypes"])
data@meta.data[common_cells_tcell, "SubCellType"] <- new_cell_types_subset

DimPlot(object = data, reduction = "umap", group.by = c("SubCellType"),
        repel = TRUE, label.size = 3, alpha = 1)

saveRDS(data, file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated.rds"))