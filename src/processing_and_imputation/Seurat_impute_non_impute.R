# Part 0: Import Packages
#devtools::install_github("KrishnaswamyLab/MAGIC", subdir='Rmagic') # install if needed
library(Rmagic)
#BiocManager::install("MAST") # download if needed # install if needed
library(MAST)
library(Seurat)
library(dplyr)
#BiocManager::install("topGO") # install if needed
library(topGO)
library(sctransform)

# Part 1: Load Integrated Non-Imputed Data and Calculate Clusters
# Load integrated data
# this may need to be changed based on where data gets downloaded and unzipped- too big for github
# may need to setwd() if using RStudio
load('../../../../../../../../Desktop/seurat_integrated_intermediate/integratedData.RData')

# Compute clusters in the non-imputed data
immune.combined = RenameIdents(object = immune.combined, `0` = "T-Cells 1", `1` = "CD8+ T Cells", 
                              `2` = "Myeloid 1", `3` = "Myeloid 2" , `5` = "B-Cells", `6` = "NK Cells", 
                              `7` = "Myeloid 3",`8` = "Plasma Cells",`10`="T-Cells 2")
original_plot = DimPlot(object = immune.combined, reduction = "tsne", do.return =TRUE, label.size = 10)
DefaultAssay(object = immune.combined) = "RNA"

# Compute markers for each cluster in the non-imputed data
all.markers = FindAllMarkers(object=immune.combined,assay='RNA',only.pos=TRUE,test.use="MAST")
all.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)
bc.cell.type.genes = unique(all.markers$gene) # Takes all the unique cell type specific genes

# Make table of number of cells per cluster by sample
original_table = table(Idents(object=immune.combined),immune.combined$sample)
write.table(original_table, '../../data/seurat_tables/original_table')

# Part 2: Do MAGIC on the integrated data object
magic_cells = magic(data = immune.combined,genes=bc.cell.type.genes)
DefaultAssay(object = magic_cells) <- "MAGIC_RNA"

# Scale data and perform PCA and t-SNE on MAGIC-imputed Seurat object
magic_cells = ScaleData(object = magic_cells, verbose = FALSE)
magic_cells = FindVariableFeatures(object = magic_cells, selection.method = "vst", nfeatures = 2000,force.recalc=TRUE)
magic_cells = RunPCA(object = magic_cells, npcs = 30, verbose = FALSE)
magic_cells = RunTSNE(object = magic_cells, reduction = "pca", dims = 1:20)

# Compute clusters in the MAGIC-imputed Seurat object
magic_cells = FindNeighbors(object = magic_cells, reduction = "pca", dims = 1:20)
magic_cells = FindClusters(magic_cells, resolution = 0.2)
imputed_plot = DimPlot(object = magic_cells, reduction = "tsne",assay="MAGIC_RNA", no.axes=TRUE)

# Make table of number of cells per cluster by sample
imputed_table = table(Idents(object=magic_cells),magic_cells$sample)
write.table(original_table, '../../data/seurat_tables/imputed_table')