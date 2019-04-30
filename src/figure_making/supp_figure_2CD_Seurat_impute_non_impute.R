# Part 0: Import Packages
#devtools::install_github("KrishnaswamyLab/MAGIC", subdir='Rmagic') # install if needed
library(Rmagic)
#BiocManager::install("MAST") # download if needed # install if needed
library(MAST)
library(Seurat)
library(dplyr)
#BiocManager::install("topGO") # install if needed
library(topGO)
library(ggplot2)
library(sctransform)
library(RColorBrewer)
library(ggpubr)

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

# Figure Making
original_table = table(Idents(object=immune.combined),immune.combined$sample)
original_frame = as.matrix(original_table)
original_frame = prop.table(original_frame, 2)
cols = RColorBrewer::brewer.pal(name = "Paired", n=12)
jpeg("../../output/figures/fig_2_nonimputed_barplot.jpg", width = 1200, height = 600, pointsize = 15, res = 75)
original_barplot = barplot(original_frame, col = cols, xlim =c(0,11), xaxt = "n", ylab = 'Fraction Cells', cex.axis=1.1, cex.lab=1.1)
original_barplot = legend('right', legend = rownames(original_frame), fill=cols, xpd=TRUE, horiz=FALSE, cex = 1, 
       x.intersp = 0.15, x = 8.5, y = 1)
original_barplot = axis(side = 1, labels = c('HIV+1,Bld', 'HIV+1,CSF', 'HIV+2,Bld', 'HIV+2,CSF', 'HIV+3,CSF', 'HIV-1,CSF', 'HIV-2,CSF'), at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8))
original_barplot = title('Non-Imputed Cell Cluster Membership by Sample')
dev.off()

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
imputed_frame = as.matrix(imputed_table)
imputed_frame = prop.table(imputed_frame, 2)
jpeg("../../output/figures/fig_2_imputed_barplot.jpg", width = 1200, height = 600, pointsize = 15, res = 75)
imputed_barplot = barplot(imputed_frame, col = cols, xlim =c(0,11), xaxt = "n", ylab = 'Fraction Cells', cex.axis=1.1, cex.lab=1.1)
imputed_barplot = legend('right', legend = rownames(original_frame), fill=cols, xpd=TRUE, horiz=FALSE, cex = 1, 
       x.intersp = 0.15, x = 8.5, y = 1)
imputed_barplot = axis(side = 1, labels = c('HIV+1,Bld', 'HIV+1,CSF', 'HIV+2,Bld', 'HIV+2,CSF', 'HIV+3,CSF', 'HIV-1,CSF', 'HIV-2,CSF'), at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8))
imputed_barplot = title('Imputed Cell Cluster Membership by Sample')
dev.off()