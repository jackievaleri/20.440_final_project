# Part 0: Import Packages
#devtools::install_github("KrishnaswamyLab/MAGIC", subdir='Rmagic') # download if needed
library(Rmagic)
#BiocManager::install("MAST") # download if needed
library(MAST)

library(Seurat)
library(dplyr)
#BiocManager::install("topGO")
library(topGO)
library(ggplot2)
library(sctransform)


# Part 1: Load Integrated Non-Imputed Data and Calculate Clusters
# Load integrated data
load('integratedDataFinal.RData')

# Compute clusters in the non-imputed data
<<<<<<< HEAD
original_plot = DimPlot(object = immune.combined, reduction = "tsne", split.by = "sample", no.axes=TRUE) # shows integrated plot
=======
DimPlot(object = immune.combined, reduction = "tsne", split.by = "sample") # shows integrated plot
>>>>>>> parent of d03eaa2... updating imputed data
DefaultAssay(object = immune.combined) <- "RNA"

# Compute markers for each cluster in the non-imputed data
all.markers<- FindAllMarkers(object=immune.combined,assay='RNA',only.pos=TRUE,test.use="MAST")
all.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)
<<<<<<< HEAD
bc.cell.type.genes <- unique(all.markers$gene) # Takes all the unique cell type specific genes

# Make table of number of cells per cluster by sample
immune.combined <- RenameIdents(object = immune.combined, `0` = "T-Cells 1", `1` = "CD8+ T Cells", 
                                `2` = "Myeloid 1", `3` = "Myeloid 2" , `5` = "B-Cells", `6` = "NK Cells", 
                                `7` = "Myeloid 3",`8` = "Plasma Cells",`10`="T-Cells 2")
original_table = table(Idents(object=immune.combined),immune.combined$sample)

# Figure Making
original_table = table(Idents(object=immune.combined),immune.combined$sample)
original_frame = as.matrix(original_table)
original_frame = prop.table(original_frame, 2)
cols = c("yellow2", "hotpink4" , "brown", 
         "rosybrown2", "seagreen", "royalblue", "deeppink1", "orange", "red", "white", 
         "darkorchid", "azure3", "cyan", "chartreuse")
barplot(original_frame, col = cols, xlim =c(0,11))
legend('right', legend = rownames(original_frame), fill=cols, xpd=TRUE, horiz=FALSE, cex = 1.3, 
       x.intersp = 0.15, x = 8.5, y = 0.8)
axis(side = 1, labels = c('HIV+1,Bld', 'HIV+1,CSF', 'HIV+2,Bld', 'HIV+2,CSF', 'HIV+3,CSF', 'HIV-1,CSF', 'HIV-2,CSF'), at = c(1, 2.2, 3.4, 4.6, 5.8, 7, 8.2))
=======
top10 <- all.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)
DoHeatmap(object = immune.combined, features = top10$gene)

#number of cells per cluster by sample
table(Idents(object=immune.combined),immune.combined$sample)

bc.cell.type.genes <- unique(all.markers$gene) # Takes all the unique cell type specific genes

immune.combined <- RenameIdents(object = immune.combined, `0` = "T-Cells 1", `1` = "CD8+ T Cells", 
                                `2` = "Myeloid 1", `3` = "Myeloid 2" , `5` = "B-Cells", `6` = "NK Cells", 
                                `7` = "Myeloid 3",`8` = "Plasma Cells",`10`="T-Cells 2")
markersplot<-c("IL7R","LTB","LEF1","CCL5","CD8A","GZMB","HLA-DRA","FCER1A","CD74","S100A9","LZ1","CX3CR1","AIF1","CD14","MS4A1","CD79A","IGHM1","IGJ","IGHG","TRBV7","IL8","FCGR3A")

# Do MAGIC
tcells<-magic(data = immune.combined,genes=bc.cell.type.genes)
>>>>>>> parent of d03eaa2... updating imputed data

DefaultAssay(object = tcells) <- "MAGIC_RNA"
tcells <- ScaleData(object = tcells, verbose = FALSE)
tcells <- FindVariableFeatures(object = tcells, selection.method = "vst", nfeatures = 2000,force.recalc=TRUE)
tcells <- RunPCA(object = tcells, npcs = 30, verbose = FALSE)
tcells <- RunTSNE(object = tcells, reduction = "pca", dims = 1:20)

tcells <- FindNeighbors(object = tcells, reduction = "pca", dims = 1:20)
tcells <- FindClusters(tcells, resolution = 0.2)

<<<<<<< HEAD
# Compute clusters in the MAGIC-imputed Seurat object
magic_cells <- FindNeighbors(object = magic_cells, reduction = "pca", dims = 1:20)
magic_cells <- FindClusters(magic_cells, resolution = 0.2)
imputed_plot = DimPlot(object = magic_cells, reduction = "tsne",assay="MAGIC_RNA", split.by = "sample", no.axes=TRUE)

# Make table of number of cells per cluster by sample
imputed_table = table(Idents(object=magic_cells),magic_cells$sample)
imputed_frame = as.matrix(imputed_table)
imputed_frame = prop.table(imputed_frame, 2)
cols = c("yellow2", "hotpink4" , "brown", 
         "rosybrown2", "seagreen", "royalblue", "deeppink1", "orange", "red", "white", 
                  "darkorchid", "azure3", "cyan", "chartreuse")
barplot(imputed_frame, col = cols)

# Part 3: Output 
=======
DimPlot(object = tcells, reduction = "tsne",assay="MAGIC_RNA", split.by = "sample",label=TRUE)

c0markers<- FindConservedMarkers(tcells,ident.1=0,grouping.var='sample')

tcellMAGICMarkers<- FindAllMarkers(object=tcells,assay='MAGIC_RNA')
>>>>>>> parent of d03eaa2... updating imputed data
