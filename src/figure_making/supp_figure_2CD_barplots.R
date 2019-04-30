# Part 0: Import Packages
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# NOTE: May need to setwd() if using RStudio, as this is different for each individual machine.

# Part 1: Table with Non-Imputed Counts
original_table = read.table('../../data/seurat_tables/original_table')
original_frame = as.matrix(original_table)
original_frame = prop.table(original_frame, 2)
cols = RColorBrewer::brewer.pal(name = "Paired", n=12)
jpeg("../../output/figures/supp_fig_2C_nonimputed_barplot.jpg", width = 1200, height = 600, pointsize = 15, res = 75)
original_barplot = barplot(original_frame, col = cols, xlim =c(0,11), xaxt = "n", ylab = 'Fraction Cells', cex.axis=1.1, cex.lab=1.1)
original_barplot = legend('right', legend = rownames(original_frame), fill=cols, xpd=TRUE, horiz=FALSE, cex = 1, 
                          x.intersp = 0.15, x = 8.5, y = 1)
original_barplot = axis(side = 1, labels = c('HIV+1,Bld', 'HIV+1,CSF', 'HIV+2,Bld', 'HIV+2,CSF', 'HIV+3,CSF', 'HIV-1,CSF', 'HIV-2,CSF'), at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8))
dev.off()

# Part 2: Table with the Imputed Counts
imputed_table = read.table('../../data/seurat_tables/imputed_table')
imputed_frame = as.matrix(imputed_table)
imputed_frame = prop.table(imputed_frame, 2)
jpeg("../../output/figures/supp_fig_2D_imputed_barplot.jpg", width = 1200, height = 600, pointsize = 15, res = 75)
imputed_barplot = barplot(imputed_frame, col = cols, xlim =c(0,11), xaxt = "n", ylab = 'Fraction Cells', cex.axis=1.1, cex.lab=1.1)
imputed_barplot = legend('right', legend = rownames(original_frame), fill=cols, xpd=TRUE, horiz=FALSE, cex = 1, 
                         x.intersp = 0.15, x = 8.5, y = 1)
imputed_barplot = axis(side = 1, labels = c('HIV+1,Bld', 'HIV+1,CSF', 'HIV+2,Bld', 'HIV+2,CSF', 'HIV+3,CSF', 'HIV-1,CSF', 'HIV-2,CSF'), at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8))
dev.off()