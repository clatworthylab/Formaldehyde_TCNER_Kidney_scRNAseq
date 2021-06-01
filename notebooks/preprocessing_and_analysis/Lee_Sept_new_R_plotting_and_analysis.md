# @Author: kt16
# @Date:   2019-09-30 20:24:14
# @Last Modified by:   kt16
# @Last Modified time: 2021-02-15 20:40:15
# @Notes:	making vectorized figures for Lee's paper

### DEG dot plot
```R
# and now to tabulate in R
library(dplyr)
library(readr)
library(tibble)
library(pbmcapply)

setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG')
files <- as.list(list.files(pattern = '.txt'))
names(files) <- gsub('.txt', '', list.files(pattern = '.txt'))
degs <- lapply(files, read_tsv)

# lfc threshold of 1
res <- pbmclapply(degs, function(x){
	Adh5_up <- x %>% filter(Adh5_pvals_adj < 0.05 & Adh5_logfoldchanges >= 1) %>% select(X1) %>% unlist %>% length
	Csb_up <- x %>% filter(Csb_pvals_adj < 0.05 & Csb_logfoldchanges >= 1) %>% select(X1) %>% unlist %>% length
	DKO_up <- x %>% filter(DKO_pvals_adj < 0.05 & DKO_logfoldchanges >= 1) %>% select(X1) %>% unlist %>% length
	Adh5_dn <- x %>% filter(Adh5_pvals_adj < 0.05 & Adh5_logfoldchanges <= -1) %>% select(X1) %>% unlist %>% length
	Csb_dn <- x %>% filter(Csb_pvals_adj < 0.05 & Csb_logfoldchanges <= -1) %>% select(X1) %>% unlist %>% length
	DKO_dn <- x %>% filter(DKO_pvals_adj < 0.05 & DKO_logfoldchanges <= -1) %>% select(X1) %>% unlist %>% length
	y <- data.frame(Adh5_up = Adh5_up, Csb_up = Csb_up, DKO_up = DKO_up, Adh5_dn = Adh5_dn, Csb_dn = Csb_dn, DKO_dn = DKO_dn)
	return(y)
}, mc.cores = 24)

results <- do.call(rbind, res)
results

# Lee wants these as a dot plot
results_up <- results[,1:3]
results_down <- results[,4:6]

results_up$id = rownames(results)
results_down$id = rownames(results)

library(reshape2)
results_up = melt(results_up)
results_down = melt(results_down)

colnames(results_up) <- c('Cell type', 'Group', 'No. of DEGs')
colnames(results_down) <- c('Cell type', 'Group', 'No. of DEGs')

results <- rbind(results_up, results_down)

g <- ggplot(results, aes(x=Group, y=`Cell type`, size=`No. of DEGs`, fill=`No. of DEGs`)) + 
	geom_point(shape = 21, color = "grey") + 
	scale_fill_gradientn(colors=c('#fff5f0', '#fee0d2', '#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d')) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), legend.key=element_blank()) + 
	scale_size(range=c(0,4))
ggsave("../../figures/dotplot/DEGs1.pdf", g, h = 6, w = 4, useDingbats = FALSE)

g <- ggplot(results, aes(x=Group, y=`Cell type`, size=`No. of DEGs`, fill=`No. of DEGs`)) + 
	geom_point(shape = 21, color = "grey") + 
	scale_fill_gradientn(colors=c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), legend.key=element_blank()) + 
	scale_size(range=c(0,4))
ggsave("../../figures/dotplot/DEGs2.pdf", g, h = 6, w = 4, useDingbats = FALSE)
# results_up$Group <- gsub("_up", "", results_up$Group)
# results_down$Group <- gsub("_dn", "", results_down$Group)

library(ggplot2)
g <- ggplot(results_up, aes(x=Group, y=`Cell type`, size=`No. of DEGs`, fill=`No. of DEGs`))+geom_point(shape = 21, color = "grey") + scale_fill_gradientn(colors=c('#fff5f0', '#fee0d2', '#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), legend.key=element_blank()) + scale_size(range=c(0,6))
ggsave("../../figures/dotplot/DEGsUp.pdf", g, h = 6, w = 4, useDingbats = FALSE)

g <- ggplot(results_down, aes(x=Group, y=`Cell type`, size=`No. of DEGs`, fill=`No. of DEGs`))+geom_point(shape = 21, color = "grey") + scale_fill_gradientn(colors=c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), legend.key=element_blank()) + scale_size(range=c(0,6))
ggsave("../../figures/dotplot/DEGsDown.pdf", g, h = 6, w = 4, useDingbats = FALSE)

```

### umaps
```R
library(Seurat) # version 3.1
library(dplyr)
library(ggplot2)
setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/')

# load in file
lee <- readRDS(file = "./scanpy/out/Lee_Sept_new.rds")

c.palette <- c('#FABD2E', '#F58700', '#296B7E', '#1C5389', '#1D3061', '#113A5B', '#80929E', '#FC5C03', '#53A592', '#B92A1F', '#3B4854', '#EC3624', '#3A91D0', '#3A91D0', '#3A91D0', '#7ADFC4', '#FB4649', '#FF837F')

g.palette <- c('#3a4854', '#51a693', '#1b538a', '#ef3723')

small_legend <- function(...){
	small_legend_theme <- 
	theme(legend.title = element_text(size = 4), 
		legend.text  = element_text(size = 4),
		legend.key.size = unit(0.01, "lines"), ...)
	return(small_legend_theme)
}

small_guide <- function(...){
	small_guide <- 
	guides(shape = guide_legend(override.aes = list(size = 0.1)), 
		color = guide_legend(override.aes = list(size = 0.1)), ...)
	return(small_guide)
}

topright_legend <- function(...){
	legend <- theme(legend.position = c(.99, .99), legend.justification = c('right', 'top'), legend.box.just = "right", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

bottomleft_legend <- function(...){
	legend <- theme(legend.position = c(.01, .01), legend.justification = c('left', 'bottom'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

bottomright_legend <- function(...){
	legend <- theme(legend.position = c(.99, .01), legend.justification = c('right', 'bottom'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}


small_axis <- function(...){
	axis <- theme(text = element_text(size=8), axis.text = element_text(size=8), axis.line = element_line(size = 0.1), axis.ticks = element_line(size = 0.1), ...)
	return(axis)
}

lee$group = factor(lee$group, levels = c('WT', 'Adh5', 'Csb', 'DKO'))

pdf('./scanpy/figures/umap/Lee_Sept_new_labels.pdf', h = 5, w = 5)
DimPlot(lee, group.by = "celltype", cols = c.palette) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

pdf('./scanpy/figures/umap/Lee_Sept_new_labels_group.pdf', h = 5, width = 20)
DimPlot(lee, group.by = "celltype", cols = c.palette, split.by = "group", ncol = 4) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

WT <- subset(lee, subset = group == 'WT')
Adh5 <- subset(lee, subset = group == 'Adh5')
Csb <- subset(lee, subset = group == 'Csb')
DKO <- subset(lee, subset = group == 'DKO')

pdf('./scanpy/figures/umap/Lee_Sept_new_labels_WT.pdf', h = 5, w = 5)
DimPlot(WT, group.by = "celltype", cols = c.palette, label = F) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

pdf('./scanpy/figures/umap/Lee_Sept_new_labels_Adh5.pdf', h = 5, w = 5)
DimPlot(Adh5, group.by = "celltype", cols = c.palette, label = F) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

pdf('./scanpy/figures/umap/Lee_Sept_new_labels_Csb.pdf', h = 5, w = 5)
DimPlot(Csb, group.by = "celltype", cols = c.palette, label = F) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

pdf('./scanpy/figures/umap/Lee_Sept_new_labels_DKO.pdf', h = 5, w = 5)
DimPlot(DKO, group.by = "celltype", cols = c.palette, label = F) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

pdf('./scanpy/figures/umap/Lee_Sept_new_group.pdf', h = 5, width = 20)
DimPlot(lee, group.by = "group", cols = g.palette, label = F, split.by = "group", ncol = 4) + NoAxes() + NoLegend()
dev.off()

pdf('./scanpy/figures/umap/Lee_Sept_new_Cdkn1a.pdf', h = 5, width = 20)
FeaturePlot(lee, "Cdkn1a", cols = c('grey', viridis::plasma(49)), split.by = "group", ncol = 4, order = TRUE, pt.size = 2) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

pdf('./scanpy/figures/umap/Lee_Sept_new_Gdf15.pdf', h = 5, width = 20)
FeaturePlot(lee, "Gdf15", cols = c('grey', viridis::plasma(49)), split.by = "group", ncol = 4, order = TRUE, pt.size = 2) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

pdf('./scanpy/figures/umap/Lee_Sept_new_Havcr1.pdf', h = 5, width = 20)
FeaturePlot(lee, "Havcr1", cols = c('grey', viridis::plasma(49)), split.by = "group", ncol = 4, order = TRUE, pt.size = 2) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

# cell numbers
table(lee$celltype, lee$group)
  #                  WT  Adh5   Csb   DKO
  # B               277   226   400   500
  # Basophil         14     4     8     6
  # CD-IC           871   242   699   855
  # CD-PC            30    12    36     2
  # CT             1335   461   989  1166
  # DCT             128    47    97    55
  # Endothelial    1341   528  1262   574
  # Erythrocyte     823   118   405   151
  # LOH            1305   337   978   974
  # MNP            1126   565   901  2184
  # Myofibroblast   368   177   459   328
  # Neutrophil       37    15    19   125
  # PT            23272 12899 22090 19673
  # PT-DKO          247   144   299  2391
  # PT-Novel1        21    25   131   151
  # Podocyte         17    10    33    29
  # T_NK            402   207   271   626
  # pDC              10     6     5    12

Idents(lee) <- 'group'
library(patchwork)

pdf('./scanpy/figures/violinplot/Lee_Sept_new_Cdkn1a_vln.pdf', h = 3.5, width = 6)
p1 <- VlnPlot(lee, features = "Cdkn1a", cols = c.palette, pt.size = 0, idents = 'WT', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('WT')
p2 <- VlnPlot(lee, features = "Cdkn1a", cols = c.palette, pt.size = 0, idents = 'Adh5', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('Adh5')
p3 <- VlnPlot(lee, features = "Cdkn1a", cols = c.palette, pt.size = 0, idents = 'Csb', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('Csb')
p4 <- VlnPlot(lee, features = "Cdkn1a", cols = c.palette, pt.size = 0, idents = 'DKO', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('DKO')
p1$layers[[1]]$aes_params$size = 0.1
p2$layers[[1]]$aes_params$size = 0.1
p3$layers[[1]]$aes_params$size = 0.1
p4$layers[[1]]$aes_params$size = 0.1
(p1 + p2) / (p3 + p4)
dev.off()

pdf('./scanpy/figures/violinplot/Lee_Sept_new_Havcr1_vln.pdf', h = 3.5, width = 6)
p1 <- VlnPlot(lee, features = "Havcr1", cols = c.palette, pt.size = 0, idents = 'WT', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('WT')
p2 <- VlnPlot(lee, features = "Havcr1", cols = c.palette, pt.size = 0, idents = 'Adh5', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('Adh5')
p3 <- VlnPlot(lee, features = "Havcr1", cols = c.palette, pt.size = 0, idents = 'Csb', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('Csb')
p4 <- VlnPlot(lee, features = "Havcr1", cols = c.palette, pt.size = 0, idents = 'DKO', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('DKO')
p1$layers[[1]]$aes_params$size = 0.1
p2$layers[[1]]$aes_params$size = 0.1
p3$layers[[1]]$aes_params$size = 0.1
p4$layers[[1]]$aes_params$size = 0.1
(p1 + p2) / (p3 + p4)
dev.off()

pdf('./scanpy/figures/violinplot/Lee_Sept_new_Gdf15_vln.pdf', h = 3.5, width = 6)
p1 <- VlnPlot(lee, features = "Gdf15", cols = c.palette, pt.size = 0, idents = 'WT', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('WT')
p2 <- VlnPlot(lee, features = "Gdf15", cols = c.palette, pt.size = 0, idents = 'Adh5', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('Adh5')
p3 <- VlnPlot(lee, features = "Gdf15", cols = c.palette, pt.size = 0, idents = 'Csb', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('Csb')
p4 <- VlnPlot(lee, features = "Gdf15", cols = c.palette, pt.size = 0, idents = 'DKO', group.by = 'celltype', lwd = 0.1) + NoLegend() + small_axis() + ggtitle('DKO')
p1$layers[[1]]$aes_params$size = 0.1
p2$layers[[1]]$aes_params$size = 0.1
p3$layers[[1]]$aes_params$size = 0.1
p4$layers[[1]]$aes_params$size = 0.1
(p1 + p2) / (p3 + p4)
dev.off()

Idents(lee) <- 'celltype'
pdf('./scanpy/figures/violinplot/Lee_Sept_new_Cdkn1a_vln_merged.pdf', h = 2.5, width = 5)
p1 <- VlnPlot(lee, features = "Cdkn1a", cols = c.palette, pt.size = 0.1, lwd = 0.1) + NoLegend() + small_axis()
p1$layers[[1]]$aes_params$size = 0.1
p1
dev.off()

pdf('./scanpy/figures/violinplot/Lee_Sept_new_Havcr1_vln_merged.pdf', h = 2.5, width = 5)
p1 <- VlnPlot(lee, features = "Havcr1", cols = c.palette, pt.size = 0.1, lwd = 0.1) + NoLegend() + small_axis()
p1$layers[[1]]$aes_params$size = 0.1
p1
dev.off()
pdf('./scanpy/figures/violinplot/Lee_Sept_new_Gdf15_vln_merged.pdf', h = 2.5, width = 5)
p1 <- VlnPlot(lee, features = "Gdf15", cols = c.palette, pt.size = 0.1, lwd = 0.1) + NoLegend() + small_axis()
p1$layers[[1]]$aes_params$size = 0.1
p1
dev.off()

```

### umap for pt
```R
library(Seurat) # version 3.1
library(dplyr)
library(ggplot2)
setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/')

# load in file
pt <- readRDS(file = "./scanpy/out/Lee_Sept_new_pt.rds")

c.palette <- c('#c7e9b4','#7fcdbb','#41b6c4','#41b6c4','#ef3723','#1d91c0', '#225ea8', '#0c2c84')
g.palette <- c('#3a4854', '#51a693', '#1b538a', '#ef3723')

small_legend <- function(...){
	small_legend_theme <- 
	theme(legend.title = element_text(size = 4), 
		legend.text  = element_text(size = 4),
		legend.key.size = unit(0.01, "lines"), ...)
	return(small_legend_theme)
}

small_guide <- function(...){
	small_guide <- 
	guides(shape = guide_legend(override.aes = list(size = 0.1)), 
		color = guide_legend(override.aes = list(size = 0.1)), ...)
	return(small_guide)
}

topright_legend <- function(...){
	legend <- theme(legend.position = c(.99, .99), legend.justification = c('right', 'top'), legend.box.just = "right", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

bottomleft_legend <- function(...){
	legend <- theme(legend.position = c(.01, .01), legend.justification = c('left', 'bottom'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

bottomright_legend <- function(...){
	legend <- theme(legend.position = c(.99, .01), legend.justification = c('right', 'bottom'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}


small_axis <- function(...){
	axis <- theme(text = element_text(size=8), axis.text = element_text(size=8), axis.line = element_line(size = 0.1), axis.ticks = element_line(size = 0.1), ...)
	return(axis)
}

small_axis <- function(...){
	axis <- theme(text = element_text(size=8), axis.text = element_text(size=8), axis.line = element_line(size = 0.1), axis.ticks = element_line(size = 0.1), ...)
	return(axis)
}

```

### combine pt annotation for plotting
```R
library(Seurat) # version 3.1
library(dplyr)
library(ggplot2)
setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/')

# load in file
lee <- readRDS(file = "./scanpy/out/Lee_Sept_new.rds")
pt <- readRDS(file = "./scanpy/out/Lee_Sept_new_pt.rds")

pt$pt_celltype <- paste0('PT-',pt$leiden_pt)
pt$pt_celltype <- gsub('PT-5', 'PT-5-DKO', pt$pt_celltype)


lee$pt_celltype <- as.character(lee$celltype)

id <- match(row.names(pt@meta.data), row.names(lee@meta.data))

lee$pt_celltype[id] <- pt$pt_celltype
lee$pt_celltype <- factor(lee$pt_celltype, levels = c('B', 'Basophil', 'CD-IC', 'CD-PC', 'CT', 'DCT', 'Endothelial', 'Erythrocyte', 'LOH', 'MNP', 'Myofibroblast', 'Neutrophil', 'Novel1', 'PT-0', 'PT-1', 'PT-2', 'PT-3', 'PT-4', 'PT-5-DKO', 'PT-6', 'PT-7', 'Podocyte', 'T_NK', 'pDC'))
pt.palette <- c('#D9D9D9','#D1D1D1','#CACACA','#C2C2C2','#BABABA','#B0B0B0','#A5A5A5','#9B9B9B','#919191','#888888','#7E7E7E','#757575','#6C6C6C',
	'#023fa5','#4a6fe3','#7d87b9','#8595e1','#b5bbe3', '#ef3723', '#bec1d4', '#051094', '#636363','#5A5A5A','#525252')
g.palette <- c('#3a4854', '#51a693', '#1b538a', '#ef3723')

small_legend <- function(...){
	small_legend_theme <- 
	theme(legend.title = element_text(size = 4), 
		legend.text  = element_text(size = 4),
		legend.key.size = unit(0.01, "lines"), ...)
	return(small_legend_theme)
}

small_guide <- function(...){
	small_guide <- 
	guides(shape = guide_legend(override.aes = list(size = 0.1)), 
		color = guide_legend(override.aes = list(size = 0.1)), ...)
	return(small_guide)
}

topright_legend <- function(...){
	legend <- theme(legend.position = c(.99, .99), legend.justification = c('right', 'top'), legend.box.just = "right", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

bottomleft_legend <- function(...){
	legend <- theme(legend.position = c(.01, .01), legend.justification = c('left', 'bottom'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

bottomright_legend <- function(...){
	legend <- theme(legend.position = c(.99, .01), legend.justification = c('right', 'bottom'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}


small_axis <- function(...){
	axis <- theme(text = element_text(size=8), axis.text = element_text(size=8), axis.line = element_line(size = 0.1), axis.ticks = element_line(size = 0.1), ...)
	return(axis)
}

lee$group = factor(lee$group, levels = c('WT', 'Adh5', 'Csb', 'DKO'))

pdf('./scanpy/figures/umap/Lee_Sept_new_detailed_labels_group.pdf', h = 5, width = 20)
DimPlot(lee, group.by = "pt_celltype", cols = pt.palette, split.by = "group", ncol = 4) + NoAxes() + small_legend() + bottomright_legend()
dev.off()

```

### bubblechart for pt
```R
library(dplyr)
library(readr)
# setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset')
setwd('~/Documents/Clatworthy_scRNAseq/Lee_analysis/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset')
files <- as.list(paste0(c(0:7), '.txt'))

degs <- lapply(files, read_tsv)
names(degs) <- gsub(".txt", "", files)

## lets start with DKO?
# extract the significant genes
gene_df <- lapply(degs, function(x){
	df <- x %>% dplyr::filter(DKO_pvals_adj < 0.01) %>% dplyr::select(X1, Adh5_logfoldchanges, Csb_logfoldchanges, DKO_logfoldchanges, Adh5_pvals_adj, Csb_pvals_adj, DKO_pvals_adj)
	colnames(df) <- c("genes", 'Adh5_logfoldchanges', 'Csb_logfoldchanges', 'DKO_logfoldchanges', 'Adh5_pvals_adj', 'Csb_pvals_adj', 'DKO_pvals_adj')
	df$ranking <- df$DKO_logfoldchanges
	df <- df[order(df$ranking), ]
	df$order <- order(df$ranking)
	top <- df %>% top_n(25) %>% dplyr::select(genes)
	bottom <- df %>% top_n(-25) %>% dplyr::select(genes)
	keep <- unlist(c(top, bottom))
	df$top_genes <- df$genes
	df$top_genes[!df$top_genes %in% keep] <- NA
	return(df)
})

gene_df <- lapply(gene_df, function(x){
	x <- x %>% tail(25)
	return(x)
})

bubblefall <- function(df){
	require(ggplot2)
	require(ggrepel)
	p <- ggplot() + geom_hline(aes(yintercept=0), colour = 'grey', linetype = 'dashed', size = 0.5) + labs(x = NULL) +
	labs(y = expression('Log'[2]*'Fold-change')) +
	theme_bw() +
	# guides(size=FALSE) +
	theme(axis.line = element_line(colour = 'black'),
     		legend.position = 'bottom', 
        	panel.grid.major = element_blank(),
        	panel.grid.minor = element_blank(),
        	panel.border = element_blank(),
        	panel.background = element_blank(),
        	axis.line.x=element_blank(),
        	axis.text.x=element_blank(),
        	axis.ticks.x=element_blank()) +     
     theme(plot.title = element_text(lineheight = 0.8, hjust = 0.5))

	p <- p + geom_point(data = df, aes(x = order, y = Adh5_logfoldchanges, size = -log10(Adh5_pvals_adj), fill = 'Adh5', colour = 'Adh5'), shape = 21, alpha = 0.7, stroke = 0.1) + 
		 geom_point(data = df, aes(x = order, y = Csb_logfoldchanges, size = -log10(Csb_pvals_adj), fill = 'Csb', colour = 'Csb'), shape = 21, alpha = 0.7, stroke = 0.1) +
		 geom_point(data = df, aes(x = order, y = DKO_logfoldchanges, size = -log10(DKO_pvals_adj), fill = 'DKO', colour = 'DKO'), shape = 21, alpha = 0.7, stroke = 0.1)
	p <- p + scale_size(range = c(1,4))

    p <- p + geom_text_repel(aes(x = df$order, y = df$DKO_logfoldchanges, label = df$top_genes, size=3), box.padding = unit(0.3, 'lines'), segment.size = 0.2, force=0.5)

	p <- p + scale_fill_manual(name = '', values = c('Adh5' = '#51a693', 'Csb' = '#1b538a', 'DKO' = '#ef3723'), breaks=c('Adh5', 'Csb', 'DKO'), guide = 'legend') +
         scale_colour_manual(name = '', values = c('Adh5' = '#FFFFFF', 'Csb' = '#FFFFFF', 'DKO' = '#FFFFFF'), breaks=c('Adh5', 'Csb', 'DKO'), guide = 'legend')

	# p <- p + theme(legend.title = element_blank())
	return(p)}

small_legend <- function(...){
	small_legend_theme <- 
	theme(legend.title = element_text(size = 4), 
		legend.text  = element_text(size = 4),
		legend.key.size = unit(0.01, "lines"), ...)
	return(small_legend_theme)
}

small_guide <- function(...){
	small_guide <- 
	guides(shape = guide_legend(override.aes = list(size = 0.1)), 
		color = guide_legend(override.aes = list(size = 0.1)), ...)
	return(small_guide)
}


small_axis <- function(...){
	axis <- theme(text = element_text(size=8), axis.text = element_text(size=8), axis.line = element_line(size = 0.1), axis.ticks = element_line(size = 0.1), ...)
	return(axis)
}

topleft_legend <- function(...){
	legend <- theme(legend.position = c(.01, .99), legend.justification = c('left', 'top'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

bottomleft_legend <- function(...){
	legend <- theme(legend.position = c(.01, .01), legend.justification = c('left', 'bottom'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

pdf("../../../figures/dotplot/PT-0_topgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['0']]) +
	ggtitle('PT-0') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-1_topgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['1']]) +
	ggtitle('PT-1') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-2_topgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['2']]) +
	ggtitle('PT-2') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-3_topgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['3']]) +
	ggtitle('PT-3') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-4-DKO_topgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['4']]) + scale_size_continuous(range = c(1, 4), breaks = c(-log10(5e-2),-log10(1e-2),-log10(1e-3),-log10(1e-4),-log10(1e-5),-log10(1e-10), -log10(1e-20))) +
	ggtitle('PT-4-DKO') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-5_topgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['5']]) +
	ggtitle('PT-5') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-6_topgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['6']]) +
	ggtitle('PT-6') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	bottomleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-7-Novel1_topgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['7']]) +
	ggtitle('PT-7-Novel1') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

```

### Lee found some cisplastin dataset and asked if i can analyse it
### I found a NFE2L2 gene set from nature biotech paper
```R
# The resulting DEG file was transferred to the hub
library(Seurat) # version 3.1
library(dplyr)
library(ggplot2)
setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/')

# load in file
lee <- readRDS(file = "./scanpy/out/Lee_Sept_new.rds")
pt <- readRDS(file = "./scanpy/out/Lee_Sept_new_pt.rds")

pt$pt_celltype <- paste0('PT-',pt$leiden_pt)
pt$pt_celltype <- gsub('PT-5', 'PT-5-DKO', pt$pt_celltype)

lee$pt_celltype <- as.character(lee$celltype)

id <- match(row.names(pt@meta.data), row.names(lee@meta.data))

lee$pt_celltype[id] <- pt$pt_celltype
lee$pt_celltype <- factor(lee$pt_celltype, levels = c('B', 'Basophil', 'CD-IC', 'CD-PC', 'CT', 'DCT', 'Endothelial', 'Erythrocyte', 'LOH', 'MNP', 'Myofibroblast', 'Neutrophil', 'Novel1', 'PT-0', 'PT-1', 'PT-2', 'PT-3', 'PT-4', 'PT-5-DKO', 'PT-6', 'PT-7', 'Podocyte', 'T_NK', 'pDC'))
c.palette <- c('#D9D9D9','#D1D1D1','#CACACA','#C2C2C2','#BABABA','#B0B0B0','#A5A5A5','#9B9B9B','#919191','#888888','#7E7E7E','#757575','#6C6C6C',
	'#006BA4', '#ef3723', '#636363','#5A5A5A','#525252')
pt.palette <- c('#D9D9D9','#D1D1D1','#CACACA','#C2C2C2','#BABABA','#B0B0B0','#A5A5A5','#9B9B9B','#919191','#888888','#7E7E7E','#757575','#6C6C6C',
	'#023fa5','#4a6fe3','#7d87b9','#8595e1','#b5bbe3', '#ef3723', '#bec1d4', '#051094', '#636363','#5A5A5A','#525252')
g.palette <- c('#3a4854', '#51a693', '#1b538a', '#ef3723')

lee$group = factor(lee$group, levels = c('WT', 'Adh5', 'Csb', 'DKO'))
pt$group = factor(pt$group, levels = c('WT', 'Adh5', 'Csb', 'DKO'))

g.palette <- c('#3a4854', '#51a693', '#1b538a', '#ef3723')

# read in the deg list
gs1 <- kelvinny::parse_gmt('./dataset/NFE2L2.V2.gmt')
gs1 <- as.vector(gs1[-1,])
gs2 <- kelvinny::parse_gmt('./dataset/SENESCENCE_UP.gmt')
gs2 <- as.vector(gs2[-1,])

# # Eoin's Kidney senescence genes
# gs3 <- c('Mdk', 'Spp1', 'Lama3', 'Pros1', 'Nrp1', 'Egfr', 'Jag1', 'Ephb2', 'Pdia3', 'Mal2', 'Lgals3')

# also enrich for a cisplatin paper
# The preprocessing was performed in GSE117167.md
cis <- readr::read_csv('./dataset/cisplastinvsctrl_deg.csv')

# filter for the top 200 up-regulated deg
cis.f <- cis %>% filter(log2FoldChange > 0 & padj < 0.05)
cis.f <- cis.f[order(cis.f$log2FoldChange, decreasing = TRUE), ]

# convert to mouse symbols
library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
m <- getBM(attributes=c("external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=mart)

# genes <- list('NFE2L2 activation' = m$external_gene_name[m$hsapiens_homolog_associated_gene_name %in% gs1], 'Senescence UP' = m$external_gene_name[m$hsapiens_homolog_associated_gene_name %in% gs2], 'Cisplatin UP' = cis.f$X1[1:20])

genes <- list('FRIDMAN - Senescence UP' = m$external_gene_name[m$hsapiens_homolog_associated_gene_name %in% gs2], 'KIM - NFE2L2 activation' = m$external_gene_name[m$hsapiens_homolog_associated_gene_name %in% gs1], 'YIMIT - Cisplatin UP' = cis.f$X1[1:20])

# run geneset test
enrichGS <- function(seu, geneset){
	require(Seurat) # version 3.1

	seu.gs <- AddModuleScore(seu, geneset, name = names(geneset))
	return(seu.gs@meta.data)
}

objects <- list(lee, pt)
enrich <- lapply(objects, enrichGS, genes)

plotMeanEnrichment <- function(scdata, geneset, enrichment, labels, subset = NULL, annotateColumn = TRUE, heat_cols = c(rep("black",4), viridis::inferno(6)), ...)
{
	require(dplyr)
	require(pbmcapply)
	cat(crayon::red("Formating table and calculating mean"), sep = "\n")

	if (class(scdata) %in% c("SingleCellExperiment", "SummarizedExperiment")) {
        cat("data provided is a SingleCellExperiment/SummarizedExperiment object", sep = "\n")
        cat("extracting colData", sep = "\n")
        require(SummarizedExperiment)
        require(SingleCellExperiment)
        metadata <- colData(scdata)
    } else if (class(scdata) == "Seurat") {
        cat("data provided is a Seurat object", sep = "\n")
        cat("extracting metadata", sep = "\n")
        metadata <- scdata@meta.data
    }

	if(length(labels) > 1){
		dat <- split(enrichment, enrichment[[labels[1]]])
		dat <- lapply(dat, function(x){
			y <- split(x, x[[labels[2]]])
			return(y)})
		dat <- pbmclapply(dat, function(z) {
			y <- lapply(z, function(x) {
				x <- x[ ,-c(1:ncol(metadata))]
				x <- colMeans(x)
				return(x)}
				)
		}, mc.cores = parallel::detectCores())
	} else {
		dat <- split(enrichment, enrichment[[labels]])
		dat <- pbmclapply(dat, function(x) {
			x <- x[ ,-c(1:ncol(metadata))]
			x <- colMeans(x)
			return(x)
		}, mc.cores = parallel::detectCores())
	}

	if(!is.null(subset)){
		cat(crayon::blue(paste0("subsetting to ", subset)), sep = "\n")
		dat <- dat[names(dat) %in% subset]
	}

	cat(crayon::green("Converting to matrix"), sep = "\n")
	if(length(labels) > 1){
		dat <- pbmclapply(dat, function(x){
			y <- do.call(rbind, x)
			return(y)
		}, mc.cores = parallel::detectCores())
		for(i in 1:length(dat)){
			rownames(dat[[i]]) <- paste0(names(dat)[i], "_", rownames(dat[[i]]))
		}
		dat <- do.call(rbind, dat)
		cat(crayon::yellow("Converting any NA values to the smallest value in each enrichment"), sep = "\n")
		dat <- apply(dat, 2, function(x){
			x[which(is.na(x))] <- min(x, na.rm = TRUE)
			return(x)
		})
		mat <- as.matrix(dat)
	} else {
		dat <- do.call(rbind, dat)
		cat(crayon::yellow("Converting any NA values to the smallest value in each enrichment"), sep = "\n")
		dat <- apply(dat, 2, function(x){
			x[which(is.na(x))] <- min(x, na.rm = TRUE)
			return(x)
		})
		mat <- as.matrix(dat)
	}

	# remove the last number from the gene set name
	colnames(mat) <- substr(colnames(mat), 1,nchar(colnames(mat))-1)

	if (annotateColumn){
		if(length(labels) > 1){
			annotation_col <- data.frame(row.names = row.names(mat), label1 = as.factor(
				gsub(
					paste(
						unique(
							paste0("_", gsub(".*_", "", row.names(mat)), "$"
								)
							), collapse = "|"), "", row.names(mat)
					)
				), 
			label2 = as.factor(gsub(".*_", "",row.names(mat))))
		} else {
			annotation_col <- data.frame(row.names = row.names(mat), label = row.names(mat))
		}
	}
	colnames(annotation_col) <- labels

	cat(crayon::cyan("Plotting"), sep = "\n")
	if (annotateColumn){
	kelvinny::plotHeat(t(mat),
			cellheight = 10,
			cellwidth = 10,
			treeheight_col = 10,
			treeheight_row = 10,
			color = heat_cols,
			annotation_col = annotation_col,
			...)
	} else {
	kelvinny::plotHeat(t(mat),
			cellheight = 10,
			cellwidth = 10,
			treeheight_col = 10,
			treeheight_row = 10,
			color = heat_cols,
			...)
	}
}

c.palette1 <- c.palette
names(c.palette1) <- levels(lee$celltype)
c.palette2 <- pt.palette
names(c.palette2) <- levels(lee$pt_celltype)
c.palette1 <- c.palette1[14:15]
c.palette2 <- c.palette2[14:21]

names(g.palette) <- c("WT", "Adh5", "Csb", "DKO")
annotation_color1 <- list(celltype = c.palette1, group = g.palette)
annotation_color2 <- list(pt_celltype = c.palette2, group = g.palette)

pdf('./scanpy/figures/heatmap/geneset_enrichment_lee.pdf', h = 10, w = 10)
plotMeanEnrichment(lee, genes, enrich[[1]], subset = c("PT", "PT-DKO"), labels = c("celltype", "group"), annotation_color = annotation_color1, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

pdf('./scanpy/figures/heatmap/geneset_enrichment_pt.pdf', h = 10, w = 10)
plotMeanEnrichment(lee, genes, enrich[[1]], subset = c('PT-0', 'PT-1', 'PT-2', 'PT-3', 'PT-4', 'PT-5-DKO', 'PT-6', 'PT-7'), labels = c("pt_celltype", "group"), annotation_color = annotation_color2, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

```

### check gene length
```R
setwd("/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG")
library(dplyr)
library(readr)
library(pbmcapply)
files <- as.list(list.files(pattern='.txt'))
degs <- lapply(files, read_tsv)

library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
m <- getBM(attributes=c("external_gene_name", "start_position", "end_position"), mart=mart)
m$geneLength <- m$end_position - m$start_position

library(dplyr)

lengths <- pbmclapply(degs, function(x) {
	Adh5_up <- x %>% dplyr::filter(Adh5_pvals_adj < 0.05 & Adh5_logfoldchanges > 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	Csb_up <- x %>% dplyr::filter(Csb_pvals_adj < 0.05 & Csb_logfoldchanges > 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	DKO_up <- x %>% dplyr::filter(DKO_pvals_adj < 0.05 & DKO_logfoldchanges > 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	Adh5_dn <- x %>% dplyr::filter(Adh5_pvals_adj < 0.05 & Adh5_logfoldchanges < 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	Csb_dn <- x %>% dplyr::filter(Csb_pvals_adj < 0.05 & Csb_logfoldchanges < 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	DKO_dn <- x %>% dplyr::filter(DKO_pvals_adj < 0.05 & DKO_logfoldchanges < 0) %>% dplyr::select(X1) %>% unlist %>% as.character

	Adh5_up_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% Adh5_up)]), error = function(e) return(NA))
	Csb_up_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% Csb_up)]), error = function(e) return(NA))
	DKO_up_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% DKO_up)]), error = function(e) return(NA))
	Adh5_dn_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% Adh5_dn)]), error = function(e) return(NA))
	Csb_dn_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% Csb_dn)]), error = function(e) return(NA))
	DKO_dn_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% DKO_dn)]), error = function(e) return(NA))

 	out <- list(Adh5_up = Adh5_up_geneLength, Csb_up = Csb_up_geneLength, DKO_up = DKO_up_geneLength, Adh5_dn = Adh5_dn_geneLength, Csb_dn = Csb_dn_geneLength, DKO_dn = DKO_dn_geneLength)
	return(out)}, mc.cores = 24)

# collapsing deg list
library(reshape2)

geneLengths <- pbmclapply(lengths, function(x){

	df_up <- tryCatch(melt(x[c(1,2,3)]), error = function(e) return(NA))
	if(length(df_up) > 1){
		colnames(df_up) <- c("log10geneLength", "Comparison")
		df_up$Comparison <- factor(df_up$Comparison, levels = c('Adh5_up', 'Csb_up', 'DKO_up'))
	}
	
	df_down <- tryCatch(melt(x[c(4,5,6)]), error = function(e) return(NA))
	if(length(df_down) > 1){
		colnames(df_down) <- c("log10geneLength", "Comparison")
		df_down$Comparison <- factor(df_down$Comparison, levels = c('Adh5_dn', 'Csb_dn', 'DKO_dn'))
	}	

	out <- list(up = df_up, down = df_down)
	return(out)
}, mc.cores = 24)

names(geneLengths) <- gsub('.txt', '', files)

kelvinny::dirCreate('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/figures/histograms')
library(ggplot2)
library(patchwork)

setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/figures/histograms')

files2 <- gsub('.txt', '', unlist(files))

outPaths <- paste0(files2, '_DEG_geneLength.png')

plotGeneLength <- function(i, g = geneLengths, o = outPaths, cols1 = c('#51a693', '#1b538a', '#ef3723'), cols2 = c('#51a693', '#1b538a', '#ef3723')){
	p1 <- ggplot(g[[i]]$up, aes(x = log10geneLength, color = Comparison)) + geom_density() + ggtitle("DEG Up") + scale_color_manual(values = cols1, drop = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(1, 7)
	p2 <- ggplot(g[[i]]$down, aes(x = log10geneLength, color = Comparison)) + geom_density() + ggtitle("DEG Down") + scale_color_manual(values = cols2, drop = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(1, 7)
	p3 <- p1 / p2
	ggsave(o[i], plot = p3, device = "png", h = 10, w = 10)
}

for (i in 1:length(outPaths)){ plotGeneLength(i)}

## repeat with PT subtybes
setwd("/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset_vsWT")
library(dplyr)
library(readr)
library(pbmcapply)
files <- as.list(list.files(pattern='.txt'))
degs <- lapply(files, read_tsv)

lengths <- pbmclapply(degs, function(x) {
	Adh5_up <- x %>% dplyr::filter(Adh5_pvals_adj < 0.05 & Adh5_logfoldchanges > 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	Csb_up <- x %>% dplyr::filter(Csb_pvals_adj < 0.05 & Csb_logfoldchanges > 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	DKO_up <- x %>% dplyr::filter(DKO_pvals_adj < 0.05 & DKO_logfoldchanges > 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	Adh5_dn <- x %>% dplyr::filter(Adh5_pvals_adj < 0.05 & Adh5_logfoldchanges < 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	Csb_dn <- x %>% dplyr::filter(Csb_pvals_adj < 0.05 & Csb_logfoldchanges < 0) %>% dplyr::select(X1) %>% unlist %>% as.character
	DKO_dn <- x %>% dplyr::filter(DKO_pvals_adj < 0.05 & DKO_logfoldchanges < 0) %>% dplyr::select(X1) %>% unlist %>% as.character

	Adh5_up_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% Adh5_up)]), error = function(e) return(NA))
	Csb_up_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% Csb_up)]), error = function(e) return(NA))
	DKO_up_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% DKO_up)]), error = function(e) return(NA))
	Adh5_dn_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% Adh5_dn)]), error = function(e) return(NA))
	Csb_dn_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% Csb_dn)]), error = function(e) return(NA))
	DKO_dn_geneLength <- tryCatch(log10(m$geneLength[which(m$external_gene_name %in% DKO_dn)]), error = function(e) return(NA))

 	out <- list(Adh5_up = Adh5_up_geneLength, Csb_up = Csb_up_geneLength, DKO_up = DKO_up_geneLength, Adh5_dn = Adh5_dn_geneLength, Csb_dn = Csb_dn_geneLength, DKO_dn = DKO_dn_geneLength)
	return(out)}, mc.cores = 24)

# collapsing deg list
library(reshape2)

geneLengths <- pbmclapply(lengths, function(x){

	df_up <- tryCatch(melt(x[c(1,2,3)]), error = function(e) return(NA))
	if(length(df_up) > 1){
		colnames(df_up) <- c("log10geneLength", "Comparison")
		df_up$Comparison <- factor(df_up$Comparison, levels = c('Adh5_up', 'Csb_up', 'DKO_up'))
	}
	
	df_down <- tryCatch(melt(x[c(4,5,6)]), error = function(e) return(NA))
	if(length(df_down) > 1){
		colnames(df_down) <- c("log10geneLength", "Comparison")
		df_down$Comparison <- factor(df_down$Comparison, levels = c('Adh5_dn', 'Csb_dn', 'DKO_dn'))
	}	

	out <- list(up = df_up, down = df_down)
	return(out)
}, mc.cores = 24)

files2 <- gsub('.txt', '', unlist(files))
setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/figures/histograms')

outPaths <- paste0(files2, '_DEG_geneLength.png')
for (i in 1:length(outPaths)){ plotGeneLength(i)}
```

### still testing
```R
setwd("/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset_vsWT")
library(dplyr)
library(readr)
library(pbmcapply)
files <- as.list(list.files(pattern='.txt'))
degs <- lapply(files, read_tsv)

lengths <- pbmclapply(degs, function(x) {
	y <- x %>% dplyr::select(X1) %>% unlist %>% as.character
	z <- log10(m$geneLength[match(y, m$external_gene_name)])
	x$geneLengths <- z
	return(x)}, mc.cores = 24)

lengths2 <- lapply(lengths, function(x) {
	z <- list()
	x <- x[order(x$geneLengths), ]
	x <- x[-grep('^Gm|^mt|[1-9]', x$X1), ]
	top10_shortest <- head(x, 20)
	top10_longest <- tail(x, 20)
	z[[1]] <- top10_shortest %>% dplyr::select(X1) %>% unlist %>% as.character
	z[[2]] <- top10_longest %>% dplyr::select(X1) %>% unlist %>% as.character
	return(z)
})


setwd("/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG")
library(dplyr)
library(readr)
library(pbmcapply)
files <- as.list(list.files(pattern='.txt'))
degs <- lapply(files, read_tsv)

lengths <- pbmclapply(degs, function(x) {
	y <- x %>% dplyr::select(X1) %>% unlist %>% as.character
	z <- log10(m$geneLength[match(y, m$external_gene_name)])
	x$geneLengths <- z
	return(x)}, mc.cores = 24)

lengths2 <- lapply(lengths, function(x) {
	z <- list()
	x <- x[order(x$geneLengths), ]
	top10_shortest <- head(x, 10)
	top10_longest <- tail(x, 10)
	z[[1]] <- top10_shortest %>% dplyr::select(X1) %>% unlist %>% as.character
	z[[2]] <- top10_longest %>% dplyr::select(X1) %>% unlist %>% as.character
	return(z)
})


```

### pathway analysis
```R
setwd("/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset")
library(dplyr)
library(readr)
library(pbmcapply)
degs <- read_tsv('PT-4-DKO_topmarkers.txt')

library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
m <- getBM(attributes=c("external_gene_name", "entrezgene_id"), mart=mart)

library(org.Mm.eg.db)
library(clusterProfiler)

topup <- degs %>% filter(pvals_adj < 0.01 & logfoldchanges > 1) %>% dplyr::select(X1) %>% unlist %>% as.character
topup <- m$entrezgene_id[m$external_gene_name %in% topup]
topup <- as.character(topup[!is.na(topup)])

universeID = degs$X1 %>% unlist %>% as.character
universeID <- m$entrezgene_id[m$external_gene_name %in% universeID]
universeID <- as.character(universeID[!is.na(universeID)])

res <- enrichGO(gene = topup, 
          OrgDb = org.Mm.eg.db,
          keyType = 'ENTREZID',
          ont = "BP",
          universe = universeID,
          pAdjustMethod = "BH",
          pvalueCutoff = .05,
          qvalueCutoff = .05,
          readable = TRUE)
res <- simplify(res)
write_tsv(as.data.frame(res), 'PT-4-DKO_GOBP.txt')

res <- enrichKEGG(gene = topup, 
          organism = "mmu",
          keyType = 'ncbi-geneid',
          universe = universeID,
          pAdjustMethod = "BH")

geneID_splitted <- strsplit(res@result$geneID, "/")
idx <- lapply(geneID_splitted, match, m$entrezgene_id)
return_name <- function(x , y = m, n = 1) y[x, n]
genename_matched_and_splitted <- lapply(idx, return_name, y = m, n = 1)
genename_matched <- lapply(genename_matched_and_splitted, paste, collapse="/")
res@result$geneID <- unlist(genename_matched)

write_tsv(as.data.frame(res), 'PT-4-DKO_KEGG.txt')
```

## get the genes used from KEGG analysis to do score gene sets
```R
setwd("/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset")
library(KEGGREST)
library(pbmcapply) 
queryList <- list(c('mmu05169', 'mmu04612', 'mmu04940', 'mmu05140', 'mmu05330', 'mmu05320', 'mmu05332', 'mmu05150', 'mmu05145', 'mmu05416'),
 c('mmu05321', 'mmu04658', 'mmu05310', 'mmu05168', 'mmu05322', 'mmu05152', 'mmu04659', 'mmu05164', 'mmu04145', 'mmu05167'),
 c('mmu05166', 'mmu04672', 'mmu04640', 'mmu03320', 'mmu04668', 'mmu04610', 'mmu05323', 'mmu04514', 'mmu04620', 'mmu00830'),
 c('mmu05170', 'mmu04218', 'mmu04115', 'mmu05133', 'mmu05160', 'mmu00071', 'mmu05163', 'mmu04622', 'mmu04621', 'mmu00980'),
 c('mmu04380'))

query <- lapply(queryList, keggGet)

queryGenes <- lapply(query, function(x){
	y <- lapply(x, function(z){
		z <- z[['GENE']]
		z <- z[c(FALSE, TRUE)]
		z <- strsplit(z, ';')
		z <- do.call(rbind, z)
		z <- unlist(z[,1])
		z <- z[!grepl(' ', z)]	
		return(z)
	})
	return(y)
})

queryGenes <- lapply(rapply(queryGenes, enquote, how="unlist"), eval)
names(queryGenes) <- c('Epstein-Barr virus infection', 'Antigen processing and presentation', 'Type I diabetes mellitus', 'Leishmaniasis', 'Allograft rejection', 'Autoimmune thyroid disease', 'Graft-versus-host disease', 'Staphylococcus aureus infection', 'Toxoplasmosis', 'Viral myocarditis', 'Inflammatory bowel disease (IBD)', 'Th1 and Th2 cell differentiation', 'Asthma', 'Herpes simplex virus 1 infection', 'Systemic lupus erythematosus', 'Tuberculosis', 'Th17 cell differentiation', 'Influenza A', 'Phagosome', 'Kaposi sarcoma-associated herpesvirus infection', 'Human T-cell leukemia virus 1 infection', 'Intestinal immune network for IgA production', 'Hematopoietic cell lineage', 'PPAR signaling pathway', 'TNF signaling pathway', 'Complement and coagulation cascades', 'Rheumatoid arthritis', 'Cell adhesion molecules (CAMs)', 'Toll-like receptor signaling pathway', 'Retinol metabolism', 'Human immunodeficiency virus 1 infection', 'Cellular senescence', 'p53 signaling pathway', 'Pertussis', 'Hepatitis C', 'Fatty acid degradation', 'Human cytomegalovirus infection', 'RIG-I-like receptor signaling pathway', 'NOD-like receptor signaling pathway', 'Metabolism of xenobiotics by cytochrome P450', 'Osteoclast differentiation')

saveRDS(queryGenes, 'significant_KEGGgenes_PT4.rds')
```

```R
library(Seurat) # version 3.1
library(dplyr)
library(ggplot2)
setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/')
pt <- readRDS(file = "scanpy/out/Lee_Sept_new_pt.rds")
queryGenes <- readRDS('scanpy/out/DEG/PT_subset/significant_KEGGgenes_PT4.rds')


pt2 <- AddModuleScore(pt, queryGenes, name = names(queryGenes))
pt4 <- pt2@meta.data %>% filter(pt_celltype == 'PT-4-DKO')
pt4 <- split(pt4, pt4$group)
pt4 <- lapply(pt4, function(x){
	x <- x[, -c(1:ncol(pt@meta.data))]
	return(x)
})
pt4 <- lapply(pt4, function(x){
	x <- colMeans(x)
	return(x)
})
pt4 <- do.call(rbind, pt4)
pt4 <- pt4[c(4,1,2,3),]
colnames(pt4) <- names(queryGenes)

library(kelvinny)

annotation_col = data.frame(row.names = c('WT', 'Adh5', 'Csb', 'DKO'), group = c('WT', 'Adh5', 'Csb', 'DKO'))
g.palette <- c('#3a4854', '#51a693', '#1b538a', '#ef3723')
names(g.palette) <- c('WT', 'Adh5', 'Csb', 'DKO')
anno_color <- list(group = g.palette)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
plotmat <- apply(pt4,2,range01)

pdf('scanpy/out/DEG/PT_subset/significant_KEGGgenes_PT4.pdf')
plotHeat(t(plotmat),
			cellheight = 10,
			cellwidth = 10,
			treeheight_col = 10,
			treeheight_row = 10,
			scale = 'none',
			color = viridis::viridis(50),
			annotation_col = annotation_col,
			annotation_color = anno_color,
			cluster_cols = FALSE)
dev.off()

```

### DEG dot plot for PT
```R
# and now to tabulate in R
library(dplyr)
library(readr)
library(tibble)
library(pbmcapply)

setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset_vsWT')
files <- as.list(list.files(pattern = '.txt'))
names(files) <- gsub('.txt', '', list.files(pattern = '.txt'))
degs <- lapply(files, read_tsv)

# lfc threshold of 1
res <- pbmclapply(degs, function(x){
	Adh5_up <- x %>% filter(Adh5_pvals_adj < 0.05 & Adh5_logfoldchanges >= 1) %>% select(X1) %>% unlist %>% length
	Csb_up <- x %>% filter(Csb_pvals_adj < 0.05 & Csb_logfoldchanges >= 1) %>% select(X1) %>% unlist %>% length
	DKO_up <- x %>% filter(DKO_pvals_adj < 0.05 & DKO_logfoldchanges >= 1) %>% select(X1) %>% unlist %>% length
	Adh5_dn <- x %>% filter(Adh5_pvals_adj < 0.05 & Adh5_logfoldchanges <= -1) %>% select(X1) %>% unlist %>% length
	Csb_dn <- x %>% filter(Csb_pvals_adj < 0.05 & Csb_logfoldchanges <= -1) %>% select(X1) %>% unlist %>% length
	DKO_dn <- x %>% filter(DKO_pvals_adj < 0.05 & DKO_logfoldchanges <= -1) %>% select(X1) %>% unlist %>% length
	y <- data.frame(Adh5_up = Adh5_up, Csb_up = Csb_up, DKO_up = DKO_up, Adh5_dn = Adh5_dn, Csb_dn = Csb_dn, DKO_dn = DKO_dn)
	return(y)
}, mc.cores = 24)

results <- do.call(rbind, res)
results

# Lee wants these as a dot plot
results_up <- results[,1:3]
results_down <- results[,4:6]

results_up$id = rownames(results)
results_down$id = rownames(results)

library(reshape2)
results_up = melt(results_up)
results_down = melt(results_down)

colnames(results_up) <- c('Cell type', 'Group', 'No. of DEGs')
colnames(results_down) <- c('Cell type', 'Group', 'No. of DEGs')

results <- rbind(results_up, results_down)

g <- ggplot(results, aes(x=Group, y=`Cell type`, size=`No. of DEGs`, colour=`No. of DEGs`)) + 
	geom_point(shape = 16) + 
	scale_colour_gradientn(colors=c('#fff5f0', '#fee0d2', '#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d')) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), legend.key=element_blank()) + 
	scale_size(range=c(0,4))
ggsave("../../../figures/dotplot/DEGs_PT1.pdf", g, h = 3, w = 4, useDingbats = FALSE)

g <- ggplot(results, aes(x=Group, y=`Cell type`, size=`No. of DEGs`, colour=`No. of DEGs`)) + 
	geom_point(shape = 16) + 
	scale_colour_gradientn(colors=c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), legend.key=element_blank()) + 
	scale_size(range=c(0,4))
ggsave("../../../figures/dotplot/DEGs_PT2.pdf", g, h = 3, w = 4, useDingbats = FALSE)


```


```R
library(Seurat) # version 3.1
library(dplyr)
library(ggplot2)
setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/')
pt <- readRDS(file = "scanpy/out/Lee_Sept_new_pt.rds")

annotation_col = data.frame(row.names = c('PT-0', 'PT-1', 'PT-2', 'PT-3', 'PT-4-DKO', 'PT-5', 'PT-6', 'PT-7-Novel'), group = c('PT-0', 'PT-1', 'PT-2', 'PT-3', 'PT-4-DKO', 'PT-5', 'PT-6', 'PT-7-Novel'))
g.palette <- c('#0570b0', '#3690c0', '#74a9cf', '#7ADFC4', '#F28E2b', '#53A592', '#02818a', '#034e7b')
names(g.palette) <- c('PT-0', 'PT-1', 'PT-2', 'PT-3', 'PT-4-DKO', 'PT-5', 'PT-6', 'PT-7-Novel')
anno_color <- list(group = g.palette)

library(kelvinny)
library(ktplots) 
# alcohol metabolic process
gs1=c('Aldh1a7', 'Aldh1a1', 'Hmgcs2')
# regulation of apoptotic process
gs2=c('Gdf15', 'Hspb1', 'Btg2', 'Cdkn1a', 'Plac8', 'Cryab', 'Phlda3', 'Sfn', 'Lgals3')
# response to interferon-beta
gs3=c('Gbp5', 'Iigp1', 'Irgm1', 'Bst2', 'Stat1', 'Tgtp1', 'Ifit1', 'Ifit3', 'Gbp3', 'Irf1', 'Xaf1', 'Gbp2')
# antigen processing and presentation of peptide antigen
gs4=c('H2-K1', 'Tap2', 'H2-DMa', 'H2-Ab1', 'H2-Aa', 'H2-Q7', 'H2-Eb1', 'Cd74', 'H2-DMb1', 'H2-D1')
# protein activation cascade
gs5=c('C3', 'Fgg', 'Fgb', 'C4b', 'Cfi', 'Fga')

# create a vector to store the cell cluster
clusters <- pt$pt_celltype

plotmat <- GetAssayData(pt)

plotmat1 <- plotmat[gs1,]
plotmat2 <- plotmat[gs2,]
plotmat3 <- plotmat[gs3,]
plotmat4 <- plotmat[gs4,]
plotmat5 <- plotmat[gs5,]

plotmat1 <- data.frame(cluster = clusters, Matrix::t(plotmat1), check.names = FALSE)
plotmat2 <- data.frame(cluster = clusters, Matrix::t(plotmat2), check.names = FALSE)
plotmat3 <- data.frame(cluster = clusters, Matrix::t(plotmat3), check.names = FALSE)
plotmat4 <- data.frame(cluster = clusters, Matrix::t(plotmat4), check.names = FALSE)
plotmat5 <- data.frame(cluster = clusters, Matrix::t(plotmat5), check.names = FALSE)

plotmats <- list(plotmat1, plotmat2, plotmat3, plotmat4, plotmat5)

plotmats <- lapply(plotmats, function(x){
	y <- split(x, x$cluster)
	y <- lapply(y, function(z){
		z <- z[,-1]
		z <- colMeans(z)
		return(z)
	})
	y <- do.call(rbind, y)
	return(y)
})

plotmats <- lapply(plotmats, function(x) apply(x,2,range01))
pathwaynames <- c('alcohol metabolic process', 'regulation of apoptotic process', 'response to interferon-beta', 'antigen processing and presentation of peptide antigen', 'protein activation cascade')

pdf('scanpy/out/DEG/PT_subset/PT_function_genes.pdf')
for (i in 1:5){
	plotHeat(t(plotmats[[i]]),
			cellheight = 10,
			cellwidth = 10,
			treeheight_col = 10,
			treeheight_row = 10,
			scale = 'none',
			color = viridis::viridis(50),
			annotation_col = annotation_col,
			annotation_color = anno_color,
			cluster_cols = FALSE,
			main = pathwaynames[i])
}
dev.off()

```
### Lee wanted down-regulated bubblecharts
### bubblechart for pt
```R
library(dplyr)
library(readr)
# setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset')
setwd('~/Documents/Clatworthy_scRNAseq/Lee_analysis/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset')
files <- as.list(paste0(c(0:7), '.txt'))

degs <- lapply(files, read_tsv)
names(degs) <- gsub(".txt", "", files)

## lets start with DKO?
# extract the significant genes
gene_df <- lapply(degs, function(x){
	df <- x %>% dplyr::filter(DKO_pvals_adj < 0.01) %>% dplyr::select(X1, Adh5_logfoldchanges, Csb_logfoldchanges, DKO_logfoldchanges, Adh5_pvals_adj, Csb_pvals_adj, DKO_pvals_adj)
	colnames(df) <- c("genes", 'Adh5_logfoldchanges', 'Csb_logfoldchanges', 'DKO_logfoldchanges', 'Adh5_pvals_adj', 'Csb_pvals_adj', 'DKO_pvals_adj')
	df$ranking <- df$DKO_logfoldchanges
	df <- df[order(df$ranking), ]
	df$order <- order(df$ranking)
	top <- df %>% top_n(25) %>% dplyr::select(genes)
	bottom <- df %>% top_n(-25) %>% dplyr::select(genes)
	keep <- unlist(c(top, bottom))
	df$top_genes <- df$genes
	df$top_genes[!df$top_genes %in% keep] <- NA
	return(df)
})

gene_df <- lapply(gene_df, function(x){
	x <- x %>% head(25)
	return(x)
})

bubblefall <- function(df){
	require(ggplot2)
	require(ggrepel)
	p <- ggplot() + geom_hline(aes(yintercept=0), colour = 'grey', linetype = 'dashed', size = 0.5) + labs(x = NULL) +
	labs(y = expression('Log'[2]*'Fold-change')) +
	theme_bw() +
	guides(size=FALSE) +
	theme(axis.line = element_line(colour = 'black'),
     		legend.position = 'bottom', 
        	panel.grid.major = element_blank(),
        	panel.grid.minor = element_blank(),
        	panel.border = element_blank(),
        	panel.background = element_blank(),
        	axis.line.x=element_blank(),
        	axis.text.x=element_blank(),
        	axis.ticks.x=element_blank()) +
     # ggtitle('Core Enrichment Genes\n\nCaproni et al. and Mosca et al. Chemical Adjuvantation gene sets\ni.m. Alum/CpG/MF59 vs. PBS') + 
     theme(plot.title = element_text(lineheight = 0.8, hjust = 0.5))

	p <- p + geom_point(data = df, aes(x = order, y = Adh5_logfoldchanges, size = -log10(Adh5_pvals_adj), fill = 'Adh5', colour = 'Adh5'), shape = 21, alpha = 0.7, stroke = 0.1) + 
		 geom_point(data = df, aes(x = order, y = Csb_logfoldchanges, size = -log10(Csb_pvals_adj), fill = 'Csb', colour = 'Csb'), shape = 21, alpha = 0.7, stroke = 0.1) +
		 geom_point(data = df, aes(x = order, y = DKO_logfoldchanges, size = -log10(DKO_pvals_adj), fill = 'DKO', colour = 'DKO'), shape = 21, alpha = 0.7, stroke = 0.1)
	p <- p + scale_size(range = c(1,4))

    p <- p + geom_text_repel(aes(x = df$order, y = df$DKO_logfoldchanges, label = df$top_genes, size=3), box.padding = unit(0.3, 'lines'), segment.size = 0.2, force=0.5)

	p <- p + scale_fill_manual(name = '', values = c('Adh5' = '#51a693', 'Csb' = '#1b538a', 'DKO' = '#ef3723'), breaks=c('Adh5', 'Csb', 'DKO'), guide = 'legend') +
         scale_colour_manual(name = '', values = c('Adh5' = '#FFFFFF', 'Csb' = '#FFFFFF', 'DKO' = '#FFFFFF'), breaks=c('Adh5', 'Csb', 'DKO'), guide = 'legend')

	# p <- p + theme(legend.title = element_blank())
	return(p)}

small_legend <- function(...){
	small_legend_theme <- 
	theme(legend.title = element_text(size = 4), 
		legend.text  = element_text(size = 4),
		legend.key.size = unit(0.01, "lines"), ...)
	return(small_legend_theme)
}

small_guide <- function(...){
	small_guide <- 
	guides(shape = guide_legend(override.aes = list(size = 0.1)), 
		color = guide_legend(override.aes = list(size = 0.1)), ...)
	return(small_guide)
}


small_axis <- function(...){
	axis <- theme(text = element_text(size=8), axis.text = element_text(size=8), axis.line = element_line(size = 0.1), axis.ticks = element_line(size = 0.1), ...)
	return(axis)
}

topleft_legend <- function(...){
	legend <- theme(legend.position = c(.01, .99), legend.justification = c('left', 'top'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

bottomleft_legend <- function(...){
	legend <- theme(legend.position = c(.01, .01), legend.justification = c('left', 'bottom'), legend.box.just = "left", legend.margin = margin(5, 5, 5, 5), ...)
	return(legend)
}

pdf("../../../figures/dotplot/PT-0_bottomgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['0']]) +
	ggtitle('PT-0') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-1_bottomgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['1']]) +
	ggtitle('PT-1') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-2_bottomgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['2']]) +
	ggtitle('PT-2') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-3_bottomgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['3']]) +
	ggtitle('PT-3') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-4-DKO_bottomgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['4']]) + 
	ggtitle('PT-4-DKO') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-5_bottomgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['5']]) +
	ggtitle('PT-5') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-6_bottomgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['6']]) +
	ggtitle('PT-6') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	bottomleft_legend()
dev.off()

pdf("../../../figures/dotplot/PT-7-Novel1_bottomgenes.pdf", h = 3, w = 4)
bubblefall(gene_df[['7']]) +
	ggtitle('PT-7-Novel1') + 
	small_legend() + theme(legend.text  = element_text(size = 7)) +
	small_guide(fill = guide_legend(override.aes = list(size = 3))) +
	small_axis() + 
	topleft_legend()
dev.off()

```

### Lee wanted a volcano plot for T/NK and MNPs
```R
library(EnhancedVolcano)
library(dplyr)
library(readr)
# setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG/PT_subset')
setwd('~/Documents/Clatworthy_scRNAseq/Lee_analysis/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG')
files <- as.list(c('MNP.txt', 'T_NK.txt'))
degs <- lapply(files, read_tsv)
names(degs) <- gsub(".txt", "", files)


EnhancedVolcano(degs[[1]],
    lab = degs[[1]]$X1,
    x = 'DKO_logfoldchanges',
    y = 'DKO_pvals_adj',
    xlim = c(-5, 5),
	title = 'DKO MNP vs WT MNP')


EnhancedVolcano(degs[[1]],
    lab = degs[[1]]$X1,
    x = 'DKO_logfoldchanges',
    y = 'DKO_pvals_adj',
    FCcutoff = .5,
    pCutoff = 10e-5,
    ylim = c(2, 70),
    xlim = c(-3, 5),
	title = 'DKO MNP vs WT MNP')

EnhancedVolcano(degs[[2]],
    lab = degs[[2]]$X1,
    x = 'DKO_logfoldchanges',
    y = 'DKO_pvals_adj',
    FCcutoff = .5,
    pCutoff = 10e-5,
    ylim = c(2, 70),
    xlim = c(-5, 5),
	title = 'DKO T_NK vs WT T_NK')

EnhancedVolcano(degs[[1]],
    lab = degs[[1]]$X1,
    x = 'DKO_logfoldchanges',
    y = 'DKO_pvals_adj',
    # selectLab = 'Ccl2',
    labSize = 2,
    FCcutoff = .5,
    pCutoff = 5e-2,
    ylim = c(0, 70),
    xlim = c(-6, 6),
	title = 'DKO MNP vs WT MNP')

EnhancedVolcano(degs[[1]],
    lab = degs[[1]]$X1,
    x = 'DKO_logfoldchanges',
    y = 'DKO_pvals_adj',
    selectLab = c('Ccl2', 'Ccl7', 'Ccl8'),
    FCcutoff = .5,
    pCutoff = 5e-2,
    ylim = c(0, 70),
    xlim = c(-6, 6),
	title = 'DKO MNP vs WT MNP')


EnhancedVolcano(degs[[1]],
    lab = degs[[1]]$X1,
    x = 'DKO_logfoldchanges',
    y = 'DKO_pvals_adj',
    selectLab = 'Ccl2',
    FCcutoff = .5,
    pCutoff = 5e-2,
    ylim = c(0, 70),
    xlim = c(-5, 5),
	title = 'DKO T_NK vs WT T_NK')


```

```R
library(Seurat)
library(reticulate)
setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/')
sc = import('scanpy')
tnk <- ReadH5AD('out/Lee_T_NK.h5ad')
adata = sc$read_h5ad('out/Lee_T_NK.h5ad')
tnk@meta.data = adata$obs
Idents(tnk) <- "celltype_T_NK"
T <- subset(tnk, idents = c("NK cell", "NK cell - cytotoxic"), invert = TRUE)


library(dplyr)
library(readr)
library(tibble)
library(pbmcapply)

setwd('/home/jovyan/Lee_KJ_Kidney_10x_sep19/scanpy/out/DEG')
files <- as.list('T.txt')
names(files) <- 'T'
degs <- lapply(files, read_tsv)

makeGeneList <- function(gl){
	gl1 <- gl %>% dplyr::select(X1, Adh5_logfoldchanges, Adh5_pvals)
	gl1$neglog10pval <- -log10(gl1$Adh5_pvals)
	rank <- unlist(gl1$neglog10pval*sign(gl1$Adh5_logfoldchanges))
	names(rank) <- gl1$X1
	rank <- rev(sort(rank))
	gl1 <- rank

	gl2 <- gl %>% dplyr::select(X1, Csb_logfoldchanges, Csb_pvals)
	gl2$neglog10pval <- -log10(gl2$Csb_pvals)
	rank <- unlist(gl2$neglog10pval*sign(gl2$Csb_logfoldchanges))
	names(rank) <- gl2$X1
	rank <- rev(sort(rank))
	gl2 <- rank

	gl3 <- gl %>% dplyr::select(X1, DKO_logfoldchanges, DKO_pvals)
	gl3$neglog10pval <- -log10(gl3$DKO_pvals)
	rank <- unlist(gl3$neglog10pval*sign(gl3$DKO_logfoldchanges))
	names(rank) <- gl3$X1
	rank <- rev(sort(rank))
	gl3 <- rank
	return(list(gl1, gl2, gl3))
}
test <- makeGeneList(degs[[1]])

library(fgsea)
h <- as.list(kelvinny::parse_gmt("/home/jovyan/Prostate_analysis/scanpy/dataset/h.all.v7.0.symbols.gmt"))
h <- lapply(h, function(x) {x <- x[-1]; x <- x[!is.na(x)]; return(x)})

# do the symbol conversion
library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
m <- getBM(attributes=c("external_gene_name", "mmusculus_homolog_associated_gene_name"), mart=mart)

h <- lapply(h , function(x){
	y <- m$mmusculus_homolog_associated_gene_name[m$external_gene_name %in% x]
	y <- y[-which(y == "")]	
	return(y)
}) 

# res <- lapply(test, function(x) fgsea(pathways = h, stats = x, nperm=10000, minSize = 0, maxSize = 1000))

# filter the genes for the gsea?	
up_cutOff = 1.5	
down_cutOff = -1.5	
geneList <- pbmclapply(degs, function(x){	
	Adh5 <- x %>% dplyr::filter(!between(Adh5_logfoldchanges, down_cutOff, up_cutOff)) %>% dplyr::select(X1, Adh5_logfoldchanges, Adh5_pvals)	
	Csb <- x %>% dplyr::filter(!between(Csb_logfoldchanges, down_cutOff, up_cutOff)) %>% dplyr::select(X1, Csb_logfoldchanges, Csb_pvals)	
	DKO <- x %>% dplyr::filter(!between(DKO_logfoldchanges, down_cutOff, up_cutOff)) %>% dplyr::select(X1, DKO_logfoldchanges, DKO_pvals)	
		
	geneList_l <- list(Adh5 = Adh5, Csb = Csb, DKO = DKO)	
	geneList_l <- pbmclapply(geneList_l, function(y) {	
		y$neglog10pval <- -log10(y[,3, drop = TRUE])	
		rank <- unlist(y$neglog10pval*sign(y[,2, drop = TRUE]))	
		rank[which(rank == Inf)] <- -log10(10^-308)	
		rank[which(rank == -Inf)] <- log10(10^-308)	
		names(rank) <- y[,1, drop = TRUE]	
		rank <- rev(sort(rank))	
		# if there's are Inf values, just need to change those to -log10(10^-308) or inverse of that	
		return(rank)	
	}, mc.cores = 4)	
	return(geneList_l)	
}, mc.cores = 4)	
library(fgsea)	
result <- list()	
for(i in 1:length(geneList)){	
	result[[i]] <- lapply(geneList[[i]], function(x) fgsea(pathways = h, stats = x, nperm=10000, minSize = 0, maxSize = 1000))		
}	

for(i in 1:length(geneList)){	
	result[[i]] <- lapply(result[[i]], function(x){	
		x <- as.data.frame(x)
		x$ranking <- -log10(x$pval)*sign(x$NES)		
		x <- x[order(x$ranking), ]
		return(x)})}

result[[1]][[1]]$group <- 'Adh5'
result[[1]][[2]]$group <- 'Csb'
result[[1]][[3]]$group <- 'DKO'

result2 <- lapply(result, function(x) {	
	y <- do.call(rbind, x)	
	y$group <- factor(y$group, levels = c('Adh5', 'Csb', 'DKO')) 	
	return(y)	
})	

plotGSEA_Hallmark <- function(gsea, group_ref = NULL, cols = NULL, newlabels = NULL) {
	require(ggplot2)
	# gsea <- result2[[1]]
	# group_ref = 'DKO'
	gsea$NES[which(is.na(gsea$NES))] <- 0
	gsea$ranking[which(is.na(gsea$ranking))] <- 0
	gsea <- gsea[order(gsea$ranking),]		
	gsea_spl <- split(gsea, gsea$group)
	if(!is.null(group_ref)){
		gsea_spl[[group_ref]] <- gsea_spl[[group_ref]][order(gsea_spl[[group_ref]]$ranking),]
		gsea_spl[[group_ref]]$ranking <- gsea_spl[[group_ref]]$ranking*999
	} else {
		gsea_spl[[2]] <- gsea_spl[[2]]$ranking*999
	}
	names(gsea_spl) <- NULL

	gsea <- do.call(rbind, gsea_spl)
	gsea <- gsea[order(gsea$ranking), ]
	gsea$pathway <- gsub("HALLMARK_|", "", gsea$pathway)
	gsea$group[which(gsea$pval >= 0.05 & gsea$padj >= 0.25)] <- "NotSig"
	
	x_lim_min <- abs(ceiling(min(-log10(gsea$pval))))
	x_lim_max <- abs(ceiling(max(-log10(gsea$pval))))

	if(x_lim_min > x_lim_max){
		xval1 <- x_lim_min * -1
		xval2 <- x_lim_min
	} else {
		xval1 <- x_lim_max * -1
		xval2 <- x_lim_max
	}

	if(!is.null(cols)){
		gg_color_hue <- function(n) {
			hues = seq(15, 375, length = n + 1)
			hcl(h = hues, l = 65, c = 100)[1:n]
		}
		cols. = gg_color_hue(dplyr::n_distinct(gsea$group, na.rm = TRUE))
	} else {
		cols. = cols
	}

	g <- ggplot(gsea, aes(x = -log10(pval)*sign(NES), y = reorder(pathway, ranking), col = group, size = NES)) + 
		geom_point() + 
		labs(x = expression(paste("Signed", " -log" ["10"], "P-value")), y = "Hallmarks") +
		theme_bw() +
		geom_vline(xintercept = 0) +
		geom_vline(xintercept = -log10(0.25)) +
		geom_vline(xintercept = -log10(0.25)*-1) +
		xlim(xval1, xval2) +
		scale_radius(range = c(1,4)) +
		theme(panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(), 
			panel.background = element_blank(), 
			axis.line = element_blank(), 
			axis.ticks = element_blank())
	if(!is.null(newlabels))
		g <- g + scale_color_manual(values = cols., na.value = 'grey90', drop = FALSE, labels = newlabels)
	else {
		g <- g + scale_color_manual(values = cols., na.value = 'grey90', drop = FALSE)
	}
	g$data <- g$data[order(g$data$group, na.last = FALSE), ]
	return(g)
}

for(i in 1:length(geneList)){	
p <- plotGSEA_Hallmark(result2[[i]], group_ref = "DKO", cols = c("#1f77b4", "#ff7f0e", "#279e68"), newlabels = c("Adh5", "Csb", "DKO", "NotSig")) + ggtitle('T')
ggsave(paste0("T_gsea_hallmarks.pdf"), plot = p, w = 8.5)
}
result3 <- result2[[1]] %>% dplyr::filter(pathway %in% c('HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY', 'HALLMARK_IL2_STAT5_SIGNALING', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_TGF_BETA_SIGNALING'))
p <- plotGSEA_Hallmark(result3, group_ref = "DKO", cols = c("#1f77b4", "#ff7f0e", "#279e68", 'grey90'), newlabels = c("Adh5", "Csb", "DKO", "NotSig")) + ggtitle('T') + scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#279e68", 'grey90', na.value = 'grey90', drop = FALSE))
ggsave(paste0("T_gsea_hallmarks.pdf"), plot = p, w = 8, h = 3)
```