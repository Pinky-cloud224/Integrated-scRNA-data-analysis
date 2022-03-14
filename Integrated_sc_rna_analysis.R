library(Seurat)
library(seuratData)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(cowplot)
library(scran)

sc_rna_data.data <- Read10X(data.dir = "/home/pinky/sc_rna_seq/filtered_gene_bc_matrices/hg19")
sc_rna_data <- CreateSeuratObject(counts=sc_rna_data.data,project="singlecellrna",min.cells=2, min.features=150)
sc_rna_data

data_1 <- "/home/pinky/sc_rna_seq/filtered_gene_bc_matrices/hg19"

data_2 <- "/home/pinky/sc_rna_seq/filtered_gene_bc_matrices/hg19"


sc_rna_data.data_1 <- Read10X(data.dir= data_1)

sample_1_control <- CreateSeuratObject(counts=sc_rna_data.data_1, project="SingleCellRNAData_1_control", min.cells=350, min.features=250)


sc_rna_data.data_2 <- Read10X(data.dir= data_2)

sample_2_test <- CreateSeuratObject(counts=sc_rna_data.data_2, project="SingleCellRNAData_2_test", min.cells=650, min.features=250)



merged.sc_rna_data <-merge (x=sample_1_control, y=sample_2_test, add.cell.ids=c("sample_1_control","sample_2_test"), project="SCRNA_Data")


# Adding metadata
library(stringr)
sample <- names(merged.sc_rna_data@active.ident)
sample_detect <- ifelse(str_detect(sample,"sample_1_control"),"sample_1_control","sample_2_test")

merged.sc_rna_data@meta.data$sample <-sample_detect
Idents(object=merged.sc_rna_data) <- "sample"


# Quality control


merged.sc_rna_data$log10GenesPerUMI <- log10(merged.sc_rna_data$nFeature_RNA)/log10(merged.sc_rna_data$nCount_RNA)

## Compute percent mito ratio
merged.sc_rna_data$mitoRatio <- PercentageFeatureSet(object=merged.sc_rna_data, pattern="^MT-")
merged.sc_rna_data$mitoRatio <- merged.sc_rna_data@meta.data$mitoRatio / 100

#mitochondrial percentage
merged.sc_rna_data[["percent.mt"]] <- PercentageFeatureSet(merged.sc_rna_data, pattern = "^MT-")
VlnPlot(merged.sc_rna_data, features = "percent.mt",split.plot=TRUE, split.by="sample")
metadata <- merged.sc_rna_data@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "sample_1_control"
metadata$sample[which(str_detect(metadata$cells, "^test_"))] <- "sample_2_test"

merged.sc_rna_data@meta.data <- metadata
save(merged.sc_rna_data, file="/home/pinky/sc_rna_seq/merged_filtered_sc_rna_data.RData")

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(merged.sc_rna_data, 
                     vars = c("ident", "sample")) %>%
  group_by(sample) %>%
  dplyr::count(ident) %>% 
  spread(ident, n) 

# View table
View(n_cells)

# Determine metrics to plot present in merged.sc_rna_data
metrics <-  c("nUMI", "nGene", "mitoRatio")

FeaturePlot(merged.sc_rna_data, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
## q10 translates to the 10% of cells

########################################################################################################

### Normalizing data
sc_rna_data.list <-SplitObject(merged.sc_rna_data, split.by="sample")

## performing standard prepocessing on each object

for (i in 1:length(sc_rna_data.list)) {
  sc_rna_data.list[[i]] <- NormalizeData(sc_rna_data.list[[i]], verbose=FALSE)
  sc_rna_data.list[[i]] <- subset(sc_rna_data.list[[i]], downsample=3000)
  sc_rna_data.list[[i]] <- FindVariableFeatures(
    sc_rna_data[[i]], selection.method="vst",
    nfeatures=2500, verbose=FALSE
  )
}

features <- SelectIntegrationFeatures(object.list=sc_rna_data.list)
sc_rna_data.list <- lapply (X = sc_rna_data.list, FUN=function(x) {
  x <- ScaleData(x, features=features, verbose=FALSE)
  x <- RunPCA(x, features= features, verbose=FALSE)
  
})


# Integration of Data
anchors <- FindIntegrationAnchors(object.list=sc_rna_data.list)


#### Data cleaning and dimentionality reduction for visualization
merged.sc_rna_data <- ScaleData(merged.sc_rna_data, verbose=FALSE)
merged.sc_rna_data <- FindVariableFeatures(merged.sc_rna_data, selection.method="vst",
                                           nfeatures=2500,
                                           verbose=FALSE)

merged.sc_rna_data <- RunPCA(merged.sc_rna_data, npcs=50, verbose=FALSE)
merged.sc_rna_data <- RunUMAP(merged.sc_rna_data, reduction="pca",dims=1:50)

columns <- c(paste0("PC_", 1:10),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(merged.sc_rna_data, 
                     vars = columns)

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:10), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide ="none", 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y))+
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)


#### Computing the nearest neighbor 

merged.sc_rna_data <- FindNeighbors(merged.sc_rna_data, reduction="pca",dims=1:50)

##### Clustering the cells

merged.sc_rna_data <- FindClusters(merged.sc_rna_data, resolution=0.5)


merged.sc_rna_data = RunTSNE(merged.sc_rna_data)
DimPlot (merged.sc_rna_data, reduction = "tsne", label = TRUE)

##################### Finding markers for the cluster
# for cluster 1
cluster1.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 0, min.pct = 0.25)
head(cluster1.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster1.markers)[1], row.names(cluster1.markers)[2]))

# for cluster 2
cluster2.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 1, min.pct = 0.25)
head(cluster2.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster2.markers)[1], row.names(cluster2.markers)[2]))

# for cluster 3
cluster3.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 2, min.pct = 0.25)
head(cluster3.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster3.markers)[1], row.names(cluster3.markers)[2]))

# for cluster 4
cluster4.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 3, min.pct = 0.25)
head(cluster4.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster4.markers)[1], row.names(cluster4.markers)[2]))


# for cluster 5
cluster5.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 4, min.pct = 0.25)
head(cluster5.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster5.markers)[1], row.names(cluster5.markers)[2]))

# for cluster 6
cluster6.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 5, min.pct = 0.25)
head(cluster6.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster6.markers)[1], row.names(cluster6.markers)[2]))

# for cluster 7
cluster7.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 6, min.pct = 0.25)
head(cluster7.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster7.markers)[1], row.names(cluster7.markers)[2]))


# for cluster 8
cluster8.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 7, min.pct = 0.25)
head(cluster8.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster8.markers)[1], row.names(cluster8.markers)[2]))


# for cluster 9
cluster9.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 8, min.pct = 0.25)
head(cluster9.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster9.markers)[1], row.names(cluster9.markers)[2]))

# for cluster 10
cluster10.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 9, min.pct = 0.25)
head(cluster10.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster10.markers)[1], row.names(cluster10.markers)[2]))

# for cluster 11
cluster11.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 10, min.pct = 0.25)
head(cluster11.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster11.markers)[1], row.names(cluster11.markers)[2]))


# for cluster 12
cluster12.markers <- FindMarkers(merged.sc_rna_data, ident.1 = 11, min.pct = 0.25)
head(cluster12.markers, n = 5)

VlnPlot(merged.sc_rna_data, features = c(row.names(cluster12.markers)[1], row.names(cluster12.markers)[2]))


#### assigning cluster name

new.cluster.ids <- c("CD3G", "RPS27", "FGR", "NKG7", "CST3", "HLA-DGB1", "CD74", "RPL10", "GNLY","MALAT1","FTL","MRPS15")
names(new.cluster.ids) <- levels(merged.sc_rna_data)

library(scRNAseq)
library(scater)
Idents(object = merged.sc_rna_data) <- "sample"

freq_table <- prop.table(x = table(merged.sc_rna_data@active.ident, merged.sc_rna_data@active.ident),
                         margin = 2)
barplot(height = freq_table)

freq_table

gc()


# find markers for every cluster compared to all remaining cells
merged.sc_rna_data.markers <- FindAllMarkers(merged.sc_rna_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

x <- merged.sc_rna_data.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
FeaturePlot(merged.sc_rna_data, features = x$gene[1:4])

FeaturePlot(merged.sc_rna_data, features = x$gene[5:8])

p <- FeaturePlot(merged.sc_rna_data, features = c("TERF2IP", "TPT1", "EIF4E2", "CD247", "FCER1A", "LITAF", "LYZ", "FTL", "FTH1"), combine = FALSE)

p <- lapply(X = p, FUN = function(x) x + 
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5)) +
              theme(axis.title.x = element_text(size = 5)) +
              theme(axis.text.y = element_text(size = 5)) +
              theme(axis.text.x = element_text(size = 5)) +
              theme(legend.position = "none")  )

CombinePlots(plots = p)
### finding out top 10 markers
top10 <- merged.sc_rna_data.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

top10

p2 <- DoHeatmap(merged.sc_rna_data, features = top10$gene, group.bar.height = 0.01,size=3,combine = FALSE) 

p2 <- lapply(X = p2, FUN = function(x) x + 
               theme(plot.title = element_text(size = 8)) +
               theme(axis.title.y = element_text(size = 5)) +
               theme(axis.title.x = element_text(size = 5)) +
               theme(axis.text.y = element_text(size = 3)) +
               theme(legend.position = "none")  )

CombinePlots(plots = p2)


# Find differentially expressed features between CD27 and all other cells
Idents(object = merged.sc_rna_data) <- "seurat_clusters"

cluster1.de.markers <- FindMarkers(merged.sc_rna_data, ident.1 = "0", ident.2 = NULL, only.pos = TRUE)
# view results
head(cluster1.de.markers)


library(ComplexHeatmap)

heatmapdf <- cluster1.de.markers[1:25,]
row_ha = rowAnnotation("CD27" = anno_barplot(heatmapdf$pct.1),
                       "Others"= anno_barplot(heatmapdf$pct.2),
                       width = unit(10, "cm"))

heat_map <- Heatmap(heatmapdf$avg_log2FC,
               name = "Log2FC",
               cluster_rows = TRUE, 
               row_labels = rownames(heatmapdf), 
               right_annotation = row_ha,
               width = unit(2, "cm"))

heat_map

### Compared heatmap

sammple.markers <- FindMarkers(merged.sc_rna_data, ident.1 = "sample_1_control", ident.2 = "sample_2_test")
# view results
head(sammple.markers)

library(ComplexHeatmap)

heatmap_data_frame <- sample.markers[1:25,]
row_ha = rowAnnotation("sample_1_control"=anno_barplot(heatmap_data_frame$pct.1),
                       width=unit(10,"cm"))
heat_map_1 <- Heatmap(heatmap_data_frame$avg_log2FC,
                    name="log2FC",
                    cluster_rows=TRUE,
                    row_labels=rownames(heatmap_data_frame),
                    right_annotation=row_ha,
                    width=unit(1,"cm"))
heat_map_1

### Compared heatmap

sammple.markers <- FindMarkers(merged.sc_rna_data, ident.1 = "sample_2_test", ident.2 = "sample_1_control")
# view results
head(sammple.markers)

library(ComplexHeatmap)

heatmap_data_frame <- sample.markers[1:25,]
row_ha = rowAnnotation("sample_2_test"=anno_barplot(heatmap_data_frame$pct.1),
                       width=unit(10,"cm"))
heat_map_2 <- Heatmap(heatmap_data_frame$avg_log2FC,
                      name="log2FC",
                      cluster_rows=TRUE,
                      row_labels=rownames(heatmap_data_frame),
                      right_annotation=row_ha,
                      width=unit(1,"cm"))
heat_map_2

heat_map_1 + heat_map_2



sessionInfo()



