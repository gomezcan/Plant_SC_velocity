
library(ggplot2)
library(ggpubr)
library(dplyr)
library(Seurat)
library(patchwork)


# Load the example dataset
test.data <- Read10X(data.dir = "../3_Counts/Counts_SRR12046049_outs/filtered_feature_bc_matrix/")
star_test.data <- Read10X(data.dir = "../3_Counts/SRR12046049_star.Counts_GeneFull/filtered/")

# Initialize the Seurat object with the raw (non-normalized data).
test <- CreateSeuratObject(counts = test.data, project = "SRR12046049", min.cells = 3, min.features = 200)

test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^ATM")

# Visualize QC metrics as a violin plot

VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# % of UMIs on mt genes vs total UMIs by cell
plot1 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt")

# sub 
dim(test)
test <- subset(test, subset = nFeature_RNA > 500  & percent.mt < 5) # & nFeature_RNA < 2500
dim(test)

plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1 + plot2

# Normalization
test <- SCTransform(test, vst.flavor="v2")

# Variable features
test <- FindVariableFeatures(test, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(test), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(test)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

########################
## Scaling the data
########################
all.genes <- rownames(test)
test <- ScaleData(test, features = all.genes)

# Dimansional reduction
test <- RunPCA(test, features = VariableFeatures(object = test))

DimPlot(test, reduction = "pca")
DimHeatmap(test, dims = 3, cells = 500, balanced = TRUE)

ElbowPlot(test)

# Cluster the cells
test <- FindNeighbors(test, dims = 1:15)
test <- FindClusters(test, resolution = 0.5)

# Non-linear dimensional reduction (UMAP/tSNE)
test <- RunUMAP(test, dims = 1:15)

DimPlot(test, reduction = "umap")

########################################################
## find markers for every cluster compared to all 
## remaining cells, report only the positive ones
## Parameter notes:
# min.pct: only test genes that are detected in a 
#   minimum fraction of min.pct cells in either of 
#   the two populations
########################################################


test.markers <- FindAllMarkers(test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

test.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# test
VlnPlot(test, features = c("AT1G50060", "AT5G54370", 
                           "AT1G28290", "AT1G54000", 
                           "AT5G64100", "AT3G09260"))
