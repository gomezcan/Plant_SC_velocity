################
## LIBRARIES
###############
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
################

args = commandArgs(trailingOnly=TRUE)

################
# Input files
################
cellrange_path = args[1] # path to filtered_feature_bc_matrix 
# cellrange_path = "../4_Velocyto/matrix_10S10/"
SampleName = args[2]     # sample name   
# SampleName ="sc_10S10"

# Read Seurat Atlas object
atlas = readRDS("GSE152766_Root_Atlas.rds")
################
print(".. reading done ... ")

# Read cellrange file
Sample.data <- Read10X(data.dir = cellrange_path)


# Initialize the Seurat object with the raw (non-normalized data).
Sample <- CreateSeuratObject(counts = Sample.data, 
                           project = SampleName, 
                           min.cells = 3, min.features = 200)

# standar normalizacion and scaling
Sample <- NormalizeData(Sample)
Sample <- FindVariableFeatures(Sample)
Sample <- ScaleData(Sample)
print(".. scaling done ... ")

# Find anchors
anchors <- FindTransferAnchors(reference = atlas, query = Sample)
print(".. anchors done ... ")

# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = atlas$celltype.anno)
Sample <- AddMetaData(object = Sample, metadata = predictions)

print(".. transfer labels done ... ")

# Normalization v2
Sample <- SCTransform(Sample, vst.flavor="v2") # vst.flavor="v2"
# Variable features
Sample <- FindVariableFeatures(Sample, selection.method = "vst", nfeatures = 2000)

# scaling after SCTransform
Sample <- ScaleData(Sample)

print(".. Second scaling done ... ")

# Dimensional reduction
Sample <- RunPCA(Sample, features = VariableFeatures(object = Sample))

## Plots exploring PCA results ##
# DimPlot(Sample, reduction = "pca")
# DimHeatmap(Sample, dims = 1, cells = 500, balanced = TRUE)
# ElbowPlot(Sample)

# Non-linear dimensional reduction (UMAP/tSNE)
Sample <- RunUMAP(Sample, dims = 1:15)
#DimPlot(test, reduction = "umap", group.by = "predicted.id")

file=paste0("Seurat_objs/", SampleName, ".rds")
saveRDS(Sample, file = file)

# Cluster the cells
# Sample <- FindNeighbors(Sample, dims = 1:15)
# Sample <- FindClusters(Sample, resolution = 0.5)
