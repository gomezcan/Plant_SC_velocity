################
## LIBRARIES
###############
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
library(tibble)


args = commandArgs(trailingOnly=TRUE)

cellrange_path = args[1] # path to filtered_feature_bc_matrix 
SampleName = args[2]     # sample name   

#cellrange_path = "../4_Velocyto/matrix_10S10/"
#SampleName ="sc_10S10"

# Load the example cellrange path
Sample.data <- Read10X(data.dir = cellrange_path)

# Initialize the Seurat object with the raw (non-normalized data).
Sample <- CreateSeuratObject(counts = Sample.data,  project = SampleName, 
                           min.cells = 3, min.features = 200)

# Percentage of mito reads
# Sample[["percent.mt"]] <- PercentageFeatureSet(Sample.data, pattern = "^ATM")

# Read Seurat Atlas object
atlas = readRDS("../3_Counts/GSE152766_Root_Atlas.rds")

Sample <- NormalizeData(Sample)
Sample <- FindVariableFeatures(Sample)
Sample <- ScaleData(Sample)
print(".. scaling done ... ")

# find anchors
anchors <- FindTransferAnchors(reference = atlas, query = Sample)

# save anchor objs
anchorsname = paste0("Seurat_objs/anchors_", SampleName, ".rds")
saveRDS(anchors, file = anchorsname)
print(".. anchors done ... ")

# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = atlas$timezone.ID, store.weights=FALSE)
predictions_2 <- TransferData(anchorset = anchors, refdata = atlas$celltype.anno, store.weights=FALSE)
predictions_3 <- TransferData(anchorset = anchors, refdata = atlas$time.celltype.anno, store.weights=FALSE)
predictions_4 <- TransferData(anchorset = anchors, refdata = atlas$consensus.time.group, store.weights=FALSE)
predictions_5 <- TransferData(anchorset = anchors, refdata = atlas$seurat_clusters, store.weights=FALSE)

# Name predicted.ID with label name
colnames(predictions)[1] <- "timezone.ID"
colnames(predictions_2)[1] <- "celltype.anno"
colnames(predictions_3)[1] <- "time.celltype.anno"
colnames(predictions_4)[1] <- "consensus.time.group"
colnames(predictions_5)[1] <- "Atlas_seurat_clusters"

# Combine labels
Labels <- left_join(rownames_to_column(predictions),  rownames_to_column(predictions_2), by=c("rowname"))
Labels <- left_join(Labels, rownames_to_column(predictions_3), by=c("rowname"))
Labels <- left_join(Labels, rownames_to_column(predictions_4), by=c("rowname"))
Labels <- left_join(Labels, rownames_to_column(predictions_5), by=c("rowname"))

# formatted "Labels" keeping only categorical prediction
col_lables <- c("rowname", "timezone.ID", "celltype.anno", 
                "time.celltype.anno", "consensus.time.group", 
                "Atlas_seurat_clusters")

Labels <- Labels[, colnames(Labels) %in% col_lables]
row.names(Labels) <- Labels$rowname
Labels <- Labels[,-c(1)] # remove bc

Sample <- AddMetaData(object = Sample, metadata = Labels)

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

filen=paste0("Seurat_objs/", SampleName, ".rds")
saveRDS(Sample, file = filen)

# Cluster the cells
# Sample <- FindNeighbors(Sample, dims = 1:15)
# Sample <- FindClusters(Sample, resolution = 0.5)
