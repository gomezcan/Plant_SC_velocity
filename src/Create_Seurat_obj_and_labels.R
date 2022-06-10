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

# Normalization v2
Sample <- SCTransform(Sample, vst.flavor="v2") # vst.flavor="v2"

# Variable features
Sample <- FindVariableFeatures(Sample, selection.method = "vst")
# scaling after SCTransform
Sample <- ScaleData(Sample)
print(".. Second scaling done ... ")

# Dimensional reduction
Sample <- RunPCA(Sample, features = VariableFeatures(object = Sample))

# Non-linear dimensional reduction (UMAP/tSNE)
Sample <- RunUMAP(Sample,  return.model = TRUE, dims = 1:20,  n_neighbors = 30, min_dist = 0.3, 
                  umap.method = "umap-learn",  metric = "correlation")

# Read Seurat Atlas object
atlas = readRDS("GSE152766_Root_Atlas.vst2.rds")
# atlas <- SCTransform(atlas, vst.flavor="v2") 
# atlas <- RunUMAP(atlas, return.model = TRUE, dims = 1:50,  n_neighbors = 30, min_dist = 0.3, 
#                  umap.method = "umap-learn",  metric = "correlation")
### saveRDS(atlas, file="../4_Velocyto/GSE152766_Root_Atlas.vst2.rds")

# Find anchors
anchors <- FindTransferAnchors(reference = atlas, query = Sample,
                               reference.reduction='pca', normalization.method='SCT')

# save anchor objs
anchorsname = paste0("Seurat_objs/anchors_", SampleNamexx, ".rds")
saveRDS(anchors, file = anchorsname)
print(".. anchors done ... ")

Sample <- MapQuery(anchorset = anchors, 
                   reference = atlas, 
                   query = Sample, 
                   refdata = list(celltype.anno = "celltype.anno"),  
                   reference.reduction = "pca", 
                   reduction.model = "umap")


# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = atlas$time.anno, store.weights=FALSE)
predictions_2 <- TransferData(anchorset = anchors, refdata = atlas$time.celltype.anno, store.weights=FALSE)
predictions_3 <- TransferData(anchorset = anchors, refdata = atlas$consensus.time.group, store.weights=FALSE)
predictions_4 <- TransferData(anchorset = anchors, refdata = atlas$seurat_clusters, store.weights=FALSE)

# Name predicted.ID with label name
colnames(predictions)[1] <- "time.anno"
colnames(predictions_2)[1] <- "time.celltype.anno"
colnames(predictions_3)[1] <- "consensus.time.group"
colnames(predictions_4)[1] <- "seurat_clusters"

# Combine labels
Labels <- left_join(rownames_to_column(predictions),  rownames_to_column(predictions_2), by=c("rowname"))
Labels <- left_join(Labels, rownames_to_column(predictions_3), by=c("rowname"))
Labels <- left_join(Labels, rownames_to_column(predictions_4), by=c("rowname"))

# formatted "Labels" keeping only categorical prediction
col_lables <- c("rowname", "time.anno",  "time.celltype.anno", 
                "consensus.time.group", "Atlas_seurat_clusters")

Labels <- Labels[, colnames(Labels) %in% col_lables]
row.names(Labels) <- Labels$rowname
Labels <- Labels[,-c(1)] # remove bc

## transfer UMAP coordenates
# DimPlot(Sample, reduction = "umap")

Sample <- AddMetaData(object = Sample, metadata = Labels)


# set cell type colors
colors_celltype <- c("aquamarine4", "deepskyblue", "forestgreen", 
                     "greenyellow", "lightpink1", "darkorange4",
                     "darkorange1", "blue", "brown1", 
                     "goldenrod1", "mediumseagreen","sandybrown",
                     "violet", "purple")
names(colors_celltype) <- as.character(unique(atlas$celltype.anno))

# set developmental colors
colors_timeAnno <- c("darkseagreen1","seagreen3", "firebrick3", 
                     "lightblue2", "forestgreen", "royalblue2",
                     "rosybrown2")
names(colors_timeAnno) <- as.character(unique(atlas$time.anno))

p1 <- DimPlot(atlas, reduction = "umap", group.by = "celltype.anno", 
              label = TRUE, label.size = 3, repel = TRUE) + 
        NoLegend() + ggtitle("Atlas Reference") + 
        scale_color_manual(values=colors_celltype)

p2 <- DimPlot(Sample, reduction = "ref.umap", group.by = "predicted.celltype.anno", 
              label = TRUE, label.size = 3, repel = TRUE, raster=FALSE) + 
      NoLegend() + ggtitle("Query transferred \nUMAP coords & Cell type") +
      scale_color_manual(values=colors_celltype)

p3 <- DimPlot(atlas, reduction = "umap", group.by = "time.anno", 
              label = TRUE, label.size = 3, repel = TRUE) + 
  NoLegend() + ggtitle("Atlas Reference") +
  scale_color_manual(values=colors_timeAnno)

p4 <- DimPlot(Sample, reduction = "ref.umap", group.by = "time.anno", 
              label = TRUE, label.size = 3, repel = TRUE, raster=FALSE) + 
  NoLegend() +ggtitle("Query transferred \nUMAP coords & time.anno") +
  scale_color_manual(values=colors_timeAnno)

Final_plot <- (p1 + p2)/(p3 + p4)

 # 10x8
print(".. ploting atlas and labels transterred .. ")

pdfname = paste0("Plot_atlas_And_", SampleName, ".pdf")

pdf(pdfname, bg = "white", width = 8, height = 9)
print(Final_plot)
dev.off() 

filen=paste0("Seurat_objs/", SampleName, ".rds")
saveRDS(Sample, file = filen)

# Cluster the cells
# Sample <- FindNeighbors(Sample, dims = 1:15)
# Sample <- FindClusters(Sample, resolution = 0.5)
## Plots exploring PCA results ##
# DimPlot(Sample, reduction = "pca")
# DimHeatmap(Sample, dims = 1, cells = 500, balanced = TRUE)
# ElbowPlot(Sample)
