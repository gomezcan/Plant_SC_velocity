################
## LIBRARIES
###############
suppressMessages(library(Seurat))
library(Matrix)

args = commandArgs(trailingOnly=TRUE)

# Load seurat obj
sobj = readRDS(args[1]) # path to rds file 
# sobj=readRDS("Seurat_objs/sc_10S10.rds")
sobj <- subset(sobj, subset=nCount_RNA>=500)

SampleName=args[2] # Sample name

# save metadata table:
sobj$barcode <- colnames(sobj)
sobj$UMAP_1 <- sobj@reductions$ref.umap@cell.embeddings[,1]
sobj$UMAP_2 <- sobj@reductions$ref.umap@cell.embeddings[,2]
# sobj$UMAP_1 <- sobj@reductions$umap@cell.embeddings[,1]
# sobj$UMAP_2 <- sobj@reductions$umap@cell.embeddings[,2]


write.csv(sobj@meta.data, file=paste0(SampleName, "_metadata.csv"), quote=F, row.names=F)


# write expression counts matrix
counts_matrix <- GetAssayData(sobj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(SampleName, "_counts.mtx"))


# write dimensional reduction matrix, in this example case pca matrix
write.csv(sobj@reductions$pca@cell.embeddings, 
          file=paste0(SampleName, "_pca.csv"),
                      quote=F, row.names=F)

# write gene names

write.table(
  data.frame('gene'=rownames(counts_matrix)), file=paste0(SampleName, "_gene_names.csv"),
  quote=F,row.names=F,col.names=F
)
