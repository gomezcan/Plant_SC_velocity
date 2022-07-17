################
## LIBRARIES
###############
suppressMessages(library(Seurat))
library(Matrix)

chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

args = commandArgs(trailingOnly=TRUE)

# Load seurat obj
sobj = readRDS(args[1]) # path to rds file 
sobj <- subset(sobj, subset=nCount_RNA>=500) # total number of molecules detected within a cell

# Define Sample name
SampleName = gsub(".rds", "", chop(args[1], "[/]", 2)) 

# Save metadata table:
barcode <- colnames(sobj)
barcode <- chop(barcode, "[-]", 1) # remove -1
barcode <- paste0(chop(SampleName, "[_-]", 1), "_", barcode) # add sample id
barcode[1:5]

#  add barcode and UMPs to metadata 
sobj$barcode <- barcode
sobj$RefUMAP_1 <- sobj@reductions$ref.umap@cell.embeddings[,1]
sobj$RefUMAP_2 <- sobj@reductions$ref.umap@cell.embeddings[,2]
sobj$UMAP_1 <- sobj@reductions$umap@cell.embeddings[,1]
sobj$UMAP_2 <- sobj@reductions$umap@cell.embeddings[,2]

# define groups of cells types to 
Samplelist <- unique(sobj$time.celltype.anno)
Samplelist <- unique(chop(Samplelist, "[_-]", 2))


SaveSubset_files <- function(Celltype){
  
  # cell to keep 
  keep_cells  <- colnames(sobj)[grepl(Celltype, sobj@meta.data$time.celltype.anno)]
  print(length(keep_cells))
  
  sobj2 <- subset(sobj, cells=keep_cells)
  
  # write metadata
  SampleName <- gsub('_', '', SampleName)
  
  filesbase = paste0("InputFilesScVelo/", SampleName, "_", Celltype)
  filesbase <- gsub(' ', '_', filesbase)
  filesbase <- gsub('Putative_', '', filesbase)
  
  write.csv(sobj2@meta.data, file=paste0(filesbase, "_metadata.csv"), quote=F, row.names=F)
  
  # write expression counts matrix
  counts_matrix <- GetAssayData(sobj2, assay='RNA', slot='counts')
  writeMM(counts_matrix, file=paste0(filesbase, "_counts.mtx"))
  
  # write dimensional reduction matrix, in this example case PCA matrix
  write.csv(sobj2@reductions$pca@cell.embeddings, 
            file=paste0(filesbase, "_pca.csv"), quote=F, row.names=F)
  
  # write gene names
  write.table(
    data.frame('gene'=rownames(counts_matrix)), file=paste0(filesbase, "_gene_names.csv"),
    quote=F,row.names=F,col.names=F
  )
  
}

lapply(Samplelist, SaveSubset_files)
