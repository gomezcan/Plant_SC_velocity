################
## LIBRARIES
###############
suppressMessages(library(Seurat))
library(Matrix)
library(CytoTRACE)

## CytoTRACE notes
# List objects: 
## 1. Numeric values for CytoTRACE (0:more differentiated, 1:less differentiated), 
## 2. Ranked CytoTRACE cells 
## 3. GCS, 
## 4. Gene counts, 
## 5. Numeric vector of the Pearson correlation between each gene and CytoTRACE, 
## 6. Numeric vector of the Pearson correlation between each gene and gene counts, 
## 7. IDs of filtered cells
## 8. Normalized gene expression table

rm(list=ls())

#######################################################
###############        Functions        ###############
#######################################################

chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}


Subset_ByCellsType <- function(Celltype){
  
  # cells to keep 
  keep_cells  <- colnames(sobj)[grepl(Celltype, sobj@meta.data$time.celltype.anno)]
  print(length(keep_cells))
  
  sobj2 <- subset(sobj, cells=keep_cells)
  
  # second filter 
  keep_cells  <- colnames(sobj2)[grepl(Celltype, sobj2@meta.data$predicted.celltype.anno)]
  print(length(keep_cells))
  
  if(length(keep_cells) > 100){
    
    sobj2 <- subset(sobj, cells=keep_cells)
    
    SampleName <- gsub('_', '', SampleName)
    
    filesbase = paste0("InputFilesScVelo/", SampleName, "_", Celltype)
    filesbase <- gsub(' ', '_', filesbase)
    filesbase <- gsub('Putative_', '', filesbase)
    
    # write metadata
    write.csv(sobj2@meta.data, file=paste0(filesbase, "_metadata.csv"), quote=F, row.names=F)
    
    # write expression counts matrix
    counts_matrix <- GetAssayData(sobj2, assay='SCT', slot='counts')
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
  
}

Subset_ByTissues <- function(tissue, tissue_name){
  tissues <- as.character(unlist(tissue))
  tissues <- paste0(tissue[1:length(tissue)], collapse = "|")
  print(tissues)
  
  # cells to keep 
  keep_cells  <- colnames(sobj)[grepl(tissues, sobj@meta.data$time.celltype.anno)]
  print(length(keep_cells))
  
  
  
  if(length(keep_cells) > 100){
    
    sobj2 <- subset(sobj, cells=keep_cells)
    # write metadata
    # SampleName <- gsub('_', '', SampleName)
    
    filesbase = paste0("InputFilesScVelo/", SampleName, "_", tissue_name)
    filesbase <- gsub(' ', '_', filesbase)
    filesbase <- gsub('Putative_', '', filesbase)
    
    metaname = paste0(filesbase, "_metadata.csv")
    print(metaname)
    
    metatem <- as.data.frame(sobj2@meta.data)
  
    write.csv(metatem, file=metaname, quote=FALSE, row.names=FALSE)
    
    # write expression counts matrix
    counts_matrix <- GetAssayData(sobj2, assay='SCT', slot='counts')
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
  
}

#######################################################

args = commandArgs(trailingOnly=TRUE)

# Load seurat obj
sobj = readRDS(args[1]) # path to rds file 
sobj <- subset(sobj, subset=nCount_RNA>=500) # total number of molecules detected within a cell

# Define Sample name
SampleName = gsub(".rds", "", chop(args[1], "[/]", 2)) 

# Save metadata table:
barcode <- colnames(sobj)
barcode <- chop(barcode, "[_-]", 1) # remove -1
#barcode <- paste0(chop(SampleName, "[_-]", 1), "_", barcode) # add sample id


#  add barcode and UMPs to metadata 
sobj$barcode <- barcode
sobj$RefUMAP_1 <- sobj@reductions$ref.umap@cell.embeddings[,1]
sobj$RefUMAP_2 <- sobj@reductions$ref.umap@cell.embeddings[,2]
sobj$UMAP_1 <- sobj@reductions$umap@cell.embeddings[,1]
sobj$UMAP_2 <- sobj@reductions$umap@cell.embeddings[,2]

ResultsCT <- GetAssayData(sobj, assay='SCT', slot='counts')
# remove -1 from cell names
colnames(ResultsCT) <- chop(colnames(ResultsCT), "[-]", 1)
# CytoTRACE job
ResultsCT <- CytoTRACE(as.matrix(ResultsCT),  ncores = 50)

# save CytoTRACE results
CytoTRACEsname = paste0("CytoTRACE_results/CytoTRACE_", SampleName, ".rds")
saveRDS(ResultsCT, file = CytoTRACEsname)

# Add CytoTRACE and CytoTRACErank to metadata
sobj$CytoTRACE <- ResultsCT$CytoTRACE[match(as.character(sobj$barcode), names(ResultsCT$CytoTRACE))] 

# Save full dataset
filesbase = paste0("InputFilesScVelo/", SampleName)

# Metadata
metatem <- as.data.frame(sobj@meta.data)
write.csv(metatem, file=paste0(filesbase, "_metadata.csv"), quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- GetAssayData(sobj, assay='SCT', slot='counts')
writeMM(counts_matrix, file=paste0(filesbase, "_counts.mtx"))

# write dimensional reduction matrix, in this example case PCA matrix
write.csv(sobj@reductions$pca@cell.embeddings,
          file=paste0(filesbase, "_pca.csv"), quote=F, row.names=F)

# Write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)), file=paste0(filesbase, "_gene_names.csv"),
  quote=F,row.names=F,col.names=F
)

##################################################
##### Define groups of cells types to write  #####
##################################################

Samplelist <- unique(sobj$time.celltype.anno)
Samplelist <- unique(chop(Samplelist, "[_-]", 2))

lapply(Samplelist, Subset_ByCellsType)

##################################################


##################################################
###  Define groups of tissues types to write  ####
##################################################

## Tissues definition
## Tissues definition
Stele <- c("Procambium", "Pericycle", "Phloem", "Xylem", "Putative Quiescent Center", "Stem Cell Niche")
GroundTissue <- c("Endodermis", "Cortex", "Putative Quiescent Center", "Stem Cell Niche")
RootCap <- c("Lateral Root Cap", "Columella", "Putative Quiescent Center", "Stem Cell Niche")
Epidermis <- c("Atrichoblast", "Trichoblast", "Putative Quiescent Center", "Stem Cell Niche")

Tissues <- list(Stele, GroundTissue, RootCap, Epidermis)
names(Tissues) <- c("Stele", "GroundTissue", "RootCap", "Epidermis")

Subset_ByTissues(Tissues$Stele, "Stele")
Subset_ByTissues(Tissues$GroundTissue, "GroundTissue")
Subset_ByTissues(Tissues$RootCap, "RootCap")
Subset_ByTissues(Tissues$Epidermis, "Epidermis")


##################################################
