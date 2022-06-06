#!/usr/bin/env Rscript

################
## LIBRARIES
###############
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library("velocyto.R"))
suppressMessages(library(argparser))

####################
## COMMAND LINE ARGUMENTS
#######################
p <- arg_parser("RNA Velocity Plots")
p <- add_argument(p, "objfile", help="an .RDS file containing a seurat object from one sample - OR - an .Rdata data file containing a seurat object called 'sobj'")
p <- add_argument(p, "loomfile", help="a .loom file produced by velocyto for the same sample")
p <- add_argument(p, "--clusters", help="name of meta.data field in seurat object containing cluster labels", default="seurat_clusters")
p <- add_argument(p, "--minumis", help="minimum number of UMIs to include a cell in the analysis", default=1000,type="numeric")
p <- add_argument(p, "--colorfile", help="a comma- or tab- separated file listing cluster names in column1 and color (for plotting) in column2", default=NULL)
p <- add_argument(p, "--cores", help="number of cores to use for parallelized functions", default=20, type="integer")
p <- add_argument(p, "--embedding", help="embedding to use for plotting", default="umap")
p <- add_argument(p, "--outfile", help=".Rdata filename to save computed results", default="<dirname objfile>/velocyto.Rdata")
p <- add_argument(p, "--plotfile", help=".pdf filename to save plots", default="<dirname objfile>/velocyto_plots.pdf")
p <- add_argument(p, "--spliced_cluster_avg", help="use only genes with an average of at least [n] spliced UMIs per cell in each cluster", default=0.1)
p <- add_argument(p, "--unspliced_cluster_avg", help="use only genes with an average of at least [n] unspliced UMIs per cell in each cluster", default=0.01)
p <- add_argument(p, "--gridsize", help="summarise per-cell velocity arrows into an [n]x[n] grid", default=40)


args <- parse_args(p)
if(grepl("<",args$outfile)){args$outfile=paste0(dirname(args$objfile),"/velocyto.Rdata")}
if(grepl("<",args$plotfile)){args$plotfile=paste0(dirname(args$objfile),"/velocyto_plots.pdf")}
if(!grepl(".Rdata$",args$outfile)){args$outfile=paste0(args$objfile,".Rdata")}
if(!grepl(".pdf$",args$plotfile)){args$plotfile=paste0(args$plotfile,".pdf")}


####################
## FUNCTIONS
###################

## chop a string by a separator and return specified field
chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}


#######################
## MAIN
####################

## load seurat object
cat("loading seurat object...")
if(grepl(".RDS$",args$objfile)){sobj=readRDS(args$objfile)}else{load(args$objfile)}
sobj<-subset(sobj,subset=nCount_RNA>=args$minumis)
## remove cells with clustersize=1... otherwise won't work with filtering step below
keep_cells=sobj@meta.data%>%as.data.frame%>%tibble::rownames_to_column("barcode")%>%
  mutate(clusters=sobj[[args$clusters]])%>%group_by(clusters)%>%mutate(cluster_size=n())%>%
  filter(cluster_size>1)%>%pull(barcode)
sobj<-subset(sobj,cells=keep_cells)
seurat_barcodes=chop(colnames(sobj),"[-_]",1)
cat(length(seurat_barcodes)," cell barcodes in seurat object\n")


## read in loom file
cat("loading velocyto loom file...")
ldat <- read.loom.matrices(file = args$loomfile)
if(sum(colnames(ldat$spliced)!=colnames(ldat$unspliced))>0){cat("Error; spliced and unspliced columns dont match.\n");quit(0)}
if(sum(rownames(ldat$spliced)!=rownames(ldat$unspliced))>0){cat("Error; spliced and unspliced rows dont match.\n");quit(0)}
loom_barcodes=gsub("x","",chop(colnames(ldat$spliced),":",2))
cat(length(loom_barcodes)," cell barcodes in velocyto object\n")

## extract count matrices
sobj$loom_ix=match(seurat_barcodes,loom_barcodes)
missing=sum(is.na(sobj$loom_ix))
if(missing>0){cat(missing," seurat cell barcodes missing from velocyto object\n");quit()}
spliced=ldat$spliced[,sobj$loom_ix];colnames(spliced)=seurat_barcodes
unspliced=ldat$unspliced[,sobj$loom_ix]; colnames(unspliced)=seurat_barcodes
amb=ldat$ambiguous[,sobj$loom_ix]; colnames(amb)=seurat_barcodes
rm(ldat)

## calculate overall splice rate
cat("overall splice rate: ",round(sum(spliced)/sum(spliced+unspliced)*100,digits=2),"%\n")

## get cluster labels for cells, and edit barcode names
clusters <- sobj[[args$clusters]]%>%tibble::rownames_to_column("barcode")%>%
    mutate(barcode=chop(barcode,"[_-]",1))%>%
    pull(2,name = barcode)

## get cell embedding coordinates
emb.plot <- Embeddings(object = sobj[[args$embedding]]); rownames(emb.plot)=seurat_barcodes
emb.pca <- Embeddings(object = sobj[["pca"]]); rownames(emb.pca)=seurat_barcodes

## calculate cell-cell distances in PCA
cell.dist <- as.dist(1-armaCor(t(emb.pca)))

## filter genes
cat(nrow(sobj)," genes in seurat object, ",nrow(spliced)," genes in velocyto object\n")
spliced <- filter.genes.by.cluster.expression(spliced,clusters,min.max.cluster.average = args$spliced_cluster_avg)
unspliced <- filter.genes.by.cluster.expression(unspliced,clusters,min.max.cluster.average = args$unspliced_cluster_avg)

use_genes=intersect(rownames(spliced),rownames(unspliced))
cat("filtered by cluster expression to ",nrow(spliced)," genes in spliced matrix, ",nrow(unspliced)," genes in unspliced matrix")
cat(" --> using ",length(use_genes)," total genes\n")

# calc velocity
cat("calculating velocity...\n")
fit.quantile <- 0.02
system.time({rvel.cd <- gene.relative.velocity.estimates(emat = spliced,nmat = unspliced, deltaT=1,kCells=50,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores=args$cores)})

## define per-cell velocity arrows
cat("drawing velocity plots...\n")
pdf(args$plotfile,width=10,height=5)
cat("saving plots to ",args$plotfile,"\n")

system.time({velocity_plot=suppressWarnings(
    show.velocity.on.embedding.cor(emb.plot,rvel.cd,n=100,scale='sqrt',cex=0.8,arrow.scale=5,arrow.lwd=.7,return.details=T,n.cores=1)
)})


### define velocity grid arrows
system.time({velocity_grid=show.velocity.on.embedding.cor(emb.plot,rvel.cd,n=200,scale='sqrt',cex=0.8,arrow.scale=5,arrow.lwd=1,
                               show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=args$gridsize,
                               do.par=F,cell.border.alpha = 0.1,cc = velocity_plot$cc,return.details=T,n.cores=1)
            })

### add arrow coords to seurat metadata
cat("adding arrow info to cell metadata...\n")
metadata <- sobj@meta.data %>% tibble::rownames_to_column("barcode") %>%
    mutate(barcode=chop(barcode,"[_-]",1)) %>%
    merge(velocity_grid$arrows%>%as.data.frame%>%tibble::rownames_to_column("barcode") ,
      by="barcode")

### use ggplot to make colored plots
gg <- metadata %>% mutate(clusters=metadata[[args$clusters]]) %>% 
ggplot(aes(x=x0,y=y0)) + geom_point(aes(color=clusters),alpha=.8) +
theme_void() + theme(legend.key.size = unit(.4,"lines"))

gg1 <- gg + ggtitle(paste0("RNA velocity + ",args$embedding," embedding")) +
geom_segment(data=velocity_grid$arrows%>%as.data.frame(),
             aes(xend = x1 , yend = y1),size=.3,alpha=.5,color="darkslategray",arrow = arrow(length = unit(0.1,"cm")))


gg2 <- gg + ggtitle(paste0(args$gridsize,"x",args$gridsize," RNA velocity grid + ",args$embedding," embedding")) +
geom_segment(data=velocity_grid$garrows%>%as.data.frame(),
             aes(xend = x1 , yend = y1),size=.3,alpha=.8,color="black",arrow = arrow(length = unit(0.1,"cm")))

if(!is.null(args$colorfile)){
    colors=fread(args$colorfile)%>%setNames(c("cluster","color"))%>%pull(color,name=cluster)
    gg1<-gg1+scale_color_manual(values=colors)
    gg2<-gg2+scale_color_manual(values=colors)
}

gg1; gg2
dev.off()


### save
cat("saving data to ",args$outfile)
save(rvel.cd,emb.plot,velocity_plot,velocity_grid,metadata,file=args$outfile)
cat("finished!\n")