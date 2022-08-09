#!/usr/bin/env python

import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

import scvelo as scv
import cellrank as cr
import seaborn as sns
import sys

####################################################
### Description:
## read metadata files create from seurat object
## with "Set_Input_scVelo.R" to calcuate scVelo velocities
## Main outputs:
## 1. Cell metadata incluting their velocities and pseudo-time: 'MetaCell_Sample.txt'
## 2. Gene metadata incluting their kinetic parameters:         'MetaGene_Sample.txt'
## 3. adata.h5ad file:                                          'adata_Sample.h5ad'     
## 4. Plot velocity field:                                      'Sample_umap.pdf'
## 5. Embeddings and velocities coordenates:                    'Vel_Embeddings_Sample.txt'
## 6. Gene's velocity vector by cell:                           'VelByCell_Sample.csv'
## 7. Gene's expression vector by cell:                         'ExpressionByCell_Sample.csv'
## 8. PAGA-cluster interactions:                                'PAGAinteractions_Sample.csv' 
####################################################


Sample = sys.argv[1]
namebase = "InputFilesScVelo/"+Sample

#Sample = "sc10S10_Endodermis"
#namebase = "InputFilesScVelo/"+'sc10S10_Endodermis'

# load sparse matrix:
X = io.mmread(namebase+"_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv(namebase+"_metadata.csv")

# load gene names:
with open(namebase+"_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv(namebase+"_pca.csv")
pca.index = adata.obs.index

# set pca
adata.obsm['X_pca'] = pca.to_numpy()
# set pca and umap
adata.obsm['X_Refumap'] = np.vstack((adata.obs['RefUMAP_1'].to_numpy(), adata.obs['RefUMAP_2'].to_numpy())).T
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# color file
colorfile = pd.read_csv('colorfile_python')
colorfile = colorfile.set_index('cluster').to_dict()["color"]

# adata = sc.read_h5ad('my_data.h5ad')

## load loom files for spliced/unspliced matrices for each sample:
namebase_loom = namebase.replace("InputFilesScVelo/", "").split("_")[0] # define loom file name
ldata1 = scv.read('Counts_'+namebase_loom+'.loom', cache=True)

# Clean barcode names
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1]  for bc in barcodes] # to remove las value

# make variable names unique
ldata1.obs.index = barcodes
ldata1.var_names_make_unique()

# merge matrices into the original adata object
adata_f = scv.utils.merge(adata, ldata1)

# Compute a neighborhood graph of observation
sc.pp.neighbors(adata_f, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata_f)

# Louvain clustering
sc.tl.louvain(adata_f, resolution=1.0)

# PAGA graph
adata_f.uns['neighbors']['distances'] = adata_f.obsp['distances']
adata_f.uns['neighbors']['connectivities'] = adata_f.obsp['connectivities']
sc.tl.paga(adata_f, groups='louvain')

# DF with connectivities
dfpaga = scv.get_df(adata_f, 'paga/connectivities', precision=2).T
dfpaga["LouvainNode1"] = dfpaga.index
dfpaga = dfpaga.melt(id_vars=['LouvainNode1'], ignore_index=True)

dfpaga = dfpaga[dfpaga.value > 0]
dfpaga = dfpaga.rename(columns={"variable":"LouvainNode2"})


# pre-process adata for scVelo

# Minimum number of counts (both unspliced and spliced) required for a gene.
# min_shared_cells : int, optional (default: None) Minimum number of cells 
# required to be expressed (both unspliced and spliced)

scv.pp.filter_and_normalize(adata_f, min_shared_counts=20, n_top_genes=2000)

scv.pp.moments(adata_f, n_pcs=30, n_neighbors=10)

scv.tl.recover_dynamics(adata_f, n_jobs=50)

## Compute velocity
scv.tl.velocity(adata_f, mode='dynamical')
scv.tl.velocity_graph(adata_f, n_jobs=50)

# Speed and coherence
scv.tl.velocity_confidence(adata_f)

# Velocity_pseudotime
scv.tl.velocity_pseudotime(adata_f)

# Latent time
scv.tl.latent_time(adata_f)

# Define df with gene counts and Ms values
gid = adata_f.var.Accession
gid = gid.reset_index(drop=True)

spliced = pd.DataFrame(adata_f.layers['spliced'].sum(axis=0).transpose(), 
                      columns=['spliced'])

unspliced = pd.DataFrame(adata_f.layers['unspliced'].sum(axis=0).transpose(), 
                      columns=['unspliced'])

Ms = pd.DataFrame(adata_f.layers['Ms'].sum(axis=0).transpose(), 
                      columns=['Ms'])

Mu = pd.DataFrame(adata_f.layers['Mu'].sum(axis=0).transpose(), 
                      columns=['Mu'])

Counts = pd.concat([gid, spliced, unspliced, Ms, Mu], axis=1)


###############################################
#######          save results           #######
###############################################

NoSave = ["predicted.celltype.anno.score", 
          "consensus.time.group",
         "predicted.celltype.anno",
          "time.anno"]

# 1. MetaCell file
adata_f.obs.drop(NoSave, axis=1).to_csv('scVeloResults/MetaCell_'+Sample+'.txt',
                                        sep='\t', 
                                        header=True,
                                        index=False 
                                       )
# 2. MetaGene file
# left join count to Gene Metatable
MetaGene = adata_f.var
MetaGene = MetaGene.merge(Counts, on='Accession', how='left')

MetaGene.rename({'Accession': 'GeneID'}, axis=1, inplace=True)

MetaGene.to_csv('scVeloResults/MetaGene_'+Sample+'.txt',
                                        sep='\t', 
                                        header=True,
                                        index=False 
                                       )

# 3. adata obj
adata.write('scVeloResults/adata_'+Sample+'.h5ad', compression='gzip')


# 4. Velocity field plots

# save velocity field coloring pseudotime and latent time
# with ref and predicted umaps as embbeddinds 

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', 
                               dpi=100, fontsize=8, frameon=False)

keys = ['velocity_pseudotime', "latent_time"]


scv.pl.velocity_embedding_grid(adata_f, 
                               basis="Refumap",
                               color=keys, 
                               title=keys, 
                               scale=0.35,
                               cmap='gnuplot', 
                               figsize= (5,5), 
                              save=Sample+'_Refumap.pdf')

scv.pl.velocity_embedding_grid(adata_f, 
                               basis="umap",
                               color=keys, 
                               title=keys, 
                               scale=0.35,
                               cmap='gnuplot', 
                               figsize= (5,5), 
                               save=Sample+'_umap.pdf',
                               fontsize=14)

# 5. Create and save DF with embbeddinds
df1 = pd.DataFrame(adata_f.obsm["X_Refumap"], columns=['x0_ref', 'y0_ref'])
df2 = pd.DataFrame(adata_f.obsm["velocity_Refumap"], columns=['x1_ref', 'y1_ref'])

df3 = pd.DataFrame(adata_f.obsm["X_umap"], columns=['x0', 'y0'])
df4 = pd.DataFrame(adata_f.obsm["velocity_umap"], columns=['x1', 'y1'])

Vel_Embeddings = pd.concat([df1,df2,df3, df4], axis=1)

Vel_Embeddings["Barcode"] = adata_f.obs.index

Vel_Embeddings.to_csv('scVeloResults/Vel_embeddings_'+Sample+'.txt',
                                        sep='\t', 
                                        header=True,
                                        index=False 
                                       )

# 6. Save velocity by cell
velocity = pd.DataFrame(adata_f.layers["velocity"], 
                            columns=adata_f.var.index)

mask_C = velocity.columns[(velocity.sum(axis=0) != 0)]
velocity = velocity[mask_C]

velocity["barcode"] = adata_f.obs.barcode.tolist()
velocity.to_csv('scVeloResults/VelByCell_'+Sample+'.csv',
                header=True, index=False)

# 7. Save expression by cell
Expression_Ms = pd.DataFrame(adata_f.layers["Ms"], 
                        columns=adata_f.var.index)

Expression_Mu = pd.DataFrame(adata_f.layers["Mu"], 
                        columns=adata_f.var.index)

Expression =  Expression_Ms + Expression_Mu

Expression["barcode"] = adata_f.obs.barcode.tolist()
Expression.to_csv('scVeloResults/ExpressionByCell_'+Sample+'.csv',
                header=True, index=False)

# 8. PAGA-cluster interactions
dfpaga.to_csv('scVeloResults/PAGAinteractions_'+Sample+'.csv',
                header=True, index=False)
