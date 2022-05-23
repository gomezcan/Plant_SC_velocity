## General notes describing data source and modification

### Genomes 
#### 1. Arabidopsis
- 1.1 Sequence: TAIR10.1_genomic.fna 
- 1.2 Annotation: [TAIR11 may2022](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release)

### Annotation modifications
#### gtf filtering (only protein-related regions):

```
- 1. cat Araport11_GTF_genes_transposons.May2022.gtf | grep -v 'ChrC' | grep -v 'ChrM' > TAIR11_GTF_No_C_M.May2022.gtf
- 2. cellranger mkgtf TAIR11_GTF_No_C_M.May2022.gtf TAIR11_protein_filtered.gtf  --attribute=gene_biotype:protein

```

### Samples: 
#### Arabidopsis root atlas: GEO GSE152766
