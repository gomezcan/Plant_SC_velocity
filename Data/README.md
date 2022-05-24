## General notes describing data source and modification

### Genomes 
#### 1. Arabidopsis
- 1.1 Sequence: TAIR10.1_genomic.fna 
- 1.2 Annotation: [TAIR11 may2022](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release)

### Annotation modifications
#### gtf filtering (only protein-related regions):

```
# tem file to match chr nomemclature
cat Araport11_GTF_genes_transposons.May2022.gtf | sed 's/Chr/chr/g' | sed 's/chrC/chrP/g' > chrTAIR11.May2022.gtf

# filter redundant coordenates 
cellranger mkgtf chrTAIR11.May2022.gtf TAIR11_protein_filtered.gtf  --attribute=gene_biotype:protein

# Cellranger genome index
cellranger mkref --genome=Index_TAIR11_CellRange --fasta=TAIR10.1_genomic.fna --genes=TAIR11_protein_filtered.gtf

# STAR genome index
STAR --runThreadN 100 --runMode genomeGenerate --genomeDir ./Index_TAIR11_STAR --genomeFastaFiles TAIR10.1_genomic.fna --sjdbGTFfile TAIR11_protein_filtered.gtf --sjdbOverhang 91

```

### Samples: 
#### Arabidopsis root atlas: GEO GSE152766
