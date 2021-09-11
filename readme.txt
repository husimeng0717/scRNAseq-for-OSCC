The code used in "Single-cell RNA sequencing reveals that TDO2+ myofibroblasts mediate immune suppression in oral squamous cell carcinoma" is available from the folder code.

Raw 10x Genomics sequencing data were processed using the CellRanger software v. 4.0.0, and the 10x human transcriptome GRCh38-2020-A was used as the reference.

Analysis of single cell RNA-seq data was done using R (v3.6.3) with publicly available packages. 
Dimensionality reduction and differential gene expression was performed using the Seurat (v3.2.2) package. 
Double cell scoring was performed using the scDblFinder (v.1.4.0) package. 
Eliminating batch effects was performed using the fastMNN algorithm. 
Plots were generated using the ggplot2 (v 3.3.2), pheatmap (v 1.0.12), ggpubr (v0.4.0) packages and Cytoscape (v3.8.2). 
Gene ontology analysis was performed using the Metascape web resource. 
