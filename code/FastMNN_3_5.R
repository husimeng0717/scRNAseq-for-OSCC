library(SeuratWrappers)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(devtools)
library(scDblFinder)
library(BiocParallel)
library(pheatmap)
library(scales)

setwd("E:/OSCC/integrate/3_5/FastMNN/OSCC")

cell_type_cols <- c(brewer.pal(9, "Set1"),
                    "#FF34B3","#CD5C5C","#BC8F8F","#20B2AA","#ADFF2F","#FFA500","#FF6A6A","#7FFFD4",
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1",
                    "#7CFC00","#708090","#8B008B","#00F5FF","#FF3030",
                    "#6AB5E4", "#DB3B65", "#E6B2E7", "#E549E6" ,"#645AA4" ,"#E8826D", "#D0B178" ,"#DFD2EB",
                    "#7176E1", "#BFEB80", "#B6E9D8", "#E24EBE", "#D5EB3F", "#477787", "#E2C5A5", "#CF7FB7", "#DEA746", "#892772", "#B362DD")
show_col(c("#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#A65628", "#999999",  "#F781BF",  "#FFFF33","#4DAF4A"))


barcoder <- function(data10x, trim="\\-(.*)",sample){
  colnames(data10x) <- gsub(trim, "", colnames(data10x))
  colnames(data10x) <- paste(colnames(data10x),sample,sep = "-")
  data10x
}

#########     3' data     #############
ChSF.Ca <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/ChSF/Ca/filtered_feature_bc_matrix/")
ChSF.Ca <- barcoder(ChSF.Ca,sample="ChSF.Ca")
ChSF.N <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/ChSF/N/filtered_feature_bc_matrix/")
ChSF.N <- barcoder(ChSF.N,sample="ChSF.N")

CZC.Ca <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/CZC/Ca/filtered_feature_bc_matrix/")
CZC.Ca <- barcoder(CZC.Ca,sample="CZC.Ca")
CZC.N <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/CZC/N/filtered_feature_bc_matrix/")
CZC.N <- barcoder(CZC.N,sample="CZC.N")

CZX.Ca <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/CZX/Ca/filtered_feature_bc_matrix/")
CZX.Ca <- barcoder(CZX.Ca,sample="CZX.Ca")
CZX.N <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/CZX/N/filtered_feature_bc_matrix/")
CZX.N <- barcoder(CZX.N,sample="CZX.N")

HNH.Ca <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/HNH/Ca/filtered_feature_bc_matrix/")
HNH.Ca <- barcoder(HNH.Ca,sample="HNH.Ca")
HNH.N <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/HNH/N/filtered_feature_bc_matrix/")
HNH.N <- barcoder(HNH.N,sample="HNH.N")

LiYX.OLK <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/LiYX/OLK/filtered_feature_bc_matrix/")
LiYX.OLK <- barcoder(LiYX.OLK,sample="LiYX.OLK")

PShY.GCa <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/PShY/GCa/filtered_feature_bc_matrix/")
PShY.GCa <- barcoder(PShY.GCa,sample="PShY.GCa")
PShY.GN <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/PShY/GN/filtered_feature_bc_matrix/")
PShY.GN <- barcoder(PShY.GN,sample="PShY.GN")

WBSh.GCa <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/WBSh/GCa/filtered_feature_bc_matrix/")
WBSh.GCa <- barcoder(WBSh.GCa,sample="WBSh.GCa")

LJH.BCa <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/LJH/BCa/filtered_feature_bc_matrix/")
LJH.BCa <- barcoder(LJH.BCa,sample="LJH.BCa")

ZhRF.OLK <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/3/ZhRF/OLK/filtered_feature_bc_matrix/")
ZhRF.OLK <- barcoder(ZhRF.OLK,sample="ZhRF.OLK")

#########     5'      #############
ChZhX.N <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/ChZhX/N/5GEX/filtered_feature_bc_matrix/")
ChZhX.N <- barcoder(ChZhX.N,sample="ChZhX.N")

HKJ.Ca <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/HKJ/Ca/5GEX/filtered_feature_bc_matrix/")
HKJ.Ca <- barcoder(HKJ.Ca,sample="HKJ.Ca")
HKJ.N <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/HKJ/N/5GEX/filtered_feature_bc_matrix/")
HKJ.N <- barcoder(HKJ.N,sample="HKJ.N")

LWH.Ca <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/LWH/Ca/5GEX/filtered_feature_bc_matrix/")
LWH.Ca <- barcoder(LWH.Ca,sample="LWH.Ca")
LWH.N <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/LWH/N/5GEX/filtered_feature_bc_matrix/")
LWH.N <- barcoder(LWH.N,sample="LWH.N")

LYX.Ca <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/LYX/Ca/5GEX/filtered_feature_bc_matrix/")
LYX.Ca <- barcoder(LYX.Ca,sample="LYX.Ca")

PHD.BCa <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/PHD/BCa/5GEX/filtered_feature_bc_matrix/")
PHD.BCa <- barcoder(PHD.BCa,sample="PHD.BCa")

WRH.GCa <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/WRH/GCa/5GEX/filtered_feature_bc_matrix/")
WRH.GCa <- barcoder(WRH.GCa,sample="WRH.GCa")

XZB.OLK <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/XZB/OLK/5GEX/filtered_feature_bc_matrix/")
XZB.OLK <- barcoder(XZB.OLK,sample="XZB.OLK")

ZhChH.Ca <- Read10X(data.dir = "E:/OSCC/cellranger-4_0_0 result/5/ZhChH/Ca/5GEX/filtered_feature_bc_matrix/")
ZhChH.Ca <- barcoder(ZhChH.Ca,sample="ZhChH.Ca")

#########     3'+5' meta data     #############
meta.data<-data.frame(barcode = c(colnames(ChSF.Ca),colnames(ChSF.N),colnames(CZC.Ca),colnames(CZC.N),
                                  colnames(CZX.Ca),colnames(CZX.N),colnames(HNH.Ca),colnames(HNH.N),
                                  colnames(LiYX.OLK),colnames(PShY.GCa),colnames(PShY.GN),colnames(WBSh.GCa),
                                  colnames(LJH.BCa),colnames(ZhRF.OLK),
                                  
                                  colnames(ChZhX.N),colnames(HKJ.Ca),colnames(HKJ.N),colnames(LWH.Ca),
                                  colnames(LWH.N),colnames(LYX.Ca),colnames(PHD.BCa),colnames(WRH.GCa),
                                  colnames(XZB.OLK),colnames(ZhChH.Ca)))
row.names(meta.data)<-meta.data$barcode

meta.data$Sample <- c(rep("ChSF.T.Ca", ncol(ChSF.Ca)), rep("ChSF.T.N", ncol(ChSF.N)), rep("CZC.T.Ca", ncol(CZC.Ca)), rep("CZC.T.N", ncol(CZC.N)),
                      rep("CZX.T.Ca", ncol(CZX.Ca)), rep("CZX.T.N", ncol(CZX.N)), rep("HNH.T.Ca", ncol(HNH.Ca)), rep("HNH.T.N", ncol(HNH.N)),
                      rep("LiYX.OLK", ncol(LiYX.OLK)), rep("PShY.G.Ca", ncol(PShY.GCa)), rep("PShY.G.N", ncol(PShY.GN)),
                      rep("WBSh.G.Ca", ncol(WBSh.GCa)),rep("LJH.B.Ca", ncol(LJH.BCa)),rep("ZhRF.OLK", ncol(ZhRF.OLK)),
                      
                      rep("ChZhX.T.N", ncol(ChZhX.N)), rep("HKJ.T.Ca", ncol(HKJ.Ca)), rep("HKJ.T.N", ncol(HKJ.N)),
                      rep("LWH.T.Ca", ncol(LWH.Ca)), rep("LWH.T.N", ncol(LWH.N)), rep("LYX.T.Ca", ncol(LYX.Ca)), rep("PHD.B.Ca", ncol(PHD.BCa)),
                      rep("WRH.G.Ca", ncol(WRH.GCa)), rep("XZB.OLK", ncol(XZB.OLK)), rep("ZhChH.T.Ca", ncol(ZhChH.Ca)))

meta.data$Patient <- c(rep("ChSF", ncol(ChSF.Ca)), rep("ChSF", ncol(ChSF.N)), rep("CZC", ncol(CZC.Ca)), rep("CZC", ncol(CZC.N)),
                       rep("CZX", ncol(CZX.Ca)), rep("CZX", ncol(CZX.N)), rep("HNH", ncol(HNH.Ca)), rep("HNH", ncol(HNH.N)),
                       rep("LiYX", ncol(LiYX.OLK)), rep("PShY", ncol(PShY.GCa)), rep("PShY", ncol(PShY.GN)),
                       rep("WBSh", ncol(WBSh.GCa)), rep("LJH", ncol(LJH.BCa)), rep("ZhRF", ncol(ZhRF.OLK)),
                       
                       rep("ChZhX", ncol(ChZhX.N)), rep("HKJ", ncol(HKJ.Ca)), rep("HKJ", ncol(HKJ.N)),
                       rep("LWH", ncol(LWH.Ca)), rep("LWH", ncol(LWH.N)), rep("LYX", ncol(LYX.Ca)), rep("PHD", ncol(PHD.BCa)),
                       rep("WRH", ncol(WRH.GCa)),rep("XZB", ncol(XZB.OLK)), rep("ZhChH", ncol(ZhChH.Ca)))

meta.data$Tissue <- c(rep("Cancer", ncol(ChSF.Ca)), rep("Normal", ncol(ChSF.N)), rep("Cancer", ncol(CZC.Ca)), rep("Normal", ncol(CZC.N)),
                      rep("Cancer", ncol(CZX.Ca)), rep("Normal", ncol(CZX.N)), rep("Cancer", ncol(HNH.Ca)), rep("Normal", ncol(HNH.N)),
                      rep("OLK", ncol(LiYX.OLK)), rep("Cancer", ncol(PShY.GCa)), rep("Normal", ncol(PShY.GN)),
                      rep("Cancer", ncol(WBSh.GCa)), rep("Cancer", ncol(LJH.BCa)),rep("OLK", ncol(ZhRF.OLK)),
                      
                      rep("Normal", ncol(ChZhX.N)), rep("Cancer", ncol(HKJ.Ca)), rep("Normal", ncol(HKJ.N)),
                      rep("Cancer", ncol(LWH.Ca)), rep("Normal", ncol(LWH.N)), rep("Cancer", ncol(LYX.Ca)), rep("Cancer", ncol(PHD.BCa)),
                      rep("Cancer", ncol(WRH.GCa)), rep("OLK", ncol(XZB.OLK)), rep("Cancer", ncol(ZhChH.Ca)))

meta.data$Site <- c(rep("Tongue", ncol(ChSF.Ca)), rep("Tongue", ncol(ChSF.N)), rep("Tongue", ncol(CZC.Ca)), rep("Tongue", ncol(CZC.N)),
                    rep("Tongue", ncol(CZX.Ca)), rep("Tongue", ncol(CZX.N)), rep("Tongue", ncol(HNH.Ca)), rep("Tongue", ncol(HNH.N)),
                    rep("Tongue", ncol(LiYX.OLK)), rep("Gingival", ncol(PShY.GCa)), rep("Gingival", ncol(PShY.GN)),
                    rep("Gingival", ncol(WBSh.GCa)), rep("Buccal", ncol(LJH.BCa)),rep("Tongue", ncol(ZhRF.OLK)),
                    
                    rep("Tongue", ncol(ChZhX.N)), rep("Tongue", ncol(HKJ.Ca)), rep("Tongue", ncol(HKJ.N)),
                    rep("Tongue", ncol(LWH.Ca)), rep("Tongue", ncol(LWH.N)), rep("Tongue", ncol(LYX.Ca)), rep("Buccal", ncol(PHD.BCa)),
                    rep("Gingival", ncol(WRH.GCa)), rep("Tongue", ncol(XZB.OLK)), rep("Tongue", ncol(ZhChH.Ca)))

meta.data$Seq.method <- c(rep("SC3P", sum(ncol(ChSF.Ca), ncol(ChSF.N), ncol(CZC.Ca), ncol(CZC.N), ncol(CZX.Ca), ncol(CZX.N), ncol(HNH.Ca), ncol(HNH.N), ncol(LiYX.OLK), ncol(PShY.GCa), ncol(PShY.GN), ncol(WBSh.GCa), ncol(LJH.BCa), ncol(ZhRF.OLK))),
                          
                          rep("SC5P", sum(ncol(ChZhX.N), ncol(HKJ.Ca), ncol(HKJ.N), ncol(LWH.Ca), ncol(LWH.N), ncol(LYX.Ca), ncol(PHD.BCa), ncol(WRH.GCa), ncol(XZB.OLK), ncol(ZhChH.Ca))))

save(meta.data,file = "E:/OSCC/cellranger-4_0_0 result/3_5_meta.data.Rda")
load("E:/OSCC/cellranger-4_0_0 result/3_5_meta.data.Rda")

##########   Initialize Seurat Object   ##########
OSCC.data<-cbind(ChSF.Ca, ChSF.N, CZC.Ca, CZC.N, CZX.Ca, CZX.N, HNH.Ca, HNH.N, LiYX.OLK, PShY.GCa, PShY.GN, WBSh.GCa, LJH.BCa,ZhRF.OLK,
                 
                 ChZhX.N, HKJ.Ca, HKJ.N, LWH.Ca, LWH.N, LYX.Ca, PHD.BCa, WRH.GCa, XZB.OLK, ZhChH.Ca)

save(OSCC.data,file = "E:/OSCC/cellranger-4_0_0 result/3_5_data.Rda")
load("E:/OSCC/cellranger-4_0_0 result/3_5_data.Rda")

OSCC <- CreateSeuratObject(counts = OSCC.data, project = "OSCC", meta.data = meta.data, min.cells = round(0.001*ncol(OSCC.data)))#round(0.001*ncol(OSCC.data))=153 

OSCC[["percent.mt"]] <- PercentageFeatureSet(OSCC, pattern = "^MT-")
OSCC[["percent.RPL"]] <- PercentageFeatureSet(OSCC, pattern = "^RPL")
OSCC[["percent.RPS"]] <- PercentageFeatureSet(OSCC, pattern = "^RPS")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
OSCC <- CellCycleScoring(OSCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

pdf(file="vlnplot_genes_umi_mt_unfilter.pdf", width = 22, height = 4)
VlnPlot(OSCC, features = c("nFeature_RNA"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("nCount_RNA"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.mt"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.RPL"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.RPS"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)

VlnPlot(OSCC, features = c("nFeature_RNA"), group.by = "Sample", pt.size = 0, y.max =1000)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("nCount_RNA"), group.by = "Sample", pt.size = 0, y.max =5000)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.mt"), group.by = "Sample", pt.size = 0, y.max =50)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.RPL"), group.by = "Sample", pt.size = 0, y.max =50)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.RPS"), group.by = "Sample", pt.size = 0, y.max =50)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
dev.off()

pdf(file="FeatureScatter_nCount_nFeatur.pdf", width = 8, height = 6)
FeatureScatter(OSCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

OSCC <- subset(OSCC, subset = nFeature_RNA >= 300 & nFeature_RNA <= 7500 & nCount_RNA <= 60000 & percent.mt <= 20 )

OSCC <- as.SingleCellExperiment(OSCC)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
OSCC <- scDblFinder(OSCC,nfeatures = 1000,samples="Sample",BPPARAM=SnowParam(8))
OSCC <- as.Seurat(OSCC)
table(OSCC@meta.data$scDblFinder.class)

OSCC <- subset(OSCC, subset = scDblFinder.class=="singlet")
save(OSCC,file = "OSCC.Rda")

pdf(file="vlnplot_genes_umi_mt_filter.pdf", width = 22, height = 4)
VlnPlot(OSCC, features = c("nFeature_RNA"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("nCount_RNA"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.mt"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.RPL"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
VlnPlot(OSCC, features = c("percent.RPS"), group.by = "Sample", pt.size = 0)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust =0.5))+labs(x =NULL)+NoLegend()+scale_fill_manual(values=cell_type_cols)
dev.off()

OSCC <- NormalizeData(OSCC)
OSCC <- FindVariableFeatures(OSCC,nfeatures = 3000)
OSCC <- ScaleData(OSCC, vars.to.regress = c("percent.mt","nCount_RNA"), features =row.names(OSCC),verbose = FALSE) 
OSCC <- RunFastMNN(object.list = SplitObject(OSCC, split.by = "Sample"),features = 3000)
OSCC <- RunUMAP(OSCC, reduction = "mnn", dims = 1:30, umap.method='umap-learn')
OSCC <- FindNeighbors(OSCC, reduction = "mnn", dims = 1:30)
OSCC <- FindClusters(OSCC,resolution = 1.5)
save(OSCC,file = "OSCC.Rda")

load("OSCC.Rda")

pdf(file="Umap_number.pdf", width = 9, height = 6)  
DimPlot(OSCC, reduction = "umap", cols = cell_type_cols,label = T)
dev.off()

pdf(file="Umap_split_clusters.pdf", width = 20, height = 50)  
DimPlot(OSCC, reduction = "umap", cols = cell_type_cols,label = T,split.by = "seurat_clusters",ncol = 5)
dev.off()

OSCC$Tissue<-factor(OSCC$Tissue,levels = c("Normal","OLK","Cancer"))

pdf(file="Umap_split_Tissue.pdf", width = 18, height = 6)
DimPlot(OSCC, reduction = "umap", group.by = "cellType",split.by = "Tissue",cols = cell_type_cols,ncol = 3)
DimPlot(OSCC, reduction = "umap", group.by = "cellType",split.by = "Tissue",cols = OSCC.cellType.cols,ncol = 3)
dev.off()

pdf(file="Umap_group_Tissue.pdf", width = 8, height = 6)
DimPlot(OSCC, reduction="umap", group.by = "Tissue",cols = c("#73A1C7", "#D67C9B", "#9D70C1"))
dev.off()

pdf(file="Umap_split_Site.pdf", width = 9, height = 3)
DimPlot(OSCC, reduction = "umap", split.by = "Site",cols = cell_type_cols,ncol = 4)+NoLegend()
dev.off()

pdf(file="Umap_group_Site.pdf", width = 8, height = 6)
DimPlot(OSCC, reduction="umap", group.by = "Site",cols = cell_type_cols)
dev.off()

pdf(file="Umap_split_Seq.method.pdf", width = 13, height = 6)
DimPlot(OSCC, reduction = "umap",group.by = "cellType", split.by = "Seq.method",cols = cell_type_cols)
DimPlot(OSCC, reduction = "umap",group.by = "cellType", split.by = "Seq.method",cols = OSCC.cellType.cols)
dev.off()

pdf(file="Umap_group_Seq.method.pdf", width = 8, height = 6)
DimPlot(OSCC, reduction = "umap", group.by = "Seq.method",cols = c("#4DAF4A", "#984EA3"))
dev.off()

pdf(file="Umap_split_Sample.pdf", width = 18, height = 18)
DimPlot(OSCC, reduction = "umap", split.by = "Sample",cols = cell_type_cols,ncol = 5)+NoLegend()
dev.off()

pdf(file="Umap_group_Sample.pdf", width = 9, height = 6)
DimPlot(OSCC, reduction = "umap", group.by = "Sample",cols = cell_type_cols)
dev.off()

pdf(file="Umap_group_Phase.pdf", width = 8, height = 6)
DimPlot(OSCC, reduction="umap", group.by = "Phase",cols = cell_type_cols)
dev.off()

#各celltype细胞数的条形图
a<-as.data.frame(table(OSCC@meta.data$cellType))
a$rate<-round(a$Freq/sum(a$Freq),4)
cellNumber<-data.frame(Cluster=a$Var1,Number=a$Freq)
pdf(file="bar_cellNumber_cellType.pdf", width = 8, height = 6)
ggplot(cellNumber,aes(x=Cluster,y=Number,fill=Cluster,))+
  geom_bar(stat="identity",width=0.8)+
  labs(x='Cell Type',y='Cell Numbers')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.ticks.y=element_blank())+
  coord_flip()+
  geom_text(mapping = aes(label = Number),size = 5, colour = 'black', vjust = 0.5, hjust = 0)+
  NoLegend()
dev.off()

#各celltype细胞数的条形图:按Tissue分组
a<-as.data.frame(table(OSCC@meta.data$cellType,OSCC@meta.data$Tissue))
colnames(a)<-c("Cluster","Tissue","Number")
a$Tissue<-factor(a$Tissue,levels = c("Normal","OLK","Cancer"))
pdf(file="bar_cellNumber_cellType_Tissue.pdf", width = 8, height = 10)
ggplot(a,aes(x=Cluster,y=Number,fill=Tissue))+
  geom_bar(stat="identity",position = "dodge",width=0.8)+
  labs(x='Cell Type',y='Cell Numbers')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.ticks.y=element_blank())+
  coord_flip()+
  geom_text(aes(label = Number),position=position_dodge(width = 0.8),size = 4, colour = 'black', vjust = 0.5, hjust = 0)
dev.off()

#Tissue
orig.ident<-OSCC@meta.data[,c("orig.ident","Tissue")]
colnames(orig.ident)<-c("cellType","Tissue")
orig.ident$cellType<-"All T cells"
cellType<-OSCC@meta.data[,c("cellType","Tissue")]
meta<-rbind(orig.ident,cellType)
pdf(file="bar_percent_cellType_group_Tissue.pdf", width = 10, height = 6)
ggplot(cellType,aes(x=Tissue,fill=cellType))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='cellType',y='Cell proportion')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16),axis.ticks.y=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

cellType<-data.frame(table(OSCC@meta.data$cellType))
Tissue<-data.frame(table(OSCC@meta.data$Tissue))
cellType_Tissue<-data.frame(table(OSCC@meta.data$cellType,OSCC@meta.data$Tissue))
colnames(cellType_Tissue)<-c("cellType","Tissue","Num_cellType_split_Tissue")
Cancer<-cellType_Tissue[which(cellType_Tissue$Tissue=="Cancer"),]
Cancer[match(Cancer$cellType,cellType$Var1),"Num_cellType"]<-cellType$Freq
Cancer$Tissue_rate<-Tissue[which(Tissue$Var1=="Cancer"),"Freq"]/sum(Tissue$Freq)
Normal<-cellType_Tissue[which(cellType_Tissue$Tissue=="Normal"),]
Normal[match(Normal$cellType,cellType$Var1),"Num_cellType"]<-cellType$Freq
Normal$Tissue_rate<-Tissue[which(Tissue$Var1=="Normal"),"Freq"]/sum(Tissue$Freq)
OLK<-cellType_Tissue[which(cellType_Tissue$Tissue=="OLK"),]
OLK[match(OLK$cellType,cellType$Var1),"Num_cellType"]<-cellType$Freq
OLK$Tissue_rate<-Tissue[which(Tissue$Var1=="OLK"),"Freq"]/sum(Tissue$Freq)
cellType_Tissue<-rbind(Cancer,Normal,OLK)
cellType_Tissue$Num_cellType_split_Tissue_rate<-cellType_Tissue$Num_cellType_split_Tissue/cellType_Tissue$Num_cellType
cellType_Tissue$rate<-cellType_Tissue$Num_cellType_split_Tissue_rate/cellType_Tissue$Tissue_rate - 1

pdf(file="bar_rate_cellType_group_Tissue.pdf", width = 8, height = 6)
ggplot(cellType_Tissue,aes(x=cellType,y=rate,fill=Tissue))+
  geom_bar(stat="identity",position='dodge',width=0.8)+
  labs(x='cellType',y='cellRate')+
  scale_fill_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A"))+
  theme(text = element_text(size=16),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.ticks.y=element_blank())+
  coord_flip()
dev.off()

#Site
orig.ident<-OSCC@meta.data[,c("orig.ident","Site")]
colnames(orig.ident)<-c("cellType","Site")
orig.ident$cellType<-"All T cells"
cellType<-OSCC@meta.data[,c("cellType","Site")]
meta<-rbind(orig.ident,cellType)
pdf(file="bar_percent_cellType_group_Site.pdf", width = 10, height = 6)
ggplot(meta,aes(x=cellType,fill=Site))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='cellType',y='Cell proportion')+
  scale_fill_manual(values=c("#B3456D" ,"#2F869E" ,"#B39C3F"))+
  theme(text = element_text(size=16),axis.ticks.y=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  coord_flip()
dev.off()

cellType<-data.frame(table(OSCC@meta.data$cellType))
Site<-data.frame(table(OSCC@meta.data$Site))
cellType_Site<-data.frame(table(OSCC@meta.data$cellType,OSCC@meta.data$Site))
colnames(cellType_Site)<-c("cellType","Site","Num_cellType_split_Site")
Buccal<-cellType_Site[which(cellType_Site$Site=="Buccal"),]
Buccal[match(Buccal$cellType,cellType$Var1),"Num_cellType"]<-cellType$Freq
Buccal$Site_rate<-Site[which(Site$Var1=="Buccal"),"Freq"]/sum(Site$Freq)
Gingival<-cellType_Site[which(cellType_Site$Site=="Gingival"),]
Gingival[match(Gingival$cellType,cellType$Var1),"Num_cellType"]<-cellType$Freq
Gingival$Site_rate<-Site[which(Site$Var1=="Gingival"),"Freq"]/sum(Site$Freq)
Tongue<-cellType_Site[which(cellType_Site$Site=="Tongue"),]
Tongue[match(Tongue$cellType,cellType$Var1),"Num_cellType"]<-cellType$Freq
Tongue$Site_rate<-Site[which(Site$Var1=="Tongue"),"Freq"]/sum(Site$Freq)
cellType_Site<-rbind(Buccal,Gingival,Tongue)
cellType_Site$Num_cellType_split_Site_rate<-cellType_Site$Num_cellType_split_Site/cellType_Site$Num_cellType
cellType_Site$rate<-cellType_Site$Num_cellType_split_Site_rate/cellType_Site$Site_rate - 1

pdf(file="bar_rate_cellType_group_Site.pdf", width = 8, height = 6)
ggplot(cellType_Site,aes(x=cellType,y=rate,fill=Site))+
  geom_bar(stat="identity",position='dodge',width=0.8)+
  labs(x='cellType',y='cellRate')+
  scale_fill_manual(values=c("#B3456D" ,"#2F869E" ,"#B39C3F"))+
  theme(text = element_text(size=16),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.ticks.y=element_blank())+
  coord_flip()
dev.off()

pdf(file = "Feature_nFeature_nCount_percent.mt.pdf", width = 14, height = 4)
p1 <- FeaturePlot(OSCC,reduction="umap",cols = c("#377EB8","#4DAF4A"),features = "nFeature_RNA")
p2 <- FeaturePlot(OSCC,reduction="umap",cols = c("#984EA3","#FF7F00"),features = "nCount_RNA")
p3 <- FeaturePlot(OSCC,reduction="umap",cols = c("#FFA200","#005DFF"),features = "percent.mt")
plot_grid(p1, p2, p3,ncol = 3)
dev.off()

pdf(file = "VlnPlot_nFeature_nCount_percent.mt.pdf", width = 18, height = 4)
p1 <- VlnPlot(OSCC,cols = cell_type_cols,features = "nFeature_RNA",pt.size = 0)
p2 <- VlnPlot(OSCC,cols = cell_type_cols,features = "nCount_RNA",pt.size = 0)
p3 <- VlnPlot(OSCC,cols = cell_type_cols,features = "percent.mt",pt.size = 0)
plot_grid(p1, p2, p3,ncol = 3)
dev.off()

genes<-c("PTPRC","CD3E","CD4","CD8A","KLRD1","FCGR3A","FOXP3",
         "CD79A","MS4A1","MZB1","LYZ","CD14", "MS4A7", "CD68","LAMP3", "CD80",
         "ACTA2","COL1A2","DCN","LILRA4","IRF7","CMA1",
         "MS4A2","ACTA1", "MYL2","PECAM1","VWF","EPCAM","JUP",
         "KRT5","KRT14","G0S2","CSF3R","FCGR3B")
pdf(file = "Feature_specific_markers.pdf", width = 18, height = 20)
FeaturePlot(OSCC, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=5,
            features =genes)
dev.off()

p1 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              #group.by = "RNA_snn_res.0.6",
              features = genes[1]) + 
  theme_bw() + 
  coord_flip() +
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank() , 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  NoLegend() + panel_border(color = "black")

p2 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[2]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.y = element_blank() , 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


p3 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[3]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.y = element_blank() , 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p4 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[4]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p5 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[5]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p6 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[6]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


p7 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[7]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p8 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[8]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p9 <- VlnPlot(OSCC, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[9]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p10 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[10]) +
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p11 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[11]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),  
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  NoLegend() + panel_border(color = "black")


p12 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[12]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p13 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[13]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p14 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[14]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p15 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[15]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p16 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[16]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p17 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[17]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p18 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[18]) + 
  theme_bw() + 
  coord_flip() +
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank() , 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  NoLegend() + panel_border(color = "black")

p19 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[19]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p20 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[20]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p21 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[21]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p22 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[22]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p23 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[23]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p24 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[24]) +
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p25 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[25]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),  
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  NoLegend() + panel_border(color = "black")


p26 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[26]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p27 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[27]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p28 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[28]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p29 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[29]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p30 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[30]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p31 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[31]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p32 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[32]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p33 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[33]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p34 <- VlnPlot(OSCC, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[34]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

pdf(file = "plot_grid_cluster_marker.pdf", width = 20, height =10)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,
          p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,ncol = 17)
dev.off()

#T:0,1,2,3,6,10,13,15,17,19; B: 11; Plasma cells: 21;  MoMacDC: 8,9,14,24?; pDC: 23; Mesenchymal cells: 7,12,18,24?,25?; 
#Myofibroblasts: 18; Mast cells: 22;  Myocytes: 20; Endothelial cells: 5,16,25?; Epithelial cells: 18; Neutrophils: 4

#res=1.5
#T cells: 0,1,2,3,4,6,9,11,12,17,22,27,28,
#Myeloid cells: 7,15,16, 21,29,30,34,35,36
#Neutrophils: 8,20,
#B cells: 10,37
#Plasma cells: 25
#Mast cells: 26
#Endothelial cells: 5,13,39
#Stromal cells: 14,18,19,23,31,32,38,40
#Myocytes: 24
#Epithelial cells: 33

cellType<-c("T cells","T cells","T cells","T cells","T cells","Endothelial cells","T cells","Myeloid cells",
            "Neutrophils","T cells","B cells","T cells","T cells","Endothelial cells","Stromal cells","Myeloid cells",
            "Myeloid cells","T cells","Stromal cells","Stromal cells","Neutrophils","Myeloid cells","T cells",
            "Stromal cells","Myocytes","Plasma cells","Mast cells","T cells","T cells","Myeloid cells",
            "Myeloid cells","Stromal cells","Stromal cells","Epithelial cells","Myeloid cells","Myeloid cells",
            "Myeloid cells","B cells","Stromal cells","Endothelial cells","Stromal cells")
OSCC$cellType<-factor(OSCC@meta.data$seurat_clusters,levels =0:40,labels =cellType)
a<-table(OSCC@meta.data$cellType)
a<-a[order(a,decreasing = T)]
OSCC@meta.data$cellType<-factor(OSCC@meta.data$cellType,levels = names(a))

OSCC.cellType.cols<-c("#9BA7E2","#F1A5C4","#A9E49D","#FBF785","#DFD2EB", "#D5EB3F", "#E8826D", "#D0B178", "#F781BF","#B6E9D8")
pdf(file="Umap_cellType.pdf", width = 8, height = 6)
DimPlot(OSCC,
        reduction="umap",
        label = F,
        label.size= 4,
        group.by = "cellType",
        cols = cell_type_cols)
DimPlot(OSCC,
        reduction="umap",
        label = F,
        label.size= 4,
        group.by = "cellType",
        cols = OSCC.cellType.cols)
dev.off()

#c24:MoMacDC/Mesenchymal cells
c24<-subset(OSCC,idents = 24)
c24 <- RunUMAP(c24, reduction = "mnn", dims = 1:30, umap.method='umap-learn')
c24 <- FindNeighbors(c24, reduction = "mnn", dims = 1:30)
c24 <- FindClusters(c24,resolution = 0.2)
pdf(file="c24_Umap.pdf", width = 9, height = 6)  
DimPlot(c24, reduction = "umap", cols = cell_type_cols, label = T)
dev.off()
pdf(file = "c24_Feature.pdf", width = 10, height = 8)
FeaturePlot(c24, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=3,
            features =c("LYZ","CD14", "MS4A7", "CD68","COL1A2","DCN","PECAM1","VWF"))
dev.off()
pdf(file = "c24_VlnPlot.pdf", width = 10, height = 8)
VlnPlot(c24, 
        cols = cell_type_cols,
        pt.size = 0,
        features = c("LYZ","CD14", "MS4A7", "CD68","COL1A2","DCN","PECAM1","VWF"))
dev.off()
#Cluster24: Mesenchymal cells/MoMacDC / Endothelial cells ; C0: double: Mesenchymal cells/MoMacDC; C1: double: MoMacDC/ Endothelial cells 
current.cluster.ids<-names(table(c24@meta.data$seurat_clusters))
new.cluster.ids<-c("Double","Double")
c24@meta.data$cellType <- plyr::mapvalues(x = c24@meta.data$seurat_clusters,
                                          from = current.cluster.ids,
                                          to = new.cluster.ids)
table(c24@meta.data$cellType)
c24.meta<-c24@meta.data[,c("orig.ident","cellType")]
c24.meta$cellType<-as.character(c24.meta$cellType)
OSCC.meta<-OSCC@meta.data
OSCC.meta$cellType<-as.character(OSCC.meta$cellType)
OSCC.meta[match(row.names(c24.meta),row.names(OSCC.meta)),"cellType"]<-c24.meta$cellType
OSCC<-AddMetaData(OSCC,OSCC.meta)

#c25:Mesenchymal cells/Endothelial cells/Neutrophils
c25<-subset(OSCC,idents = 25)
c25 <- RunUMAP(c25, reduction = "mnn", dims = 1:30, umap.method='umap-learn')
c25 <- FindNeighbors(c25, reduction = "mnn", dims = 1:30)
c25 <- FindClusters(c25,resolution = 0.2)
pdf(file="c25_Umap.pdf", width = 9, height = 6)  
DimPlot(c25, reduction = "umap", cols = cell_type_cols, label = T)
dev.off()
pdf(file = "c25_Feature.pdf", width = 10, height = 5)
FeaturePlot(c25, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=3,
            features =c("PECAM1","VWF","COL1A2","DCN","G0S2","FCGR3B"))
dev.off()
pdf(file = "c25_VlnPlot.pdf", width = 10, height = 5)
VlnPlot(c25, 
        cols = cell_type_cols,
        pt.size = 0,
        features = c("PECAM1","VWF","COL1A2","DCN","G0S2","FCGR3B"))
dev.off()
#Cluster25:C0/2: Mesenchymal cells / Endothelial cells ;C1: Mesenchymal cells /  Neutrophils
current.cluster.ids<-names(table(c25@meta.data$seurat_clusters))
new.cluster.ids<-c("Double","Double","Double")
c25@meta.data$cellType <- plyr::mapvalues(x = c25@meta.data$seurat_clusters,
                                          from = current.cluster.ids,
                                          to = new.cluster.ids)
table(c25@meta.data$cellType)
c25.meta<-c25@meta.data[,c("orig.ident","cellType")]
c25.meta$cellType<-as.character(c25.meta$cellType)
OSCC.meta<-OSCC@meta.data
OSCC.meta$cellType<-as.character(OSCC.meta$cellType)
OSCC.meta[match(row.names(c25.meta),row.names(OSCC.meta)),"cellType"]<-c25.meta$cellType
OSCC<-AddMetaData(OSCC,OSCC.meta)

Idents(OSCC)<-"cellType"
OSCC<-subset(OSCC,idents = c("T cells", "Neutrophils", "Endothelial cells", "Mesenchymal cells", "MoMacDC", "B cells", 
                                        "Epithelial cells" , "Myocytes", "Plasma cells",  "Mast cells", "pDC"))
save(OSCC,file = "OSCC_del.double.Rda")
load("OSCC_del.double.Rda")

T_cells<-subset(OSCC,idents = "T cells")
save(T_cells,file = "T/T_cells.Rda")

pdf(file="Umap_Celltype_final.pdf", width = 14, height = 6)
DimPlot(OSCC,
        reduction="umap",
        #label = T,
        #label.size= 4,
        #label.box=T,
        #repel =T,
        cols = c(cell_type_cols,cell_type_cols),
        group.by = "cellType")
dev.off()

a<-table(OSCC@meta.data$cellType,OSCC@meta.data$Patient.numberID)
write.csv(a,file = "cellType.Patient.numberID.csv")
b<-table(OSCC@meta.data$Patient.numberID)
write.csv(b,file = "number.Patient.numberID.csv")
c<-OSCC@meta.data[,c("barcode","Patient.numberID")]
write.csv(c,file = "barcode_Patient.numberID.csv")

VlnPlot(OSCC, 
        pt.size = 0,
        group.by = "cellType",
        ncol = 1,
        features = c("CDH5"))
FeaturePlot(OSCC, 
            cols = c("grey80","red"),
            max.cutoff = 3, 
            pt.size = 0.1,
            ncol=1,
            features =c("CD4"))


OSCC@meta.data$Tissue<-factor(OSCC@meta.data$Tissue,levels = c("Normal","OLK","Cancer"))

CD8T@meta.data$Celltype<-as.character(CD8T@meta.data$Celltype)
OSCC@meta.data[which(OSCC@meta.data$cellType=="CD8T cells"),"Celltype"]<-CD8T@meta.data$Celltype
CD4T@meta.data$Celltype<-as.character(CD4T@meta.data$Celltype)
OSCC@meta.data[which(OSCC@meta.data$cellType=="CD4T cells"),"Celltype"]<-CD4T@meta.data$Celltype

OSCC.meta<-OSCC@meta.data
CD4.8T.cell<-OSCC.meta[which(OSCC.meta$cellType=="CD8T cells"|
                               OSCC.meta$cellType=="CD4T cells"|
                               OSCC.meta$cellType=="CD4/8T cells"),]
CD4T.cell<-CD4T@meta.data
CD8T.cell<-CD8T@meta.data

other.cell<-OSCC.meta[which(OSCC.meta$cellType!="CD8T cells"&
                              OSCC.meta$cellType!="CD4T cells"&
                              OSCC.meta$cellType!="CD4/8T cells"),]
cell<-c(row.names(CD4T.cell),row.names(CD8T.cell),row.names(other.cell))

Idents(OSCC)<-row.names(OSCC@meta.data)
OSCC<-subset(OSCC,idents = cell)
OSCC@meta.data[row.names(CD4T.cell),"cellType"]<-"CD4T cells"
OSCC@meta.data[row.names(CD8T.cell),"cellType"]<-"CD8T cells"

CD4T.cell$Celltype<-as.character(CD4T.cell$Celltype)
CD8T.cell$Celltype<-as.character(CD8T.cell$Celltype)

OSCC@meta.data[row.names(CD4T.cell),"Celltype"]<-CD4T.cell$Celltype
OSCC@meta.data[row.names(CD8T.cell),"Celltype"]<-CD8T.cell$Celltype

T_cells.meta<-read.csv("E:/OSCC/integrate/3_5/FastMNN/T/T.del.double.dead/T_cells.del.meta.csv",row.names = 1)
current.cluster.ids<-names(table(T_cells.meta$cellType))
new.cluster.ids<-c("CD8T_CRTAM","CD8T_IFNG","CD8T_GZMK","CD4/8T_RPL34",         
                   "gdT1", "gdT2","NKT","Tfh",             
                   "CD4T_exhausted","CD8T_ISG15","CD8T_mitotic","CD4T_mitotic",         
                   "CD4T_memory_IL7R","gdT_mitotic","CD4Treg_FOXP3_lo","CD8T_CCL5",    
                   "CD4/8T_HSPA1A","CD8T_exhausted","CD4T_naive_CCR7", "Th17" ,           
                   "CD4Treg_FOXP3_hi")
T_cells.meta$cellType <- plyr::mapvalues(x = T_cells.meta$cellType,
                                                 from = current.cluster.ids,
                                                 to = new.cluster.ids)

other.meta<-read.csv("E:/OSCC/integrate/3_5/FastMNN/OSCCsgl_mnnOSCC_singlecell_metadata.csv",row.names = 1)
OSCC.meta<-other.meta[,c("samples","Celltype")]
OSCC.meta[row.names(OSCC.meta) %in% row.names(T_cells.meta),"Celltype"]<-T_cells.meta$cellType
OSCC<-AddMetaData(OSCC,OSCC.meta)

table(other.meta$Sample,other.meta$samples)
current.cluster.ids<-c("HNH.T.Ca","HNH.T.N","CZX.T.Ca","CZX.T.N","CZC.T.Ca","CZC.T.N","ChSF.T.Ca","ChSF.T.N","XZB.OLK","LWH.T.Ca","LWH.T.N",
                       "LYX.T.Ca","WRH.G.Ca","ZhChH.T.Ca","HKJ.T.Ca","HKJ.T.N","ChZhX.T.N","PHD.B.Ca","WBSh.G.Ca","PShY.G.Ca", "PShY.G.N",
                       "LiYX.OLK","LJH.B.Ca","ZhRF.OLK")
new.cluster.ids<-c("Pt01_Ca","Pt01_N","Pt02_Ca","Pt02_N","Pt03_Ca","Pt03_N","Pt04_Ca","Pt04_N","Pt05_OLK","Pt06_Ca","Pt06_N",
                   "Pt07_Ca","Pt08_Ca","Pt09_Ca","Pt10_Ca","Pt10_N","Pt11_N","Pt12_Ca","Pt13_Ca","Pt14_Ca","Pt14_N",
                   "Pt15_OLK","Pt16_Ca","Pt17_OLK")
OSCC@meta.data$Patient.numberID <- plyr::mapvalues(x = OSCC@meta.data$Sample,
                                         from = current.cluster.ids,
                                         to = new.cluster.ids)
pdf(file="Umap_group_Patient.numberID.pdf", width = 9, height = 6)
DimPlot(OSCC, reduction = "umap", group.by = "Patient.numberID",cols = cell_type_cols)
dev.off()

OSCC@meta.data$cellType<-as.character(OSCC@meta.data$cellType)
OSCC@meta.data[which(OSCC@meta.data$cellType=="DC"),"cellType"]<-"DCs"
OSCC@meta.data[which(OSCC@meta.data$cellType=="pDC"),"cellType"]<-"pDCs"
OSCC@meta.data[which(OSCC@meta.data$cellType=="Melanocyte"),"cellType"]<-"Melanocytes"
a<-table(OSCC@meta.data$cellType)
a<-a[order(a,decreasing = T)]
OSCC@meta.data$cellType<-factor(OSCC@meta.data$cellType,levels = names(a))


orig.ident<-OSCC@meta.data[,c("orig.ident","Patient.numberID")]
colnames(orig.ident)<-c("cellType","Patient.numberID")
orig.ident$cellType<-"All T cells"
cellType<-OSCC@meta.data[,c("cellType","Patient.numberID")]
meta<-rbind(orig.ident,cellType)
a<-table(meta$cellType)
a<-a[order(a,decreasing = T)]
meta$cellType<-factor(meta$cellType,levels = names(a))
pdf(file="bar_percent_cellType_group_Patient.numberID.pdf", width = 8, height = 4)
ggplot(meta,aes(x=cellType,fill=Patient.numberID))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='cellType',y='Cell proportion')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16),axis.ticks.y=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  coord_flip()
dev.off()

pdf(file="bar_percent_Patient.numberID_group_cellType.pdf", width = 8, height = 4)
meta<-meta[which(meta$cellType!="All T cells"),]
ggplot(meta,aes(x=Patient.numberID,fill=cellType))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='Cell proportion (%)')+
  scale_fill_manual(values=cell_type_cols)+
  theme_bw()+
  theme(text = element_text(size=12),axis.text.x=element_text(angle =45,hjust =1,vjust =1,colour = "black"),
        axis.text.y=element_text(colour = "black"),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

current.cluster.ids<-names(table(OSCC@meta.data$Celltype))
Idents(OSCC)<-"Celltype"
OSCC<-subset(OSCC,idents = current.cluster.ids)
MET<-subset(OSCC,idents = c("Mac-CCL18","Mac-CCL3","Mac-LYVE1","Mac-SPP1","Mac-TREM2","Mono-FCN1","Mono-S100A9","cDC1", "cDC2-CD1C","cDC2-FCER1A","cDC3",
                            "ADSC-MFAP5", "ADSC-SFRP1", "ADSC-VEGFD","Fibro-CCL19","Fibro-IGF1","Fibro-IGFBP2","Fibro-THSD4","MF-MMP11","MF-TDO2","Epithelial cell"))
save(MET,file = "E:/OSCC/integrate/3_5/FastMNN/MET/infercnv_Myeloid/MET.Rda")

Epi.MF<-subset(OSCC,idents = c("ADSC-MFAP5", "ADSC-SFRP1", "ADSC-VEGFD","Fibro-CCL19","Fibro-IGF1","Fibro-IGFBP2","Fibro-THSD4","MF-MMP11","MF-TDO2","Epithelial cell"))
save(Epi.MF,file = "E:/OSCC/integrate/3_5/FastMNN/Epi.MF/Epi.MF.Rda")


current.cluster.ids<-names(table(OSCC@meta.data$Celltype))
new.cluster.ids<-c("ADSC","ADSC","ADSC","Arteriole","B","B",
                   "B","B","B","Basal cell","CD4/8T","CD4/8T",
                   "CD4T","CD4T","CD4T","CD4T","CD4Treg","CD4Treg",
                   "CD8T","CD8T","CD8T","CD8T","CD8T","CD8T",
                   "CD8T","cDC1","cDC2","cDC2","cDC3","Epithelial cell",
                   "Fibro","Fibro","Fibro","Fibro","gdT","gdT",
                   "gdT","LEC","LEC","LEC","Mac","Mac",
                   "Mac","Mac","Mac","Mast cells","Melanocyte","MyoFib",
                   "MyoFib","Mono","Mono","Myoblast","Myocytes","Neuron",
                   "Neutro","Neutro","Neutro","Neutro","Neutro","NKT",
                   "PCV","PCV","PCV","pDC","Pericyte","Plasma",
                   "Plasma","Plasma","Plasma","Secretory cell","Smooth muscle","Tfh",
                   "Th17","Tip")
OSCC@meta.data$Celltype.short <- plyr::mapvalues(x = OSCC@meta.data$Celltype,
                                                       from = current.cluster.ids,
                                                       to = new.cluster.ids)

pdf(file="Umap_Celltype.short_del.dead.double.pdf", width = 10, height = 6)
DimPlot(OSCC,
        reduction="umap",
        label = T,
        label.size= 4,
        label.box=T,
        repel =T,
        cols = c(cell_type_cols,cell_type_cols),
        group.by = "Celltype.short")
DimPlot(OSCC,
        reduction="umap",
        label = T,
        label.size= 4,
        repel =T,
        cols = c(cell_type_cols,cell_type_cols),
        group.by = "Celltype.short")
dev.off()

current.cluster.ids<-names(table(OSCC@meta.data$Celltype))
new.cluster.ids<-c("Mesenchymal cells","Mesenchymal cells","Mesenchymal cells","BECs","B cells","B cells",
                   "B cells","B cells","B cells","Mesenchymal cells","T cells","T cells",
                   "T cells","T cells","T cells","T cells","T cells","T cells",
                   "CD8T cells","CD8T cells","CD8T cells","CD8T cells","CD8T cells","CD8T cells",
                   "CD8T cells","DC","DC","DC","DC","Epithelial cells",
                   "Mesenchymal cells","Mesenchymal cells","Mesenchymal cells","Mesenchymal cells","gdT cells","gdT cells",
                   "gdT cells","LECs","LECs","LECs","Macrophages","Macrophages",
                   "Macrophages","Macrophages","Macrophages","Mast cells","Melanocyte","Mesenchymal cells",
                   "Mesenchymal cells","Monocytes","Monocytes","Myoblasts","Myocytes",
                   "Neutrophils","Neutrophils","Neutrophils","Neutrophils","Neutrophils","NKT cells",
                   "BECs","BECs","BECs","pDC","BECs","Plasma cells",
                   "Plasma cells","Plasma cells","Plasma cells","Epithelial cells","Smooth muscle cells","CD4T cells",
                   "CD4T cells","BECs")
OSCC@meta.data$cellType <- plyr::mapvalues(x = OSCC@meta.data$Celltype,
                                                 from = current.cluster.ids,
                                                 to = new.cluster.ids)
OSCC@meta.data$cellType<-as.character(OSCC@meta.data$cellType)
OSCC@meta.data[which(OSCC@meta.data$Celltype=="Basal cell"),"cellType"]<-"Epithelial cells"
a<-table(OSCC@meta.data$cellType)
a<-a[order(a,decreasing = T)]
OSCC@meta.data$cellType<-factor(OSCC@meta.data$cellType,levels = names(a))
save(OSCC,file = "OSCC_del.dead.double.Rda")

Idents(OSCC)<-"cellType"
OSCC<-subset(OSCC,idents = c("CD8T cells","CD4T cells","Mesenchymal cells","Macrophages","BECs"               
                                         ,"Neutrophils","gdT cells","DC","B cells","Monocytes"          
                                         ,"LECs","NKT cells","Plasma cells","Mast cells","Myocytes"           
                                         ,"pDC","Smooth muscle cells","Myoblasts","Epithelial cells"   
                                         ,"Melanocyte"))
pdf(file="Umap_cellType_final.pdf", width = 8, height = 6)
DimPlot(OSCC,
        reduction="umap",
        label = T,
        label.size= 4,
        #label.box=T,
        repel =T,
        cols = cell_type_cols,
        group.by = "cellType",
        raster=FALSE)
DimPlot(OSCC,
        reduction="umap",
        label = F,
        cols = cell_type_cols,
        group.by = "cellType",
        raster=FALSE)
dev.off()

OSCC@meta.data[which(OSCC@meta.data$cellType=="T cells" |
                             OSCC@meta.data$cellType=="Myeloid cells" |
                             OSCC@meta.data$cellType=="Neutrophils" |
                             OSCC@meta.data$cellType=="B cells" |
                             OSCC@meta.data$cellType=="Plasma cells" |
                             OSCC@meta.data$cellType=="Mast cells"),"CD45"]<-"Immune cells"
OSCC@meta.data[which(OSCC@meta.data$cellType=="Stromal cells" |
                             OSCC@meta.data$cellType=="Endothelial cells" |
                             OSCC@meta.data$cellType=="Myocytes" |
                             OSCC@meta.data$cellType=="Epithelial cells"),"CD45"]<-"Non-immune cells"

pdf(file="Umap_CD45.pdf", width = 7, height = 6)
DimPlot(OSCC,
        reduction="umap",
        label = T,
        label.size= 4,
        repel =T,
        cols = c("#40AD74", "#3F7AAB"),
        group.by = "CD45")+NoLegend()
dev.off()

meta<-OSCC@meta.data[,c("CD45","Tissue")]
pdf(file="bar_percent_CD45_Tissue.pdf", width = 6, height = 5)
ggplot(meta,aes(x=Tissue,fill=CD45))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='cellType',y='Cell proportion')+
  scale_fill_manual(values=c("#40AD74", "#3F7AAB"))+
  theme(text = element_text(size=16),axis.ticks.y=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

Mac<-subset(OSCC,idents = "Macrophages")
save(Mac,file="Macrophages/Macrophages.Rda")

DC<-subset(OSCC,idents = "DC")
save(DC,file="DC/DC.Rda")

Mesenchymal<-subset(OSCC,idents = "Mesenchymal cells")
save(Mesenchymal,file="Mesenchymal/Mesenchymal.Rda")

Idents(OSCC)<-"Celltype.short"
MyoFib<-subset(OSCC,idents = "MyoFib")
save(MyoFib,file="MyoFib/MyoFib.Rda")

pdf(file = "Feature_IFN_markers.pdf", width = 12, height = 4)
FeaturePlot(OSCC, 
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=3,
            features =c("IFNB1","IFNG","IFNL1"))
dev.off()

pdf(file = "Feature_EPCAM.pdf", width = 8, height = 6)
FeaturePlot(OSCC, 
            cols = c("grey80","red"),
            max.cutoff = 3, 
            pt.size = 0.1,
            ncol=1,
            features =c("EPCAM"))
dev.off()
pdf(file = "VlnPlot_EPCAM.pdf", width = 8, height = 4)
VlnPlot(OSCC, 
        pt.size = 0.1,
        group.by = "cellType",
        cols = cell_type_cols,
        features = c("EPCAM"))
dev.off()

p1<-VlnPlot(OSCC, 
        pt.size = 0,
        group.by = "Celltype",
        features = c("IDO1"))+
  theme(axis.text.x = element_blank(),axis.title = element_blank(),plot.title = element_text(hjust = 0.5))+NoLegend() 
p2<-VlnPlot(OSCC, 
            pt.size = 0,
            group.by = "Celltype",
            features = c("IDO2"))+
  theme(axis.text.x = element_blank(), axis.title = element_blank(),plot.title = element_text(hjust = 0.5))+NoLegend()
p3<-VlnPlot(OSCC, 
            pt.size = 0,
            group.by = "Celltype",
            features = c("TDO2"))+
  theme(axis.text.x = element_text(size =11,angle =40),axis.title = element_blank(),plot.title = element_text(hjust = 0.5))+ NoLegend()
pdf(file = "IDO_TDO.pdf", width = 20, height =8)
plot_grid(p1,p2,p3,ncol = 1)
dev.off()

load("OSCC_cellType.markers.Rda")
write.csv(OSCC.markers,file = "OSCC_cellType.markers.csv")

genes<-c("CD3E","CD3D","CD8A","CD8B","CD4","ICOS","DCN","COL1A2","CD68","C1QA", 
         "PECAM1","VWF","G0S2","CSF3R","TRDC","TRGC1","CST3", "CD83","MS4A1","CD79A","FCN1","CD14", 
         "PROX1","LYVE1","KLRD1","FCGR3A","MZB1","IGKC","MS4A2","TPSB2","ACTA1", "MYL2","LILRA4","IRF4",
         "KRT14","KRT5","ACTA2","TAGLN","MYF5","PAX7","PMEL","DCT")

genes<-c("PTPRC","CD3E","CD3D","CD3G","LYZ","CD14","C1QB","C1QA",
         "DCN","COL1A1","COL3A1","COL1A2","TM4SF1","PECAM1","VWF","SELE",
         "CXCL8","G0S2","CSF3R","FCGR3B","MS4A1","CD79A","CD79B", "CD19",
         "ACTA1", "MYL1","MYH2","MYL2","MZB1","DERL3","IGKC","IGLC2",
         "TPSB2","TPSAB1","CPA3","MS4A2","KRT14","KRT5","KRT17","JUP")
pdf(file="DotPlot_cellType.pdf", width = 18, height = 6)
DotPlot(OSCC,
        features=genes,
        group.by ="cellType",
        cols =c("#4B60E1","#FF555E"))+
  RotatedAxis() + xlab(NULL) + ylab(NULL)+
  theme(text = element_text(size = 11))
dev.off()

load("OSCC_cellType.markers.Rda")
head(OSCC.markers)
OSCC.markers<-OSCC.markers[which(OSCC.markers$cluster!="Neuron"),]
table(OSCC.markers$cluster)
OSCC.markers$cluster<-droplevels(OSCC.markers$cluster)
top <- OSCC.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
OSCC@meta.data$cellType<-droplevels(OSCC@meta.data$cellType)
pdf(file="DoHeatmap.cellType.pdf", width = 10, height = 8)
DoHeatmap(OSCC, features = a,group.by="cellType", label = F,group.colors =cell_type_cols,slot = "data") +
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                       mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")
dev.off()

a<-top$gene
a<-a[!duplicated(a)]
DotPlot(OSCC,
        features=a,
        group.by ="cellType",
        cols =c("#4B60E1","#FF555E"))+
  RotatedAxis() + xlab(NULL) + ylab(NULL)


