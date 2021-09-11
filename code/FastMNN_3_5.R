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
