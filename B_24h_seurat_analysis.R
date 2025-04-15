####################################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
# Run on R.4.3.1
#####################################
#WORKING_DIR = "<enter your working dir here>" 
#setwd(WORKING_DIR)
library(cowplot)
library(dplyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)
library(pheatmap)
library(data.table)

##### color definition ######
pluri_color = "#FFCC99"
te_color = "#339999"
pre_color = "#CC6666"
twoc_color = "#9933CC"

##### Load data after CellRanger #####
B_24h.input <- Read10X(data.dir = "<path_to_data>/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
B_24h.data <- CreateSeuratObject(counts = B_24h.input, project = "SEM2", min.cells = 3, min.features = 200)

##### Filter ######
B_24h.data$sample = "24h"
B_24h.data[["percent.mt"]] <- PercentageFeatureSet(B_24h.data, pattern = "^mt-")
B_24h.data_sub <- subset(B_24h.data, subset = nFeature_RNA > 3000 & nFeature_RNA <8000 & percent.mt < 10)
B_24h.data = B_24h.data_sub

######### Normalization, UMAP and Clustering ##############
dims=1:10
B_24h.data <- NormalizeData(B_24h.data)
B_24h.data <- FindVariableFeatures(B_24h.data, selection.method = "vst", nfeatures = 2000)
B_24h.data <- ScaleData(B_24h.data, verbose = FALSE)
B_24h.data <- RunPCA(B_24h.data, npcs = 30, verbose = FALSE)
B_24h.data <- RunUMAP(B_24h.data, reduction = "pca", dims = dims,min.dist = 0.1)
B_24h.data <- FindNeighbors(B_24h.data, reduction = "pca", dims = dims)
B_24h.data <- FindClusters(B_24h.data, resolution = 0.1)
ElbowPlot(B_24h.data)

##### start here ############
#saveRDS(B_24h.data, file = "B_24h.rds")
#B_24h.data<-readRDS("B_24h.rds")

###### Visualization
DimPlot(B_24h.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)

##### count visualization #######
VlnPlot(B_24h.data, features = c("nCount_RNA"), ncol = 1,pt.size = 0,same.y.lims=F,log = T,group.by="seurat_clusters")+theme(legend.position = "none",plot.title = element_blank())+xlab("")
FeaturePlot(B_24h.data,"nCount_RNA")
FeaturePlot(B_24h.data,"nFeature_RNA")

##### Find cluster markers #######
B_24h_markers = FindAllMarkers(B_24h.data,logfc.threshold = 0.5,only.pos =T )
write.csv(B_24h_markers,"B_24h_markers.csv")

######### Marker annotation ###############
EPI = c("Nanog","Sox2","Pou5f1","Klf4","Dppa5a","Zfp42","Tfcp2l1","Fgf4","Utf1","Sall4","Esrrb","Dppa3","Tbx3","Klf5")
TE = c("Cdx2","Id2","Elf5","Lgals1","Krt8","Eomes","Id2","Cldn7","Tfap2c","Gata3","Krt18")
PrE = c("Sox17","Gata6","Pdgfra","Dab2","Ttr","Apoa1","Foxa2","Gata4","Serpinh1","Amn","Srgn")
pluri = c("Sox2","Pou5f1")
Naive = c("Sox2","Pou5f1","Nanog","Klf4","Dppa5a","Zfp42","Tfcp2l1","Fgf4","Utf1","Sall4","Esrrb","Dppa3","Tbx3","Klf5")
Primed_Formative = c("Pou5f1","Sox2","Otx2")
zscans = c("Zscan4-ps2","Zscan4-ps1","Zscan4c","Zscan4d","Zscan4e","Zscan4-ps3","Zscan4b","Zscan4f")

FeaturePlot(B_24h.data,features = Naive,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(B_24h.data,features = TE,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(B_24h.data,features = PrE,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(B_24h.data,features = zscans,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()

groups<- B_24h.data@meta.data$RNA_snn_res.0.1
type<- B_24h.data@meta.data$RNA_snn_res.0.1
type<-ifelse(groups==0,type[groups]<-"EPI",type[groups]<-type)
type<-ifelse(groups==1,type[groups]<-"2CL",type[groups]<-type)
type<-ifelse(groups==2,type[groups]<-"PrE",type[groups]<-type)
B_24h.data$type<-type

DimPlot(B_24h.data, reduction = "umap", group.by="type",label=F,pt.size = 1,cols=c(twoc_color,pluri_color,pre_color)) & NoAxes()

######### Cell cycle analysis ###############
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
B_24h.data <- CellCycleScoring(B_24h.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(object = B_24h.data,reduction = "umap",repel = T,group.by = "Phase",label=T,label.size = 2)

