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
library(Seurat) # V5.0.0
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
C_48h.input <- Read10X(data.dir = "<path_to_data>/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
C_48h.data <- CreateSeuratObject(counts = C_48h.input, project = "SEM2", min.cells = 3, min.features = 200)

##### filter ######
C_48h.data$sample = "48h"
C_48h.data[["percent.mt"]] <- PercentageFeatureSet(C_48h.data, pattern = "^mt-")
C_48h.data_sub <- subset(C_48h.data, subset = nFeature_RNA > 3000 & nFeature_RNA <8000 & percent.mt < 10)
C_48h.data = C_48h.data_sub

######### Normalization, UMAP and Clustering ##############
dims=1:10
C_48h.data <- NormalizeData(C_48h.data)
C_48h.data <- FindVariableFeatures(C_48h.data, selection.method = "vst", nfeatures = 2000)
C_48h.data <- ScaleData(C_48h.data, verbose = FALSE)
C_48h.data <- RunPCA(C_48h.data, npcs = 30, verbose = FALSE)
C_48h.data <- RunUMAP(C_48h.data, reduction = "pca", dims = dims,min.dist = 0.1)
C_48h.data <- FindNeighbors(C_48h.data, reduction = "pca", dims = dims)
C_48h.data <- FindClusters(C_48h.data, resolution = 0.2) 

ElbowPlot(C_48h.data)

#saveRDS(C_48h.data, file = "C_48h.rds")
#C_48h.data<-readRDS("C_48h.rds")

###### Visualization
DimPlot(C_48h.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)

##### count visualization #######
VlnPlot(C_48h.data, features = c("nCount_RNA"), ncol = 1,pt.size = 0,same.y.lims=F,log = T,group.by="seurat_clusters")+theme(legend.position = "none",plot.title = element_blank())+xlab("")
FeaturePlot(C_48h.data,"nCount_RNA")
FeaturePlot(C_48h.data,"nFeature_RNA")

##### Find cluster markers #######
C_48h_markers = FindAllMarkers(C_48h.data,logfc.threshold = 0.5,only.pos =T )
write.csv(C_48h_markers,"C_48h_markers.csv")

######### Marker annotation ###############
pluri = c("Sox2","Pou5f1","Nanog","Klf4","Dppa5a","Zfp42","Tfcp2l1","Fgf4","Utf1","Esrrb","Dppa3","Tbx3","Klf5")
twoC = c("Zscan4c","Zscan4d","Zbed3","Zscan4c","Zscan4d","Zscan4e","Obox1","Tcstv1","Tcstv3")
TE = c("Cdx2","Id2","Elf5","Lgals1","Krt8","Eomes","Id2","Cldn7","Tfap2c","Gata3","Krt18")
PrE = c("Sox17","Gata6","Pdgfra","Dab2","Apoa1","Foxa2","Gata4","Serpinh1","Amn","Srgn")
Primed = c("Fgf5","Otx2","Zic2")

FeaturePlot(C_48h.data,features = pluri,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(C_48h.data,features = TE,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(C_48h.data,features = PrE,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(C_48h.data,features = twoC,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(C_48h.data,features = Primed,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()

groups<- C_48h.data@meta.data$RNA_snn_res.0.2
type<- C_48h.data@meta.data$RNA_snn_res.0.2
type<-ifelse(groups==0,type[groups]<-"EPI",type[groups]<-type)
type<-ifelse(groups==1,type[groups]<-"PrE",type[groups]<-type)
type<-ifelse(groups==2,type[groups]<-"2CL",type[groups]<-type)
type<-ifelse(groups==3,type[groups]<-"TE",type[groups]<-type)
C_48h.data$type<-type

DimPlot(C_48h.data, reduction = "umap", group.by="type",label=F,pt.size = 1,cols=c(twoc_color,pluri_color,pre_color,te_color)) & NoAxes()
#saveRDS(C_48h.data, file = "C_48h.rds")

######### Cell cycle analysis ###############
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# view cell cycle scores and phase assignments
C_48h.data <- CellCycleScoring(C_48h.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(object = C_48h.data,reduction = "umap",repel = T,group.by = "Phase",label=T,label.size = 2)


