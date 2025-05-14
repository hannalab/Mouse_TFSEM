####################################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
# Run on R.4.3.1
#####################################
#WORKING_DIR = "<enter your working dir here>" 
#setwd(WORKING_DIR)
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
D_priming.input <- Read10X(data.dir = "<path_to_data>/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
D_priming.data <- CreateSeuratObject(counts = D_priming.input, project = "SEM2", min.cells = 3, min.features = 200)

##### Filter ######
D_priming.data$sample = "96h"
D_priming.data[["percent.mt"]] <- PercentageFeatureSet(D_priming.data, pattern = "^mt-")
D_priming.data_sub <- subset(D_priming.data, subset = nFeature_RNA > 3000 & nFeature_RNA <8000 & percent.mt < 10)
D_priming.data = D_priming.data_sub

######### Normalization, UMAP and Clustering ##############
dims=1:10
D_priming.data <- NormalizeData(D_priming.data)
D_priming.data <- FindVariableFeatures(D_priming.data, selection.method = "vst", nfeatures = 2000)
D_priming.data <- ScaleData(D_priming.data, verbose = FALSE)
D_priming.data <- RunPCA(D_priming.data, npcs = 30, verbose = FALSE)
D_priming.data <- RunUMAP(D_priming.data, reduction = "pca", dims = dims,min.dist = 0.1)
D_priming.data <- FindNeighbors(D_priming.data, reduction = "pca", dims = dims)
D_priming.data <- FindClusters(D_priming.data, resolution = 0.1)#seq(0,1.2,0.1))

ElbowPlot(D_priming.data)

#saveRDS(D_priming.data, file = "D_96h.rds")
#D_priming.data<-readRDS("D_96h.rds")

###### Visualization
DimPlot(D_priming.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)

##### count visualization #######
VlnPlot(D_priming.data, features = c("nCount_RNA"), ncol = 1,pt.size = 0,same.y.lims=F,log = T,group.by="seurat_clusters")+theme(legend.position = "none",plot.title = element_blank())+xlab("")
FeaturePlot(D_priming.data,"nCount_RNA")
FeaturePlot(D_priming.data,"nFeature_RNA")

##### Find cluster markers #######
D_priming_markers = FindAllMarkers(D_priming.data,logfc.threshold = 0.5,only.pos =T)
write.csv(D_priming_markers,"D_96h_markers.csv")

######### Marker annotation ###############
pluri = c("Sox2","Pou5f1","Nanog","Klf4","Dppa5a","Zfp42","Tfcp2l1","Fgf4","Utf1","Esrrb","Dppa3","Tbx3","Klf5")
twoC = c("Zscan4c","Zscan4d","Zbed3","Zscan4c","Zscan4d","Zscan4e","Obox1","Tcstv1","Tcstv3")
TE = c("Cdx2","Id2","Elf5","Lgals1","Krt8","Eomes","Id2","Cldn7","Tfap2c","Gata3","Krt18")
PrE = c("Sox17","Gata6","Pdgfra","Dab2","Apoa1","Foxa2","Gata4","Serpinh1","Amn","Srgn")
Primed = c("Fgf5","Otx2","Zic2")

FeaturePlot(D_priming.data, features = pluri,min.cutoff = 0) & NoAxes()
FeaturePlot(D_priming.data, features = twoC,min.cutoff = 0) & NoAxes()
FeaturePlot(D_priming.data, features = TE,min.cutoff = 0) & NoAxes()
FeaturePlot(D_priming.data, features = PrE,min.cutoff = 0) & NoAxes()
FeaturePlot(D_priming.data, features = Primed,min.cutoff = 0) & NoAxes()

groups<- D_priming.data@meta.data$RNA_snn_res.0.1
type<- D_priming.data@meta.data$RNA_snn_res.0.1
type<-ifelse(groups==0,type[groups]<-"TE",type[groups]<-type)
type<-ifelse(groups==1,type[groups]<-"EPI",type[groups]<-type)
type<-ifelse(groups==2,type[groups]<-"PrE",type[groups]<-type)
type<-ifelse(groups==3,type[groups]<-"TE",type[groups]<-type)
D_priming.data$type<-type

DimPlot(D_priming.data, reduction = "umap", group.by="type",label=F,pt.size = 1,cols=c(pluri_color,pre_color,te_color)) & NoAxes()

######### Cell cycle analysis ###############
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
D_priming.data <- CellCycleScoring(D_priming.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(object = D_priming.data,reduction = "umap",repel = T,group.by = "Phase",label=T,label.size = 2)

#saveRDS(D_priming.data, file = "D_96h.rds")



