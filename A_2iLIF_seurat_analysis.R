####################################
# Written by Noa Novershtern 
# 
#### Run on R.4.3.1 ################
#WORKING_DIR = "<enter your working dir here" 
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
library(plotly)

##### color definition ######
pluri_color = "#FFCC99"
te_color = "#339999"
pre_color = "#CC6666"
twoc_color = "#9933CC"

##### Load data after CellRanger #####
A_2i_LIF.input <- Read10X(data.dir = "/home/labs/hanna/noanov/SynEmb2.0/CellRanger/2i_LIF/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
A_2i_LIF.data <- CreateSeuratObject(counts = A_2i_LIF.input, project = "TFSEM", min.cells = 3, min.features = 200)
 
##### Filter ######
A_2i_LIF.data$sample = "A_2iLIF"
A_2i_LIF.data[["percent.mt"]] <- PercentageFeatureSet(A_2i_LIF.data, pattern = "^mt-")
A_2i_LIF.data_sub <- subset(A_2i_LIF.data, subset = nFeature_RNA > 3000 & nFeature_RNA <8000 & percent.mt < 10)
A_2i_LIF.data = A_2i_LIF.data_sub

######### Normalization, UMAP and Clustering ##############
dims=1:10

A_2i_LIF.data <- NormalizeData(A_2i_LIF.data)
A_2i_LIF.data <- FindVariableFeatures(A_2i_LIF.data, selection.method = "vst", nfeatures = 2000)
A_2i_LIF.data <- ScaleData(A_2i_LIF.data, verbose = FALSE)
A_2i_LIF.data <- RunPCA(A_2i_LIF.data, npcs = 30, verbose = FALSE)
A_2i_LIF.data <- RunUMAP(A_2i_LIF.data, reduction = "pca", dims = dims,min.dist = 0.1)
A_2i_LIF.data <- FindNeighbors(A_2i_LIF.data, reduction = "pca", dims = dims)
A_2i_LIF.data <- FindClusters(A_2i_LIF.data, resolution = 0.1)#seq(0,1.2,0.1))

ElbowPlot(A_2i_LIF.data)

# saveRDS(A_2i_LIF.data, file = "A_2i_LIF.rds")
# A_2i_LIF.data<-readRDS("A_2i_LIF.rds")

###### Visualization ########
DimPlot(A_2i_LIF.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)

##### count visualization #######
VlnPlot(A_2i_LIF.data, features = c("nCount_RNA"), ncol = 1,pt.size = 0,same.y.lims=F,log = T,group.by="seurat_clusters")+theme(legend.position = "none",plot.title = element_blank())+xlab("")
FeaturePlot(A_2i_LIF.data,"nCount_RNA")
FeaturePlot(A_2i_LIF.data,"nFeature_RNA")

##### Find cluster markers #######
A_2i_LIF_markers = FindAllMarkers(A_2i_LIF.data,logfc.threshold = 0.5,only.pos =T )
write.csv(A_2i_LIF_markers,"A_2i_LIF_markers.csv")


######### Marker annotation ###############
pluri = c("Sox2","Pou5f1","Nanog","Klf4","Dppa5a","Zfp42","Tfcp2l1","Fgf4","Utf1","Esrrb","Dppa3","Tbx3","Klf5")
twoC = c("Zscan4c","Zscan4d","Zbed3","Zscan4c","Zscan4d","Zscan4e","Obox1","Tcstv1","Tcstv3")
TE = c("Cdx2","Id2","Elf5","Lgals1","Krt8","Eomes","Id2","Cldn7","Tfap2c","Gata3","Krt18")
PrE = c("Sox17","Gata6","Pdgfra","Dab2","Apoa1","Foxa2","Gata4","Serpinh1","Amn","Srgn")
Primed = c("Fgf5","Otx2","Zic2")

FeaturePlot(A_2i_LIF.data,features = pluri,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(A_2i_LIF.data,features = TE,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(A_2i_LIF.data,features = PrE,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()
FeaturePlot(A_2i_LIF.data,features = twoC,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0) & NoAxes()

groups<- A_2i_LIF.data@meta.data$RNA_snn_res.0.1
type<- A_2i_LIF.data@meta.data$RNA_snn_res.0.1
type<-ifelse(groups==0,type[groups]<-"pluri",type[groups]<-type)
type<-ifelse(groups==1,type[groups]<-"pluri",type[groups]<-type)
type<-ifelse(groups==2,type[groups]<-"2CL",type[groups]<-type)
A_2i_LIF.data$type<-type

DimPlot(A_2i_LIF.data, reduction = "umap", group.by="type",label=F,pt.size = 1,cols=c(twoc_color,pluri_color)) & NoAxes()
saveRDS(A_2i_LIF.data, file = "A_2i_LIF.rds")


######### Cell cycle analysis ###############
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# view cell cycle scores and phase assignments
A_2i_LIF.data <- CellCycleScoring(A_2i_LIF.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(object = A_2i_LIF.data,reduction = "umap",repel = T,group.by = "Phase",label=T,label.size = 2)

