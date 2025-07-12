####################################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
# Run on R.4.3.1
#####################################
#WORKING_DIR = "<enter your working dir here>" 
setwd(WORKING_DIR)
library(cowplot)
library(dplyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)
library(pheatmap)
library(data.table)
library(ggtern)

###### color definition ######
pluri_color = "#FFCC99"
te_color = "#339999"
pre_color = "#CC6666"
early_pre_color = "#FF9999"
twoc_color = "#9933CC"
fibro_color = "#3366FF"
morula_color="blue"

####### Load data after CellRanger #####
ACPC12.input <- Read10X(data.dir = "<path to data>/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
ACPC12.data <- CreateSeuratObject(counts = ACPC12.input, project = "SEM2", min.cells = 3, min.features = 200)

####### Filter ######
ACPC12.data$sample = "ACPC12"
ACPC12.data[["percent.mt"]] <- PercentageFeatureSet(ACPC12.data, pattern = "^mt-")
ACPC12.data_sub <- subset(ACPC12.data, subset = nFeature_RNA > 5000 & nFeature_RNA <10000 & percent.mt < 10)
ACPC12.data = ACPC12.data_sub

######### Normalization, UMAP and Clustering ##############
dims=1:10
ACPC12.data <- NormalizeData(ACPC12.data)
ACPC12.data <- FindVariableFeatures(ACPC12.data, selection.method = "vst", nfeatures = 2000)
ACPC12.data <- ScaleData(ACPC12.data, verbose = FALSE)
ACPC12.data <- RunPCA(ACPC12.data, npcs = 30, verbose = FALSE)
ACPC12.data <- RunUMAP(ACPC12.data, reduction = "pca", dims = dims,min.dist = 0.1)
ACPC12.data <- FindNeighbors(ACPC12.data, reduction = "pca", dims = dims)
ACPC12.data <- FindClusters(ACPC12.data, resolution = 0.22)#seq(0,1.2,0.1))

ElbowPlot(ACPC12.data)

#saveRDS(ACPC12.data, file = "ACPC12.rds")
#ACPC12.data<-readRDS("ACPC12.rds")

####### Visualization ########
DimPlot(ACPC12.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)

####### Count visualization #######
VlnPlot(ACPC12.data, features = c("nFeature_RNA"), ncol = 1,pt.size = 0,same.y.lims=F,log = T,group.by="seurat_clusters")+theme(legend.position = "none",plot.title = element_blank())+xlab("")
FeaturePlot(ACPC12.data,"nCount_RNA")
FeaturePlot(ACPC12.data,"nFeature_RNA")

####### Find cluster markers #######
ACPC12_markers = FindAllMarkers(ACPC12.data,logfc.threshold = 0.5,only.pos =T )
write.csv(ACPC12_markers,"ACPC12_markers.csv")

######### Marker annotation ###############
EPI = c("Nanog","Sox2","Pou5f1","Klf4","Dppa5a","Zfp42","Tfcp2l1","Fgf4","Utf1","Sall4","Esrrb","Dppa3","Tbx3","Klf5")
PrE = c("Sox17","Gata6","Pdgfra","Dab2","Ttr","Apoa1","Foxa2","Gata4","Serpinh1","Amn","Srgn")
TE = c("Cdx2","Id2","Elf5","Lgals1","Krt8","Eomes","Cldn7","Tfap2c","Gata3","Krt18")
MEF = c("Fn1","Vim","Col1a1","Acta2","Thy1","Tagln")

FeaturePlot(ACPC12.data, features = TE,min.cutoff = 0.5,max.cutoff=2,ncol=3,order=T,cols=c("lightgrey",te_color)) & NoAxes() #& NoLegend()
FeaturePlot(ACPC12.data, features = PrE,min.cutoff = 0.5,max.cutoff=2,ncol=3,order=T,cols=c("lightgrey",pre_color)) & NoAxes() & NoLegend()
FeaturePlot(ACPC12.data, features = EPI,min.cutoff = 0.5,ncol=4,order=T,cols=c("lightgrey",pluri_color)) & NoAxes() #& NoLegend()
FeaturePlot(ACPC12.data,features = MEF,pt.size = 0.05,order=T,cols=c("lightgrey",fibro_color), ncol = 3,min.cutoff = 0.5) & NoAxes()

groups<- ACPC12.data$seurat_clusters
type<- ACPC12.data$seurat_clusters
type<-ifelse(groups==0,type[groups]<-"E_MEF",type[groups]<-type)
type<-ifelse(groups==1,type[groups]<-"D_TE",type[groups]<-type)
type<-ifelse(groups==2,type[groups]<-"A_EPI",type[groups]<-type)
type<-ifelse(groups==3,type[groups]<-"B_Early_PrE",type[groups]<-type)
type<-ifelse(groups==4,type[groups]<-"C_PrE",type[groups]<-type)
unique(type)
ACPC12.data$type<-type

######## Cell cycle analysis #########
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# view cell cycle scores and phase assignments
ACPC12.data <- CellCycleScoring(ACPC12.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(object = ACPC12.data,reduction = "umap",repel = T,group.by = "Phase",label=T,label.size = 2)

#########################################################################
######### Filter Seurat object not to include fibroblast cluster ########
#########################################################################
clean_ACPC12 =  subset(ACPC12.data,idents=c("1","2","3","4"))

######### Normalization, UMAP and Clustering ##############
dims=1:10
clean_ACPC12 <- NormalizeData(clean_ACPC12)
clean_ACPC12 <- FindVariableFeatures(clean_ACPC12, selection.method = "vst", nfeatures = 2000)
clean_ACPC12 <- ScaleData(clean_ACPC12, verbose = FALSE)
clean_ACPC12 <- RunPCA(clean_ACPC12, npcs = 30, verbose = FALSE)
clean_ACPC12 <- RunUMAP(clean_ACPC12, reduction = "pca", dims = dims,min.dist = 0.1)
clean_ACPC12 <- FindNeighbors(clean_ACPC12, reduction = "pca", dims = dims)
clean_ACPC12 <- FindClusters(clean_ACPC12, resolution = 0.22)#seq(0,1.2,0.1))

######## Visualization ########
DimPlot(clean_ACPC12, reduction = "umap", group.by="type",label=T,pt.size = 1,cols=c(pluri_color,early_pre_color,pre_color,te_color))
VlnPlot(clean_ACPC12, features = c("nFeature_RNA"), ncol = 1,pt.size = 0,same.y.lims=F,log = T,group.by="seurat_clusters")+theme(legend.position = "none",plot.title = element_blank())+xlab("")

FeaturePlot(clean_ACPC12,features = TE,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0.5) & NoAxes()
FeaturePlot(clean_ACPC12,features = PrE,pt.size = 0.05,order=T, ncol = 3,min.cutoff = 0.5) & NoAxes()
FeaturePlot(clean_ACPC12, features = EPI,pt.size = 0.05,ncol=3,order=T,min.cutoff = 0.5) & NoAxes()

######### DotPlot analysis #######
dot_features = c(EPI,PrE,TE)
as.data.frame(dot_features)[1]
DotPlot(clean_ACPC12,features = dot_features,group.by = c("type"),cluster.idents = F,assay = "RNA",scale =F,scale.min =0) + FontSize(7) + NoLegend() + scale_y_discrete(expand = expansion(mult = c(0.4, 0.4)))

######### Oct4/Cdx2 analysis######
clean_ACPC12$Cdx2_Pou5f1 = 0
data_subset <- GetAssayData(clean_ACPC12, layer = "counts")[c("Pou5f1","Cdx2"), , drop = FALSE]
clean_ACPC12$Cdx2_Pou5f1[data_subset["Pou5f1",]>1]=1
clean_ACPC12$Cdx2_Pou5f1[data_subset["Cdx2",]>1]=2
clean_ACPC12$Cdx2_Pou5f1[data_subset["Cdx2",]>1 & data_subset["Pou5f1",]>1]=3
DimPlot(clean_ACPC12, group.by = "Cdx2_Pou5f1",cols=c("lightgrey",pluri_color,te_color,"brown"),order=T) & NoAxes()
table(clean_ACPC12$Cdx2_Pou5f1)

####### Morula/twoC markers analysis ########
morula_markers=c("Pdzd3","Ifi27","Rab29","Oasl2","Sult5a1","Dnah2os")
twoC_markers = c("Zscan4c","Zscan4d","Zbed3","Omt2a","Obox1","Tcstv1")
FeaturePlot(clean_ACPC12, features = morula_markers,max.cutoff=3,min.cutoff = 0.5,ncol=3,order=T,cols=c("lightgrey",morula_color),keep.scale = "all") & NoAxes() & NoLegend()
FeaturePlot(clean_ACPC12, features = twoC_markers,max.cutoff=3,min.cutoff = 0.5,ncol=3,order=T,cols=c("lightgrey",twoc_color),keep.scale = "all") & NoAxes() & NoLegend()

######## Triangle analysis using average expression #####

TE = c("Cdx2","Id2","Elf5","Lgals1","Krt8","Eomes","Cldn7","Tfap2c","Gata3","Krt18")
PrE = c("Sox17","Gata6","Pdgfra","Dab2","Apoa1","Foxa2","Gata4","Serpinh1","Amn","Srgn")
pluri = c("Sox2","Pou5f1","Nanog")

calculate_avg_signature <- function(model,layer="counts")
{
  count_data <- GetAssayData(model, layer = layer)
  total_counts = sum(count_data)
  
  TE_matrix <- count_data[TE, ] 
  pluri_matrix <- count_data[pluri, ] 
  PrE_matrix <- count_data[PrE, ] 
  
  model  <- AddMetaData(model, metadata = colMeans(TE_matrix)*10E6/total_counts, col.name = "TE_avg")
  model  <- AddMetaData(model, metadata = colMeans(pluri_matrix)*10E6/total_counts, col.name = "pluri_avg")
  model  <- AddMetaData(model, metadata = colMeans(PrE_matrix)*10E6/total_counts, col.name = "PrE_avg")
  
  return(model)
}

clean_ACPC12<- calculate_avg_signature(clean_ACPC12, layer="counts")
columns = c("pluri_avg","TE_avg","PrE_avg")
df <- clean_ACPC12@meta.data[, columns]
df <- df / (rowSums(df[])+1e-5)
df$type <- clean_ACPC12$type
ggtern(data = df, aes(x = pluri_avg, y = TE_avg, z = PrE_avg, color=type)) + #, c("purple","blue","brown","darkgreen"))) +
  geom_point(size=1.5) +
  scale_color_manual(values = c("D_TE" = te_color, "A_EPI" = pluri_color, "C_PrE" = pre_color,"B_Early_PrE"=early_pre_color)) + 
  theme_minimal() +
  theme(legend.position = "none")  

