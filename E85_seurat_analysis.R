####################################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
# Run on R.4.3.1
#####################################
#WORKING_DIR = "<enter your working dir here>" 
#setwd(WORKING_DIR)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)
library(cowplot)
set.seed(1234)

######define my colors ########
# neural
cluster5_col = "#E58700"
cluster26_col = "#FFCC99"
# cardiac
cluster15_col="#FFFF00"
# mesoderm
cluster2_col = "#33CC00"
cluster4_col = "#6BB100"
cluster6_col = "#99CC99"
cluster8_col = "#33FF00"
cluster17_col = "#CCCC00"
cluster19_col = "#00BA38"
cluster22_col = "#A3A500"
cluster23_col = "#99FF00"
# Endoderm
cluster0_col = "#F8766D"
cluster20_col = "#CC6666"
cluster30_col = "brown1"
cluster31_col = "#FF0000"
cluster32_col = "red3"
cluster12_col = "pink3"
cluster14_col = "#660000"
cluster16_col="pink"
cluster18_col = "#330000"
cluster27_col = "#990000"
cluster33_col = "brown"
# endothelial
cluster7_col = "#B983FF"
cluster28_col = "#663399"
cluster29_col = "#9966CC"
# Chorion
cluster9_col="#CC00CC"
cluster10_col="#FF33CC"
cluster21_col = "#FFCCFF"
cluster25_col="violet"
# Blood
cluster1_col = "#619CFF"
cluster3_col = "#0066FF"
cluster11_col = "#00B0F6"
cluster13_col = "#33FFFF"
cluster24_col = "blue"

my_colors_ct = c(cluster0_col,cluster1_col,cluster2_col,cluster3_col,cluster4_col,cluster5_col,cluster6_col,cluster7_col,cluster8_col,cluster9_col,cluster10_col,cluster11_col,cluster12_col,cluster13_col,cluster14_col,cluster15_col,cluster16_col,cluster17_col,
                 cluster18_col,cluster19_col,cluster20_col,cluster21_col,cluster22_col,cluster23_col,cluster24_col,cluster25_col,cluster26_col,cluster27_col,cluster28_col,cluster29_col,cluster30_col,cluster31_col,cluster32_col,cluster33_col)

##### Load data after CellRanger #####
TFSEM_vs_Nat_E85.data <- Read10X(data.dir = "<path to data>/outs/count/filtered_feature_bc_matrix/")
TFSEM.data <- CreateSeuratObject(counts = TFSEM_vs_Nat_E85.data, project = "TFSEM", min.cells = 3, min.features = 200)

#####Add sample names and type ######
groups<-as.numeric(str_extract(colnames(TFSEM.data),"\\d+"))
sample <- groups
sample<-ifelse(groups==1,"n_TFSEM_E8_5",sample)
sample<-ifelse(groups==2,"o_TFSEM_E8_5-1EMB",sample)
sample<-ifelse(groups==3,"p_TFSEM_E8_5-5EMB",sample)
sample<-ifelse(groups==4,"q_TFSEM_E8_5_V2",sample)
sample<-ifelse(groups==5,"a_InUtE85singleembryo1",sample)
sample<-ifelse(groups==6,"b_InUtE85singleembryo2",sample)
sample<-ifelse(groups==7,"c_InUtE85singleembryo4",sample)
sample<-ifelse(groups==8,"d_InUtE85twoembryos",sample)
sample<-ifelse(groups==9,"i_Cdx2-0d-2pooled",sample)
sample<-ifelse(groups==10,"j_Cdx2-10d-2pooled",sample)
sample<-ifelse(groups==11,"k_Cds2-3d-1emb",sample)
sample<-ifelse(groups==12,"l_Cds2-3d-2pooled",sample)
sample<-ifelse(groups==13,"m_TSCbl-2pooled",sample)
sample<-ifelse(groups==14,"g_ExUt_E85_1",sample)
sample<-ifelse(groups==15,"h_ExUt_E85_pool",sample)
TFSEM.data$sample<-sample

type <- groups
type<-ifelse(groups %in% c(1,2,3,4),"TFSEM",type)
type<-ifelse(groups %in% c(5,6,7,8),"InUt",type)
type<-ifelse(groups %in% c(9,10,11,12,13),"SEM",type)
type<-ifelse(groups %in% c(14,15),"ExUt",type)
TFSEM.data$type<-type

###### QC ###### 
TFSEM.data[["percent.mt"]] <- PercentageFeatureSet(TFSEM.data, pattern = "^mt-")
VlnPlot(TFSEM.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,same.y.lims=F,log = T,split.by = "sample",group.by="sample")+theme(legend.position = "none")+xlab("")
plot1 <- FeatureScatter(TFSEM.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TFSEM.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

##### Filter ######
TFSEM.data_sub <- subset(TFSEM.data, subset = nFeature_RNA > 1000 & nFeature_RNA <8000 & percent.mt < 10 & nCount_RNA>3000)
rm(TFSEM.data)

##### Normalization ######
new_sub.list <- SplitObject(TFSEM.data_sub, split.by = "sample")
new_sub.list <- lapply(X = new_sub.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selesampleion.method = "vst", nfeatures = 2000)
})

##### Combined UMAP analysis, with anchoring ######
dims=1:15
anchors <- FindIntegrationAnchors(object.list = new_sub.list, dims = dims,anchor.features = 2000)
combined <- IntegrateData(anchorset = anchors, dims = dims)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = dims,min.dist = 0.1)
combined <- FindNeighbors(combined, reduction = "pca", dims = dims)
combined <- FindClusters(combined, resolution = 0.8)

#saveRDS(combined,"combined_TFSEM_Nat_E85.rds")

###### Visualization ######
# To upload the exact seurat object used in the paper: 
# combined <- readRDS("combined_TFSEM_Nat_E85.rds")
combined$seurat_clusters = combined$integrated_snn_res.0.8
DimPlot(combined, reduction = "umap", group.by = "seurat_clusters",label=T,pt.size = 1,cols = my_colors_ct)
DimPlot(combined, cols=my_colors_ct,reduction = "umap", group.by="seurat_clusters",split.by="type",label=T,pt.size = 1)& NoAxes()
VlnPlot(combined, features = c("nFeature_RNA"), ncol = 1,pt.size = 0,same.y.lims=F,log = T,group.by="sample")+theme(legend.position = "none",plot.title = element_blank())+xlab("")
FeaturePlot(combined,"nCount_RNA")
FeaturePlot(combined,"nFeature_RNA")

###### Find markers ######
TFSEM_85_markers = FindAllMarkers(combined,logfc.threshold = 0.5,only.pos =T )
write.csv(TFSEM_85_markers,"TFSEM85_huge_markers.csv")

###### Dotplot analysis ########
endoderm_dotblot = c("Apoa1","Apoa2","Ttr","Afp")
exend_dotblot = c("Lama1","Pdgfra")
hindgut_dotblot=c("Sox17", "Klf5","Cyp26a1")
midgut_dotblot=c("Nepn","Shh","Onecut2")
foregut_dotblot = c("Tbx1","Otx2","Trh")
blood_dotblot = c("Hbb-y","Hba-x","Hba-a2")
hpc_dotblot = c("Cd33","Ptprc","Runx1")
endo_dotblot = c("Cd34","Cdh5","Kdr")
brain_dotblot = c("Pax6","Six3","Sox2")
neural_dotblot = c("Pax3","Fgf15","Tfap2a")
cardiac_dotblot = c("Myh6","Ttn","Actc1")
placodes_dotblot = c("Utf1","Rhox9","Tfap2c")
allantois_dotblot=c("Tbx4","Vcam1","Rspo3")
somites_dotblot = c("Foxc2","Tbx18","Meox1")
exem_dotbplot = c("Bnc2","Col5a2","Lum")
meso_dotblot = c("T","Isl1","Hoxa9")
pgc_dotblot = c("Dppa3","Prdm1","Pou5f1")
tailbud_dotblot = c("Fgf8","Cdx2")
notochord_dotblot = c("Noto")
dotblot_features = c(endoderm_dotblot,exend_dotblot,hindgut_dotblot,midgut_dotblot,foregut_dotblot,
                     blood_dotblot,hpc_dotblot,endo_dotblot,tailbud_dotblot,brain_dotblot,neural_dotblot,
                     cardiac_dotblot,placodes_dotblot,exem_dotbplot,meso_dotblot,allantois_dotblot,somites_dotblot)

DotPlot(combined,features = dotblot_features,group.by = c("order"),cluster.idents = F,assay = "RNA",scale =T,scale.min =0.5) + 
  scale_y_discrete(expand = expansion(mult = c(0.04, 0.04))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size =8))


######### EPC analysis ######
TGC_progenitors =c("Adm","Fosl1","Hand1","Trpm5","Maged2","Prl5a1")
Uncommited_EPC =c("Chsy1","Gjb3","Krt19","Lgals1","Cald1","Ctsl")
Intermediate_Chorion =c("Ascl2","Fgfr2","Cited1","Gjb3","Ndrg1","Irx2","Irx3")
Chorion_progenitors =c("Sox3","Dusp6","Nat8l","Bmp4","Sox2","Esrrb","Eomes")
general_placenta =c("Rhox5")
SpA_TGC =c("Ctla2a","Pecam1","Ramp3","Igfbp7","Nos3")
SpT_Gly =c("Dlx3","Car2","Ncam1","Pcdh12","Tpbpa")
p_TGC =c("Star","Serpinb9d","Hsd3b6","Rhox6","Cts7","Prl2c3")
Chorion =c("Irx4","Esx1","Id1","Id3","Phlda2","Klhl13")

combined$TGC_progenitors_mean <- rowMeans(FetchData(combined, vars = TGC_progenitors))
combined$Uncommited_EPC_mean <- rowMeans(FetchData(combined, vars = Uncommited_EPC))
combined$Intermediate_Chorion_mean <- rowMeans(FetchData(combined, vars = Intermediate_Chorion))
combined$Chorion_mean <- rowMeans(FetchData(combined, vars = Chorion))
combined$Chorion_progenitors_mean <- rowMeans(FetchData(combined, vars = Chorion_progenitors))
combined$General_placenta_mean <- rowMeans(FetchData(combined, vars = general_placenta))
combined$SpA_TGC_mean <- rowMeans(FetchData(combined, vars = SpA_TGC))
combined$p_TGC_mean <- rowMeans(FetchData(combined, vars = p_TGC))
combined$SpT_Gly_mean <- rowMeans(FetchData(combined, vars = SpT_Gly))

par(mfrow = c(3, 3))
EPC_groups_mean = c("Uncommited_EPC_mean","General_placenta_mean","Chorion_progenitors_mean","Intermediate_Chorion_mean",
                    "Chorion_mean","TGC_progenitors_mean","SpA_TGC_mean","p_TGC_mean","SpT_Gly_mean")

FeaturePlot(combined,features = EPC_groups_mean,cols=c("grey","red","brown3","brown"),min.cutoff = 0) & NoAxes() & NoLegend()

####### EPC barplot analysis ########
# Count how many of the placental cells (Rhox5+) express 
# the majority of the gene markers of each cell type.
#######
count_data<- GetAssayData(TFSEM.data_sub,layer="counts")
row_index <- which(rownames(count_data) == "Rhox5")
selected_columns <- count_data[row_index, ] > 0
EPC <- count_data[, selected_columns]

TGC_progenitors_matrix <- EPC[TGC_progenitors,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(TGC_progenitors_matrix>0)>dim(TGC_progenitors_matrix)[1]/2, col.name = "TGC_progenitors")
Uncommited_EPC_matrix <- EPC[Uncommited_EPC,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(Uncommited_EPC_matrix>0)>dim(Uncommited_EPC_matrix)[1]/2, col.name = "Uncommited_EPC")
Intermediate_Chorion_matrix <- EPC[Intermediate_Chorion,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(Intermediate_Chorion_matrix>0)>dim(Intermediate_Chorion_matrix)[1]/2, col.name = "Intermediate_Chorion")
Chorion_progenitors_matrix <- EPC[Chorion_progenitors,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(Chorion_progenitors_matrix>0)>dim(Chorion_progenitors_matrix)[1]/2, col.name = "Chorion_progenitors")
Chorion_progenitors_matrix <- EPC[Chorion_progenitors,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(Chorion_progenitors_matrix>0)>dim(Chorion_progenitors_matrix)[1]/2, col.name = "Chorion_progenitors")
SpA_TGC_matrix <- EPC[SpA_TGC,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(SpA_TGC_matrix>0)>dim(SpA_TGC_matrix)[1]/2, col.name = "SpA_TGC")
SpT_Gly_matrix <- EPC[SpT_Gly,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(SpT_Gly_matrix>0)>dim(SpT_Gly_matrix)[1]/2, col.name = "SpT_Gly")
p_TGC_matrix <- EPC[p_TGC,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(p_TGC_matrix>0)>dim(p_TGC_matrix)[1]/2, col.name = "p_TGC")
Chorion_matrix <- EPC[Chorion,]
TFSEM.data_sub <- AddMetaData(TFSEM.data_sub, metadata = colSums(Chorion_matrix>0)>dim(Chorion_matrix)[1]/2, col.name = "Chorion")

table(TFSEM.data_sub$TGC_progenitors,TFSEM.data_sub$type)
table(TFSEM.data_sub$Uncommited_EPC,TFSEM.data_sub$type)
table(TFSEM.data_sub$Intermediate_Chorion,TFSEM.data_sub$type)
table(TFSEM.data_sub$Chorion_progenitors,TFSEM.data_sub$type)
table(TFSEM.data_sub$SpA_TGC,TFSEM.data_sub$type)
table(TFSEM.data_sub$SpT_Gly,TFSEM.data_sub$type)
table(TFSEM.data_sub$p_TGC,TFSEM.data_sub$type)
table(TFSEM.data_sub$Chorion,TFSEM.data_sub$type)


####### TGC markers in chorion analysis #########
TGC_markers = c("Prl2c2","Prl3d1","Prl4a1","Cdh5")
chorion_clusters <- subset(combined, subset = (seurat_clusters==9|seurat_clusters==10))
FeaturePlot(chorion_clusters, features = TGC_markers,min.cutoff = 0,split.by = "type",order=T,combine = T) & NoAxes() & NoLegend()

###### PGC analysis #########
PGC = c("Pou5f1","Nanog","Klf4","Sox2","Dppa3","Prdm1","Bmp4","Prdm14","Alpl","Nanos3")
cluster21 = subset(combined, subset = (seurat_clusters==21))
cluster21$PGC_mean <- rowMeans(FetchData(cluster21, vars = PGC))
FeaturePlot(cluster21,features="PGC_mean",min.cutoff = 0,order=T,cols=c("grey","red","brown3","brown"))#& NoLegend() & NoAxes() #


###### DimPlot TFSEM vs InUtero / vs SEM #########
combined$fig = "NO"
combined$fig[combined$type %in% c("TFSEM","InUt")] = "YES"
DimPlot(combined, reduction = "umap", group.by="type",split.by="fig",label=F,raster=T,order="InUt",label.size = 2,cols=c("grey","black","pink","brown2"))& NoAxes() 

combined$fig = "NO"
combined$fig[combined$type %in% c("TFSEM","SEM")] = "YES"
DimPlot(combined, reduction = "umap", group.by="type",split.by="fig",label=F,raster=T,order="InUt",label.size = 2,cols=c("black","brown2","grey","pink"))& NoAxes() 

