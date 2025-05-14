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

combined<-readRDS("combined_TFSEM_Nat_E85.rds")

##### preudo-bulk analysis, group each cluster and sample type #####
bulk <- AggregateExpression(combined, group.by = c("type", "seurat_clusters"), return.seurat = TRUE)
bulk_matrix = bulk@assays$RNA$data
df = as.data.frame(bulk_matrix)

##### For each cluster, calculate r-squared between TFSEM and InUt pseudo-bulk counts, assuming linear correlation #####
for (c in 0:33) {
  InUt_name = paste0("InUt_", c)
  x = df[[InUt_name]]
  TFSEM_name = paste0("TFSEM_", c)
  y = df[[TFSEM_name]]
  model <- lm(y ~ x)
  # Get the model summary and extract R²
  print(summary(model)$r.squared)
}

##### Print correlation scatter plot and r-squared for each cluster #####
cowplot::plot_grid(ncol = 5,
                   ggplot(df,aes(x=InUt_0,y=TFSEM_0)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.918", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_1,y=TFSEM_1)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.915", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_2,y=TFSEM_2)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.911", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_3,y=TFSEM_3)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.909", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_4,y=TFSEM_4)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.927", size = 4, color = "black"),
                   ggplot(df,aes(x=InUt_5,y=TFSEM_5)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.937", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_6,y=TFSEM_6)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.906", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_7,y=TFSEM_7)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.914", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_8,y=TFSEM_8)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.915", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_9,y=TFSEM_9)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.869", size = 4, color = "black"),
                   ggplot(df,aes(x=InUt_10,y=TFSEM_10)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.910", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_11,y=TFSEM_11)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.893", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_12,y=TFSEM_12)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.922", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_13,y=TFSEM_13)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.902", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_14,y=TFSEM_14)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.900", size = 4, color = "black"),
                   ggplot(df,aes(x=InUt_15,y=TFSEM_15)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.921", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_16,y=TFSEM_16)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.922", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_17,y=TFSEM_17)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.919", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_18,y=TFSEM_18)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.889", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_19,y=TFSEM_19)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.904", size = 4, color = "black"),
                   ggplot(df,aes(x=InUt_20,y=TFSEM_20)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.914", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_21,y=TFSEM_21)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.886", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_22,y=TFSEM_22)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.922", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_23,y=TFSEM_23)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.904", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_24,y=TFSEM_24)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.910", size = 4, color = "black"),
                   ggplot(df,aes(x=InUt_25,y=TFSEM_25)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.851", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_26,y=TFSEM_26)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.906", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_27,y=TFSEM_27)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.897", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_28,y=TFSEM_28)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.911", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_29,y=TFSEM_29)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.878", size = 4, color = "black"),
                   ggplot(df,aes(x=InUt_30,y=TFSEM_30)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.915", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_31,y=TFSEM_31)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.907", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_32,y=TFSEM_32)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.896", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_33,y=TFSEM_33)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+annotate("text", x = 4, y = 0.5, label = "0.818", size = 4, color = "black"), 
                   ggplot(df,aes(x=InUt_33,y=TFSEM_33)) + geom_point() + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") # space holder
)

