##############################################
# Written by Dr. Noa Novershtern & Nathan Levy 
# https://github.com/hannalab/Mouse_TFSEM
# This code trains scVI and scANVI models
##############################################
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvi  # type: ignore
from pprint import pprint
from utils import plot_dist
import sys
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
import json
import os

args = sys.argv

######## Read parameters from a config file ########
with open(args[1], "r") as file:
    config = json.load(file)
    
WORKING_DIR = config["WORKING_DIR"]
OUTPUT_DIR = config["OUTPUT_DIR"]

TRAIN_FILE = config["TRAIN_FILE"]

if (not os.path.exists(OUTPUT_DIR)):
  os.mkdir(OUTPUT_DIR)
if (not os.path.exists(OUTPUT_DIR+"/figures")):
   os.mkdir(OUTPUT_DIR+"/figures")
   
N_HIDDEN = int(config["N_HIDDEN"])
N_LAYERS = int(config["N_LAYERS"])
N_LATENT = int(config["N_LATENT"])
EPOCHS_TRAIN = int(config["EPOCHS_TRAIN"])
EPOCHS_TEST = int(config["EPOCHS_TEST"])
CHECK_VAL_EVERY_EPOCH = int(config["CHECK_VAL_EVERY_EPOCH"])
LR = float(config["LR"])
SCANVI_LABELS_KEY = config["SCANVI_LABELS_KEY"]
BATCH_CATEGORY = config["BATCH_CATEGORY"]
    
######### Read train data ###########
print("Reading data..............................")
adata = sc.read(TRAIN_FILE)

######### Process and filter train data ###########
adata.layers["counts"] = adata.X
adata.obs["total_counts"] = adata.layers["counts"].sum(axis=1)

print("Filtering data...............................")
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_counts=5)
sc.pp.filter_genes(adata, max_counts=100000)
sc.pp.filter_cells(adata, max_counts=100000)

print(adata)

######### Plot count distribution ###########
print("Plotting count distribution...............................")
plot_dist(adata, group_by=BATCH_CATEGORY, query="total_counts")
plt.savefig(OUTPUT_DIR + "figures/train_counts_per_batch.png", dpi=300, bbox_inches="tight")

######### Normalize ###########
print("Normalizing...............................")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

######### Select highly variable genes ###########
print("Selecting high var genes ...............................")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="cell_ranger", subset=True)
adata.var["highly_variable"].to_csv(OUTPUT_DIR + "train_top_var_genes.csv")

######### PCA ###########
print("PCA..................................")
sc.pp.pca(adata, use_highly_variable=True)

######## UMAP ###########
print("Calculating neighbors..................................")
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)

print("Calculating umap..................................")
sc.tl.umap(adata, min_dist=0.3)

sc.pl.umap(
    adata,
    color=["new_celltype"],
    frameon=False,
)
plt.savefig(
    OUTPUT_DIR + "figures/umap_before_learning_cell_type.png", dpi=300, bbox_inches="tight"
)

sc.pl.umap(adata, color=[BATCH_CATEGORY], ncols=2, frameon=False)
plt.savefig(
    OUTPUT_DIR + "figures/umap_before_learning_batch.png", dpi=300, bbox_inches="tight"
)

######## Define an scVI model ###########
print("Defining and training scvi model..................................")
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH_CATEGORY)
vae = scvi.model.SCVI(adata, n_hidden=N_HIDDEN, n_layers=N_LAYERS, n_latent=N_LATENT)

######## Train the model ###########
vae.train(
    max_epochs=EPOCHS_TRAIN,
    check_val_every_n_epoch=CHECK_VAL_EVERY_EPOCH,
    plan_kwargs={"lr": LR},
    early_stopping=True, # will stop if validation loss is not improving 
)

with open(OUTPUT_DIR+"scvi_history.txt", "w") as f:
    pprint(vae.history, stream=f)
    
print("Finished_training......scvi ELBO is", vae.get_elbo().item())

######## Record loss ###########
vae.history["kl_local_validation"].plot()
plt.savefig(OUTPUT_DIR + "figures/scVI_kl_local_validation.png", dpi=300, bbox_inches="tight")
vae.history["reconstruction_loss_validation"].plot()
plt.savefig(
    OUTPUT_DIR + "figures/scVI_recon_loss_validation.png", dpi=300, bbox_inches="tight"
)

######## Present embeddings on UMAP ###########
print("Umap after scvi learning..................................")
adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.1)

sc.pl.embedding(
    adata,
    basis="X_umap",
    color=["new_celltype"],
    frameon=False,
    ncols=1,
)
plt.savefig(
    OUTPUT_DIR + "figures/scVI_umap_after_learning_cell_type.png", dpi=300, bbox_inches="tight"
)

sc.pl.embedding(
    adata,
    basis="X_umap",
    color=[BATCH_CATEGORY],
    frameon=False,
    ncols=1,
)
plt.savefig(
    OUTPUT_DIR + "figures/scVI_umap_after_learning_batch.png", dpi=300, bbox_inches="tight"
)

######## Save scVI model ###########
print("Saving vae model..................................")
#adata.write_h5ad(OUTPUT_DIR + "trainData_after_scVI.h5ad")
vae.save(dir_path=OUTPUT_DIR+"scVI_model",save_anndata=True,)


######## Define and train scANVI model ###########
print("Defining and training scanvi model..................................")
adata.obs[SCANVI_LABELS_KEY] = adata.obs["new_celltype"].values.tolist()
adata.obs["labels_scanvi"].fillna("NaN", inplace=True)

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key=SCANVI_LABELS_KEY,
    unlabeled_category="NaN",
)

lvae.train(
    max_epochs=EPOCHS_TRAIN,
    check_val_every_n_epoch=CHECK_VAL_EVERY_EPOCH,
    plan_kwargs={"lr": LR},
)

with open(OUTPUT_DIR+"scanvi_history.txt", "w") as f:
    pprint(lvae.history, stream=f)
    
print("finished_training......sancvi ELBO is", lvae.get_elbo().item())

######## Record loss ###########
lvae.history["kl_local_validation"].plot()
plt.savefig(OUTPUT_DIR + "figures/kl_local_val_after_scanVI.png")

lvae.history["reconstruction_loss_validation"].plot()
plt.savefig(OUTPUT_DIR + "figures/recon_loss_val_after_scanVI.png")


######## Present embeddings on UMAP ###########
adata.obsm["X_scanVI"] = lvae.get_latent_representation()

print("Umap after scanvi training..................................")
sc.pp.neighbors(adata, use_rep="X_scanVI")
sc.tl.umap(adata, min_dist=0.1)

sc.pl.embedding(
    adata,
    basis="X_umap",
    color=["new_celltype"],
    frameon=False,
    ncols=1,
)
plt.savefig(
    OUTPUT_DIR + "figures/scanVI_umap_after_learning_cell_type.png",
    dpi=300,
    bbox_inches="tight",
)

sc.pl.embedding(
    adata,
    basis="X_umap",
    color=["new_celltype"],
    frameon=False,
    ncols=1,
    legend_loc = 'on data',
)

plt.savefig(
    OUTPUT_DIR + "figures/scanVI_umap_after_learning_cell_type_w_labels.png",
    dpi=300,
    bbox_inches="tight",
)

sc.pl.embedding(
    adata,
    basis="X_umap",
    color=[BATCH_CATEGORY],
    frameon=False,
    ncols=1,
)
plt.savefig(
    OUTPUT_DIR + "figures/scanVI_umap_after_learning_batch.png", dpi=300, bbox_inches="tight"
)

######## Save scVI model ###########
print("saving scanvi model ..................................")
lvae.save(dir_path=OUTPUT_DIR+"lvea_model",save_anndata=True,)


print("...........................................................DONE")

