#######################################################
# Written by Dr. Noa Novershtern & Nathan Levy 
# https://github.com/hannalab/Mouse_TFSEM
# This code trains with query data (transfer learning)
#######################################################
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvi # type: ignore
import sys
import anndata as ad
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
import json
import os

args = sys.argv

with open(args[1], "r") as file:
    config = json.load(file)
    
######## Read parameters from a config file ########
WORKING_DIR = config["WORKING_DIR"]
OUTPUT_DIR = config["OUTPUT_DIR"]
MODEL_FILE = config["MODEL_FILE"]
ADATA_FILE = config["ADATA_FILE"]
TEST_FILE = config["TEST_FILE"]
TEST_PREDICTIONS = int(config["TEST_PREDICTIONS"])
EPOCHS_TEST = int(config["EPOCHS_TEST"])
QUERY_ANNOTATIONS = config["QUERY_ANNOTATIONS"]
CHECK_VAL_EVERY_EPOCH = int(config["CHECK_VAL_EVERY_EPOCH"])

if (not os.path.exists(OUTPUT_DIR)):
  os.mkdir(OUTPUT_DIR)
if (not os.path.exists(OUTPUT_DIR+"/figures")):
   os.mkdir(OUTPUT_DIR+"/figures")
   
    
    
######### Read train and query data ###########    
print("reading data..............................")

adata_query = sc.read(TEST_FILE)
adata = ad.read_h5ad(ADATA_FILE)
lvae = scvi.model.SCANVI.load(MODEL_FILE,
                              adata = adata)

adata_query.layers["counts"] = adata_query.X
adata_query = adata_query[:, adata_query.var_names.isin(adata.var_names)].copy()

######### Define model and train with query data ###########    
print("define model with query data ..................................")
scvi.model.SCANVI.prepare_query_anndata(adata_query, lvae)

vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    lvae,
)

print("train model with query data ..................................")
vae_q.train(
    max_epochs=EPOCHS_TEST,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=CHECK_VAL_EVERY_EPOCH,
)

######## Record loss ###########
vae_q.history["kl_local_validation"].plot()
plt.savefig(
    OUTPUT_DIR + "figures/scVI_kl_local_validation_query.png", dpi=300, bbox_inches="tight"
)
vae_q.history["reconstruction_loss_validation"].plot()
plt.savefig(
    OUTPUT_DIR + "figures/scVI_recon_loss_validation_query.png", dpi=300, bbox_inches="tight"
)

######## Get and record predictions #########
print("get predictions ..................................")
adata_query.obsm["X_scanVI"] = vae_q.get_latent_representation()
adata_query.obs["predictions"] = vae_q.predict()

df = adata_query.obs.groupby([QUERY_ANNOTATIONS, "predictions"]).size().unstack(fill_value=0)
df.to_csv(OUTPUT_DIR+"prediction_matrix.csv")

norm_df = df.T/df.T.sum(axis=0)

plt.figure(figsize=(12, 12))
_ = plt.pcolor(norm_df, edgecolors='grey', linewidths=1)
_ = plt.xticks(np.arange(0.5, len(df.T.columns), 1), df.T.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.T.index), 1), df.T.index)
plt.xlabel("Predicted")
plt.ylabel("Clusters")

plt.savefig(OUTPUT_DIR + "figures/predictions.png", dpi=300, bbox_inches="tight")



if (TEST_PREDICTIONS):
    predictions = adata_query.obs["predictions"].tolist()
    true_labels = adata_query.obs[QUERY_ANNOTATIONS].tolist()
    print("accuracy score.................")
    print(accuracy_score(true_labels,predictions))
    print(classification_report(true_labels, predictions, digits=4))

print("saving predictions to file............................................")
adata_query.obs["predictions"].to_csv(OUTPUT_DIR+"predictions.tsv",sep="\t")
adata_query.obs_names.to_csv(OUTPUT_DIR+"cell_names.tsv",sep="\t")


######## plot data with query #########
print("ploting data with query..............................................")
adata_full = ad.concat([adata,adata_query], join='outer')
adata_full.obsm['X_scanVI'] = np.vstack([adata.obsm['X_scanVI'], adata_query.obsm['X_scanVI']])
neighbors = sc.pp.neighbors(adata_full, use_rep="X_scanVI")
sc.tl.umap(adata_full, min_dist=0.1)

sc.pl.embedding(
    adata_full,
    basis="X_umap",
    color=["dataset"],
    frameon=False,
    ncols=1,
)
plt.savefig(
    OUTPUT_DIR + "figures/umap_data_query_batch.png", dpi=300, bbox_inches="tight"
)


print("...........................................................DONE")

