import os
import tempfile

import matplotlib.pyplot as plt
import scanpy as sc

# Reproducibility
import scvi
from scvi.external import SysVI

import pickle 

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

save_dir = '/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/batchEffect/rna/sysvi/'
print("save dir is", save_dir)

# load data
adata = sc.read('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/sc_RNA/merge_all/rna_lbsb.h5ad')
print('adata is')
adata

# Setup adata for training
print('---------set up data for training------------')
SysVI.setup_anndata(
    adata=adata,
    batch_key="lbsb",
    categorical_covariate_keys=["library"],
)

# Initialise the model
print('-----------------train model-------------')
model = SysVI(adata=adata)
# Train
max_epochs = 50
model.train( plan_kwargs={
        "loss_weights": {
            "kl_weight": 0,
            "z_distance_cycle_weight": 0
            # Add additional parameters, such as number of epochs
        }
    },
    max_epochs=max_epochs,
    # Parameters used for checking losses
    log_every_n_steps=1,
    check_val_every_n_epoch=1,
    val_check_interval=1.0,
)

# save model 
with open(save_dir+'rnalbsb_50epo_kl0_model.pkl', 'wb') as file:
    pickle.dump(model, file)

# with open(save_dir+'rna_model.pkl', 'rb') as file:
#     model = pickle.load(file)

print(model)
# Plot loses
# The plotting code below was specifically adapted to the above-specified model and its training
# If changing the model or training the plotting functions may need to be adapted accordingly

# Make detailed plot after N epochs
print('---------------plot losses-----------')
epochs_detail_plot = 100
steps_detail_plot = epochs_detail_plot * int(
    model.trainer.logger.history["loss_validation"].shape[0] / max_epochs
)
detail_plot = epochs_detail_plot

# Losses to plot
losses = [
    "loss_train",
    "reconstruction_loss_train",
    "kl_local_train",
    "z_distance_cycle_train",
]
fig, axs = plt.subplots(2, len(losses), figsize=(len(losses) * 3, 4))
for ax_i, l_train in enumerate(losses):
    l_val = l_train.replace("_train", "_validation")
    l_name = l_train.replace("_train", "")
    # Change idx of epochs to start with 1 so that below adjustment when
    # train on step which only works for val leads to appropriate multiplication
    l_val_values = model.trainer.logger.history[l_val].copy()
    l_val_values.index = l_val_values.index + 1
    l_train_values = model.trainer.logger.history[l_train].copy()
    l_train_values.index = l_train_values.index + 1
    # This happens if log on step as currently this works only for val loss
    if l_train_values.shape[0] < l_val_values.shape[0]:
        l_train_values.index = l_train_values.index * int(
            l_val_values.shape[0] / l_train_values.shape[0]
        )
    for l_values, c, alpha, dp in [
        # Train loss logged on epoch in either case now
        (l_train_values, "tab:blue", 1, epochs_detail_plot),
        (l_val_values, "tab:orange", 0.5, detail_plot),
    ]:
        axs[0, ax_i].plot(l_values.index, l_values.values.ravel(), c=c, alpha=alpha)
        axs[0, ax_i].set_title(l_name)
        axs[1, ax_i].plot(
            l_values.index[dp:], l_values.values.ravel()[dp:], c=c, alpha=alpha
        )

# fig.tight_layout()
fig.savefig(save_dir+'rnalbsb_50epo_kl0_losses.png')

# integrated embedding
print('----integrate embedding-----')

# Get embedding - save it into X of new AnnData
embed = model.get_latent_representation(adata=adata)
embed = sc.AnnData(embed, obs=adata.obs)

unique_patients = embed.obs["lbsb"].unique()
patient_mapping = dict(zip(range(0, len(unique_patients)), unique_patients))


# Make system categorical for plotting below
embed.obs["lbsb"] = embed.obs["lbsb"].map(patient_mapping)

# Compute UMAP
print('----compute umap------')
sc.pp.neighbors(embed, use_rep="X")
sc.tl.umap(embed)
with open(save_dir+'rnalbsb_50epo_kl0_embed.pkl', 'wb') as file:
    pickle.dump(embed, file)

# Plot UMAP embedding

# Obs columns to color by
cols = ["Individual.ID", "singler.pruned.labels", "Subtype", 'lbsb']

# One plot per obs column used for coloring
print('-----plot umap----')
fig, axs = plt.subplots(len(cols), 1, figsize=(3, 3 * len(cols)))
for col, ax in zip(cols, axs):
    sc.pl.embedding(
        embed,
        "X_umap",
        color=col,
        s=10,
        ax=ax,
        show=False,
        sort_order=False,
        frameon=False,
    )
fig.savefig(save_dir+'lbsb_50epo_kl0_umap.png')

print('finished!')