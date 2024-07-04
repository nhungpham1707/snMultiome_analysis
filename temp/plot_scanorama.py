import scanorama
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import matplotlib.pyplot as plt

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

sc.settings.set_figure_params(dpi=80)

save_dir = '/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/batchEffect/atac/scanorama/'
print("save dir is", save_dir)
save_name=save_dir+'atac'


with open(save_name+'_scanorama_tsne_umap.pkl', 'rb') as file:
    adata = pickle.load(file)



fig, axs = plt.subplots(2, 2, figsize=(8,8),constrained_layout=True)
sc.pl.tsne(adata, color="library", title="lib tsne", ax=axs[0,0], show=False)
sc.pl.umap(adata, color="library", title="lib umap", ax=axs[0,1], show=False)
sc.pl.tsne(adata, color="Subtype", title="Subtype tsne", ax=axs[1,0], show=False)
sc.pl.umap(adata, color="Subtype", title="Subtype umap", ax=axs[1,1], show=False)

fig.savefig(save_name+'_scanorama_umap.png')