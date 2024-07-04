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

# adata = sc.read('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/sc_RNA/merge_all/rna.h5ad')
# adata_ori = adata
# adata2 = sc.AnnData(X=adata.raw.X, var=adata.raw.var, obs = adata.obs)

# sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
# sc.pp.log1p(adata2)
# #variable genes for the full dataset
# sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5)

# sc.pl.highly_variable_genes(adata2)

# print("Highly variable genes: %d"%sum(adata2.var.highly_variable))

# var_genes_all = adata2.var.highly_variable
# sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'library')

# print("Highly variable genes intersection: %d"%sum(adata2.var.highly_variable_intersection))

# print("Number of batches where gene is variable:")
# print(adata2.var.highly_variable_nbatches.value_counts())

# var_genes_batch = adata2.var.highly_variable_nbatches > 0

# print("Any batch var genes: %d"%sum(var_genes_batch))
# print("All data var genes: %d"%sum(var_genes_all))
# print("Overlap: %d"%sum(var_genes_batch & var_genes_all))
# print("Variable genes in all batches: %d"%sum(adata2.var.highly_variable_nbatches ==3))
# print("Overlap batch instersection and all: %d"%sum(var_genes_all & adata2.var.highly_variable_intersection))

# var_select = adata2.var.highly_variable_nbatches > 1
# var_genes = var_select.index[var_select]
# len(var_genes)

# # split per batch into new objects.
# batches = adata2.obs['library'].unique()
# alldata = {}
# for batch in batches:
#     alldata[batch] = adata2[adata2.obs['library'] == batch,]

# alldata    

# #subset the individual dataset to the same variable genes as in MNN-correct.
# alldata2 = dict()
# for ds in alldata.keys():
#     print(ds)
#     alldata2[ds] = alldata[ds][:,var_genes]

# #convert to list of AnnData objects
# adatas = list(alldata2.values())

# # run scanorama.integrate
# scanorama  = scanorama.integrate_scanpy(adatas, dimred = 50,)
# x_scanoramas =[anndata_obj.obsm['X_scanorama'] for anndata_obj in adatas]
# all_s = np.concatenate(x_scanoramas)
# print(all_s.shape)
# # add to the AnnData object
# adata.obsm["SC"] = all_s

# with open(save_name+'_scanorama.pkl', 'wb') as file:
#     pickle.dump(adata, file)


# # tsne and umap
# sc.pp.neighbors(adata, n_pcs =50, use_rep = "SC")
# sc.tl.umap(adata)
# sc.tl.tsne(adata, n_pcs = 50, use_rep = "SC")

# with open(save_name+'_scanorama_tsne_umap.pkl', 'wb') as file:
#     pickle.dump(adata, file)

with open(save_name+'_scanorama_tsne_umap.pkl', 'rb') as file:
    adata = pickle.load(file)


adata.write_h5ad(save_name+'_scanorama')

fig, axs = plt.subplots(2, 2, figsize=(8,8),constrained_layout=True)
sc.pl.tsne(adata, color="library", title="lib tsne", ax=axs[0,0], show=False)
sc.pl.umap(adata, color="library", title="lib umap", ax=axs[0,1], show=False)
sc.pl.tsne(adata, color="Subtype", title="Subtype tsne", ax=axs[1,0], show=False)
sc.pl.umap(adata, color="Subtype", title="Subtype umap", ax=axs[1,1], show=False)

sc.pl.umap(adata_ori, color="Subtype", title="Subtype umap")

fig.savefig(save_name+'_scanorama_umap.png')