import urllib.request
import anndata
import pandas as pd
from pathlib import Path
from ikarus import classifier, utils, data

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn import metrics


signatures_path = 'ikarus/signatures.gmt'
pd.read_csv(signatures_path, sep="\t", header=None)

# load pretrained model
model_path = Path("ikarus/core_model.joblib")
model = classifier.Ikarus(signatures_gmt=signatures_path, out_dir="ikarus")
model.load_core_model(model_path)

adata = anndata.read_h5ad("ikarus/lx049.h5ad")

# add gene_symbol metadata for model.predict 
adata.var['gene_symbol']= adata.var.index


_ = model.predict(adata, "lx049", save=True)

print('get umap on predicted model')
_ = model.get_umap(adata, "lx049", save=True)


print('loadn predicted result')
path = Path("ikarus/lx049")
results = pd.read_csv(path / "prediction.csv", index_col=0)
adata = anndata.read_h5ad(path / "adata_umap.h5ad")

print('plot confusion matrix')
y = adata.obs.loc[:, "tier_0"]
y_pred_lr = results["final_pred"]
acc = metrics.balanced_accuracy_score(y, y_pred_lr)
print(metrics.classification_report(y, y_pred_lr, labels=["Normal", "Tumor"]))
fig, ax = plot_confusion_matrix(
    y,
    y_pred_lr,
    classes=["Normal", "Tumor"],
    title=f"confusion matrix \n train on laughney+lee, test on lx049 \n balanced acc: {acc:.3f}",
)
fig.tight_layout()

fig.savefig('ikarus/lx049/matrix.png')

print('plot umap')
adata.obs.loc[:, "tier_0_pred_correctness"] = "wrongly assigned"
adata.obs.loc[
    adata.obs["tier_0"] == adata.obs["final_pred"],
    "tier_0_pred_correctness"
] = "correctly assigned"
adata.obs.loc[:, "tier_0_pred_wrong"] = pd.Categorical(
    adata.obs["tier_0"].copy(),
    categories=np.array(["Normal", "Tumor", "wrongly assigned"]),
    ordered=True
)
adata.obs.loc[
    adata.obs["tier_0_pred_correctness"] == "wrongly assigned",
    "tier_0_pred_wrong"
] = "wrongly assigned"

plt.rcParams["figure.figsize"] = [9, 6]

colors = [
    ["major"],
    ["tier_0_pred_wrong"]
    ]
titles = [
    ["major types"],
    ["prediction"]
    ]
palettes = [
    ["#7f7f7f", "#98df8a", "#aec7e8", "#9467bd", "#ff0000"],
    ["#aec7e8", "#ff0000", "#0b559f"], 
]
for color, title, palette in zip(colors, titles, palettes):
    ax = sc.pl.umap(
        adata, ncols=1, size=20, 
        color=color,
        title=title,
        wspace=0.25,
        vmax="p99",
        legend_fontsize=12,
        palette=palette,
        show=False
    )
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + 
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()

plt.savefig('ikarus/lx049/umap.png')

print('finished!')