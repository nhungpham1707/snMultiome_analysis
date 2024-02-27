import os
from scipy.io import mmread
import pandas as pd
import anndata
import scanpy as sc
import numpy as np
from pybedtools import BedTool
from scipy.sparse import coo_matrix

data_dir = /hpc/pmc_drost/PROJECTS/cell_origin_NP/data/analyses/analyses
#make a master peak set
!cat data/atac_v1_pbmc_10k_peaks.bed data/atac_pbmc_10k_nextgem_peaks.bed | bedtools sort | bedtools merge > data/masterpeaks.bed


atac_peaks.bed    
# load the master peaks as BedTool object 
peak =pd.read_csv('data/masterpeaks.bed', sep='\t', header=None)
peak.columns=['chrom','start','end']

peak.loc[:, "ridx"] = range(peak.shape[0])
# remove sex chroms and chrom M
dfpeak = peak[~peak.chrom.isin(['chrX', 'chrY','chrM'])].copy()
dfpeak.loc[:,'idx'] = dfpeak.apply(lambda row: f'{row.chrom}:{row.start}-{row.end}', axis=1)
dfpeak.set_index('idx', inplace=True)
peak = BedTool.from_dataframe(dfpeak)

batches =  {
            'pbmc_10k': {'frag': 'data/atac_v1_pbmc_10k_fragments.tsv.gz',
                         'cells': 'data/atac_v1_pbmc_10k_singlecell.csv'},
            'pbmc_10k_nextgem': {'frag': 'data/atac_pbmc_10k_nextgem_fragments.tsv.gz',
                         'cells': 'data/atac_pbmc_10k_nextgem_singlecell.csv'},
           }

def get_fragment_bedtool(fragments, barcodes):
    """ Load and filter fragments
    
    Only valid cells (defined by is__cell_barcode) are used.
    """
    df = pd.read_csv(fragments,sep='\t', header=None)
    df.columns = ['chr','start','end','barcode', 'count']
    #bcf = pd.read_csv(keepbarcodes)
    barcodes = barcodes[barcodes.is__cell_barcode==1]
    barcodes.loc[:,'idx'] = range(barcodes.shape[0])
    df = pd.merge(df, barcodes, on='barcode', how='inner')[['chr','start','end','barcode', 'idx']]
    return BedTool.from_dataframe(df)

adatas = []
for batchname in batches:
    barcodes = pd.read_csv(batches[batchname]['cells'])
    barcodes = barcodes[barcodes.is__cell_barcode==1]
    barcodes.loc[:,"batch"] = batchname
    barcodes.set_index('barcode', inplace=True)
    
    frags = get_fragment_bedtool(batches[batchname]['frag'], barcodes)
    
    peakcounts = peak.intersect(frags,
                             wa=True,
                             wb=True).to_dataframe()
    sparse_data = np.asarray([np.ones(peakcounts.shape[0]),
                             peakcounts.name, peakcounts.itemRgb]).T
    sparse_data = np.unique(sparse_data, axis=0)
    mat = coo_matrix((sparse_data[:,0], (sparse_data[:,1], sparse_data[:,2])),
                     shape=(len(peak), len(barcodes)))
    adata = anndata.AnnData(mat.T.tocsr(), obs=barcodes, var=dfpeak)
    adatas.append(adata)

adata = anndata.concat(adatas, axis=0)
adata.obs_names_make_unique()

# remove regions waith < %1 coverage across cells
regioncover = np.asarray(adata.X.sum(0)).flatten()
adata = adata[:, regioncover>=0.01*adata.shape[0]].copy()

adata.write('data/test.h5ad')
adata
adata.obs.batch.value_counts()