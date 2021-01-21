from anndata import AnnData
import numpy as np
from scipy.io import mmread
import pandas as pd


prefix = "/Users/negretn/postdoc/code/devo_scseq/velocity/epi/data/epi_"
adata = AnnData(mmread(prefix + "exp_data.mtx.gz").transpose().tocsc())
umap = np.recfromcsv(prefix + "umap.csv")

umap_list = []
for row in umap:
    umap_list.append([row[1], row[2]])


idents = pd.read_csv(prefix + "idents.csv", index_col = 0)
timepoint = pd.read_csv(prefix + "timepoint.csv", index_col = 0)

bulk_celltype = pd.read_csv(prefix + "bulk_celltype.csv", index_col = 0)
cell_subtype = pd.read_csv(prefix + "cell_subtype.csv", index_col = 0)

obs_data = pd.concat([timepoint, bulk_celltype, cell_subtype], axis=1)
obs_data.columns = ['timepoint', 'bulk_cellype', 'cell_subtype']
obs_data = obs_data.astype('category')

genes = pd.read_csv(prefix + "genes.csv").iloc[:,1].values
adata.var = adata.var.set_index(genes)


adata.obsm['X_umap'] = np.array(umap_list)
adata.obsm['X_pca'] = np.array(umap_list)
adata.obs = obs_data

adata.write(prefix + "meso_data.h5ad", compression = 'gzip')
