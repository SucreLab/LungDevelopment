import scvelo as scv
import anndata
import scachepy
import numpy as np
import cellrank as cr
import scanpy as sc

print("### Starting! Working on Meso")
print("# Loading E12")
e12_loom = scv.read("../loom_files/E12_JS14.loom", cache=False)
e12_loom.var_names_make_unique()
print("# Loading E15")
e15_loom = scv.read("../loom_files/E15_JS12.loom", cache=False)
e15_loom.var_names_make_unique()
print("# Loading E18")
e18_loom = scv.read("../loom_files/E18_JS8.loom", cache=False)
e18_loom.var_names_make_unique()
print("# Loading P0_1")
p0_1_loom = scv.read("../loom_files/P0_1_JS1.loom", cache=False)
p0_1_loom.var_names_make_unique()
print("# Loading P0_2")
p0_2_loom = scv.read("../loom_files/P0_2_JS6.loom", cache=False)
p0_2_loom.var_names_make_unique()
print("# Loading P3")
p3_loom = scv.read("../loom_files/P3_JS2.loom", cache=False)
p3_loom.var_names_make_unique()
print("# Loading P5_1")
p5_1_loom = scv.read("../loom_files/P5_1_JB1.loom", cache=False)
p5_1_loom.var_names_make_unique()
print("# Loading P5_1")
p5_2_loom = scv.read("../loom_files/P5_2_JB3.loom", cache=False)
p5_2_loom.var_names_make_unique()
print("# Loading P7_1")
p7_1_loom = scv.read("../loom_files/P7_1_JS13.loom", cache=False)
p7_1_loom.var_names_make_unique()
print("# Loading P7_2")
p7_2_loom = scv.read("../loom_files/P7_2_JS3.loom", cache=False)
p7_2_loom.var_names_make_unique()
print("# Loading P7_3")
p7_3_loom = scv.read("../loom_files/P7_3_JS9.loom", cache=False)
p7_3_loom.var_names_make_unique()
print("# Loading P14_1")
p14_1_loom = scv.read("../loom_files/P14_1_JS10.loom", cache=False)
p14_1_loom.var_names_make_unique()
print("# Loading P14_2")
p14_2_loom = scv.read("../loom_files/P14_2_JS4.loom", cache=False)
p14_2_loom.var_names_make_unique()
print("# AnnData")
adata = scv.read("data/meso_data.h5ad")


print("# Merging")
adata_w_velo = scv.utils.merge(adata, anndata.concat([e12_loom,
                                                      e15_loom,
                                                      e18_loom,
                                                      p0_1_loom,
                                                      p0_2_loom,
                                                      p3_loom,
                                                      p5_1_loom,
                                                      p5_2_loom,
                                                      p7_1_loom,
                                                      p7_2_loom,
                                                      p7_3_loom,
                                                      p14_1_loom,
                                                      p14_2_loom]))


adata_w_velo.write("data/meso_data_w_velo.h5ad")

print("### Debug info:")
print(str(np.size(adata.X,0)) + " cells in adata")
print(str(np.size(adata_w_velo.X,0)) + " cells in adata_w_velo")


print("## Starting scvelo")
scv.pp.moments(adata_w_velo, n_pcs=30, n_neighbors=30)


c = scachepy.Cache('./cache/cellrank', compression='lzma')
c.tl.recover_dynamics(adata_w_velo, force=False) # scv


c.tl.velocity(adata_w_velo, mode='dynamical') # scv
c.tl.velocity_graph(adata_w_velo) # scv

n_var_genes = 200
top_genes = adata_w_velo.var['fit_likelihood'].sort_values(ascending=False).index[:n_var_genes]
scv.tl.differential_kinetic_test(adata_w_velo, var_names=top_genes, groupby='cell_subtype')

c.tl.velocity(adata_w_velo, diff_kinetics=True, mode='dynamical') # scv
c.tl.velocity_graph(adata_w_velo, fname='velo_graph_diff_kinetics') # scv

adata_w_velo.uns['cell_subtype_colors'] = ["#2ca02c", # Dark green
                                           "#d62728", # Dark red
                                           "#ffbb78", # Light orange
                                           "#9467bd", # Light purple
                                           "#98df8a", # Light green
                                           "#ff7f0e", # Dark orange
                                           "#1f77b4", # Dark blue
                                           "#ff9896", # Light red
                                           "#aec7e8", # Light blue
                                           ]

from cellrank.tl.kernels import VelocityKernel
from cellrank.tl.kernels import ConnectivityKernel

vk = VelocityKernel(adata_w_velo).compute_transition_matrix(seed = 42)
ck = ConnectivityKernel(adata_w_velo).compute_transition_matrix()

combined_kernel = 0.7 * vk + 0.3 * ck
print(combined_kernel)

from cellrank.tl.estimators import GPCCA
g = GPCCA(combined_kernel)
print(g)

g.compute_schur(n_components=100)

g.compute_macrostates(n_states=10, cluster_key='cell_subtype', n_cells = 120)
g.set_terminal_states_from_macrostates()
g.compute_absorption_probabilities()
cr.tl.lineages(adata_w_velo)
cr.pl.lineages(adata_w_velo, same_plot=False)


vk_b = VelocityKernel(adata_w_velo, backward = True).compute_transition_matrix(seed = 42)
ck_b = ConnectivityKernel(adata_w_velo, backward = True).compute_transition_matrix()
combined_kernel_b = 0.2 * vk_b + 0.8 * ck_b
g_b = GPCCA(combined_kernel_b)
g_b.compute_schur(n_components=90)
g_b.plot_spectrum()
g_b.compute_macrostates(n_states=9, cluster_key='cell_subtype', n_cells = 120)
macrostates_list = list(g_b.macrostates.cat.categories)
g_b.set_terminal_states_from_macrostates(names = macrostates_list)
g_b.compute_absorption_probabilities()



adata_lite = anndata.AnnData(X=None,
                             obs = adata_w_velo.obs,
                             var = adata_w_velo.var,
                             uns = adata_w_velo.uns,
                             obsm = adata_w_velo.obsm,
                             varm = adata_w_velo.varm)


keys = ['to_terminal_states_names', 'from_initial_states_names', 'initial_states_names', 'terminal_states_names']
out_dict = {}
for key in keys:
    out_dict[key] = list(adata_lite.uns[key])

print(out_dict)
adata_lite.uns = out_dict

adata_lite.write("./only_metadata.h5ad", compression = 'gzip')

print("Recovering latent time..")
c.tl.recover_latent_time(
    adata_w_velo, root_key="initial_states_probs", end_key="terminal_states_probs"
)

root_idx = np.where(adata_w_velo.obs['initial_states'] == 'Prenatal Wnt2+')[0][0]
adata_w_velo.uns['iroot'] = root_idx
sc.tl.dpt(adata_w_velo)

scv.tl.paga(
    adata_w_velo,
    groups="cell_subtype",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime",
)

scv.tl.score_genes_cell_cycle(adata_w_velo)

scv.tl.velocity_confidence(adata_w_velo)

scv.tl.velocity_pseudotime(adata_w_velo)

cr.tl.lineage_drivers(adata_w_velo)

adata_w_velo.write("./after_cellrank.h5ad", compression = 'gzip') # gzip brings this from 20gb to 3gb

import scipy.io
from scipy.sparse import csr_matrix
scipy.io.mmwrite('meso_expression_Ms.mtx', csr_matrix(adata_w_velo.layers['Ms']))