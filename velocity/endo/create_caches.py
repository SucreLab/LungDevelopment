import scvelo as scv
import anndata
import scachepy
import numpy as np

print("### Starting! Working on Endo")
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
adata = scv.read("data/endo_data.h5ad")


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

adata_w_velo.write("data/endo_data_w_velo.h5ad")

print("### Debug info:")
print(str(np.size(adata.X,0)) + " cells in adata")
print(str(np.size(adata_w_velo.X,0)) + " cells in adata_w_velo")


print("## Starting scvelo")
scv.pp.moments(adata_w_velo, n_pcs=30, n_neighbors=30)


c = scachepy.Cache('./cache/cellrank')
c.tl.recover_dynamics(adata_w_velo, force=False) # scv


c.tl.velocity(adata_w_velo, mode='dynamical') # scv
c.tl.velocity_graph(adata_w_velo) # scv