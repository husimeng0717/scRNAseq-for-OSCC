#Integrating Loom File and Meta-data
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
#We now can import our loom file(s) and all of our Seurat meta-data using anndata
OSCC_all = anndata.read_loom("/Bailab7/PROJECT/husimeng/OSCC/velocyto/merged.loom/OSCC_all.loom")
OSCC_all.var_names_make_unique()
sample_obs = pd.read_csv("/Bailab7/PROJECT/husimeng/OSCC/FastMNN/3_5/T/CD4T.velocyto/cellID_obs.csv")
umap_cord = pd.read_csv("/Bailab7/PROJECT/husimeng/OSCC/FastMNN/3_5/T/CD4T.velocyto/umap_embeddings.csv")
cell_Celltype = pd.read_csv("/Bailab7/PROJECT/husimeng/OSCC/FastMNN/3_5/T/CD4T.velocyto/Celltype.csv")
#With our extracted Cell IDs from Seurat, we'll need to filter our uploaded loom (now as an anndata object) based upon them.
CD4T = OSCC_all[np.isin(OSCC_all.obs.index,sample_obs['x'])]
#add UMAP coordinates
CD4T.obs.index
CD4T_index = pd.DataFrame(CD4T.obs.index)
CD4T_index = CD4T_index.rename(columns = {0:'Cell ID'})
umap_cord = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = CD4T_index.merge(umap_cord,on="Cell ID")
umap_ordered
umap_ordered = umap_ordered.iloc[:,1:]
CD4T.obsm['X_umap'] = umap_ordered.values

cell_Celltype = cell_Celltype.rename(columns = {'barcode':'Cell ID'})
cell_Celltype_ordered = CD4T_index.merge(cell_Celltype,on="Cell ID")
cell_Celltype_ordered
cell_Celltype_ordered = cell_Celltype_ordered.iloc[:,2:]
CD4T.obs['Celltype'] = cell_Celltype_ordered.values
#Running RNA Velocity
scv.pp.filter_and_normalize(CD4T)
scv.pp.moments(CD4T)
##Estimate RNA velocity
scv.tl.velocity(CD4T, mode = "stochastic")
scv.tl.velocity_graph(CD4T)
scv.pl.velocity_embedding_stream(CD4T, color='Celltype',  save = 'scv_plot_stream.pdf')
scv.pl.velocity_embedding_stream(CD4T, color='Celltype',  save = 'scv_plot_stream.svg')
scv.pl.velocity_embedding_stream(CD4T, color='Celltype',  save = 'scv_plot_stream.eps')
CD4T.write_loom("CD4T.loom", write_obsm_varm=True)
