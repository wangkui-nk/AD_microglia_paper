
# PAGA analysis


```python
import numpy as np
import pandas as pd
%matplotlib inline
import scanpy as sc
from sklearn.metrics import confusion_matrix
from IPython.display import display
import matplotlib.pyplot as plt
import seaborn as sns

```


```python
sc.settings.set_figure_params(dpi=600, frameon=False) 
```


```python
adata=sc.read("microglia.h5ad")
```


```python
figure,axes=plt.subplots(1,1,figsize=[10,5])
ax=sns.violinplot(x=adata.obs['mictype'],y=adata.obs.n_counts,ax=axes,inner=None)
ax.set_xlabel("Subpopulation")
ax.set_ylabel('Total counts')
xlabels = ax.get_xticklabels()
ax.set_xticklabels(xlabels,rotation=45,ha='right',size=10)
```




    [Text(0, 0, 'Homeostatic'),
     Text(0, 0, 'Motile'),
     Text(0, 0, 'ARM'),
     Text(0, 0, 'Dystrophic'),
     Text(0, 0, 'Oli_microglia'),
     Text(0, 0, 'Ast_microglia'),
     Text(0, 0, 'Ex_microglia')]




![png](output_4_1.png)



```python
figure,axes=plt.subplots(1,1,figsize=[10,5])
ax=sns.violinplot(x=adata.obs['mictype'],y=adata.obs.n_genes,ax=axes,inner=None)
ax.set_xlabel("Subpopulation")
ax.set_ylabel('Expressed genes')
xlabels = ax.get_xticklabels()
ax.set_xticklabels(xlabels,rotation=45,ha='right',size=10)
```




    [Text(0, 0, 'Homeostatic'),
     Text(0, 0, 'Motile'),
     Text(0, 0, 'ARM'),
     Text(0, 0, 'Dystrophic'),
     Text(0, 0, 'Oli_microglia'),
     Text(0, 0, 'Ast_microglia'),
     Text(0, 0, 'Ex_microglia')]




![png](output_5_1.png)



```python

figure,axes=plt.subplots(1,1,figsize=[5,5])
ax=sc.pl.scatter(adata, x='n_genes', y='n_counts',color='mictype',title='Subpopulation',ax=axes,show=False)
ax.set_xlabel("Expressed genes")
ax.set_ylabel('Total counts')
xlabels = ax.get_xticklabels()

```


![png](output_6_0.png)



```python
adata
```




    AnnData object with n_obs × n_vars = 3982 × 5000 
        obs: 'cellname', 'sample', 'groupid', 'n_genes', 'percent_mito', 'n_counts', 'n_genes_e', 'n_counts_e', 'desc_0.1', 'desc_0.15', 'desc_0.2', 'desc_0.25', 'desc_0.3', 'desc_0.35', 'desc_0.4', 'desc_0.45', 'desc_0.5', 'desc_0.55', 'desc_0.6', 'desc_0.65', 'desc_0.7', 'desc_0.75', 'desc_0.8', 'trem2', 'atscore', 'apoe', 'orgdesc', 'mictype'
        var: 'Ensembl', 'genename', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
        uns: 'ProjectName', 'prob_matrix0.1', 'prob_matrix0.15', 'prob_matrix0.2', 'prob_matrix0.25', 'prob_matrix0.3', 'prob_matrix0.35', 'prob_matrix0.4', 'prob_matrix0.45', 'prob_matrix0.5', 'prob_matrix0.55', 'prob_matrix0.6', 'prob_matrix0.65', 'prob_matrix0.7', 'prob_matrix0.75', 'prob_matrix0.8', 'mictype_colors'
        obsm: 'X_Embeded_z0.1', 'X_Embeded_z0.15', 'X_Embeded_z0.2', 'X_Embeded_z0.25', 'X_Embeded_z0.3', 'X_Embeded_z0.35', 'X_Embeded_z0.4', 'X_Embeded_z0.45', 'X_Embeded_z0.5', 'X_Embeded_z0.55', 'X_Embeded_z0.6', 'X_Embeded_z0.65', 'X_Embeded_z0.7', 'X_Embeded_z0.75', 'X_Embeded_z0.8', 'X_tsne', 'X_tsne0.1', 'X_tsne0.15', 'X_tsne0.2', 'X_tsne0.25', 'X_tsne0.3', 'X_tsne0.35', 'X_tsne0.4', 'X_tsne0.45', 'X_tsne0.5', 'X_tsne0.55', 'X_tsne0.6', 'X_tsne0.65', 'X_tsne0.7', 'X_tsne0.75', 'X_tsne0.8'




```python
submic=sc.read("submic.h5ad")
```


```python
data=submic.copy()
ires='0.8'
desc='mictype'
root='Homeostatic'
neighbor=20
rand=np.random.randint(low=0,high=229,size=1)[0]

sc.pp.neighbors(data,n_neighbors=neighbor,use_rep='X_Embeded_z'+ires)
# sc.tl.draw_graph(data)
# sc.pl.draw_graph(data, color=desc, legend_loc='on data')
sc.tl.diffmap(data)#,n_comps=10)
sc.pp.neighbors(data,n_neighbors=neighbor,use_rep='X_diffmap')
# sc.tl.draw_graph(data)
# sc.pl.draw_graph(data, color=desc, legend_loc='on data')

sc.tl.paga(data, groups=desc)
sc.pl.paga(data,color=[desc])
# sc.pl.draw_graph(data, color=desc, edges=True, size=20, legend_loc='on data')

data.uns['iroot'] = np.flatnonzero((data.obs[desc] == root))[rand]
sc.tl.dpt(data,n_dcs=15)
sc.pl.diffmap(data, color=['dpt_pseudotime', desc])
sc.pl.tsne(data, color=['dpt_pseudotime',desc])
sc.tl.draw_graph(data)
sc.pl.draw_graph(data, color='dpt_pseudotime')
sc.pl.draw_graph(data, color='dpt_pseudotime', edges=True, size=20,title='DPT_pseudotime')
sc.pl.draw_graph(data, color='trem2', edges=True, size=10,title='Trem2')
sc.pl.draw_graph(data, color='atscore', edges=True, size=10,title='AT Score')
sc.pl.draw_graph(data, color='apoe', edges=True, size=10,title='APOE')
```


![png](output_9_0.png)



![png](output_9_1.png)



![png](output_9_2.png)



![png](output_9_3.png)



![png](output_9_4.png)



![png](output_9_5.png)



![png](output_9_6.png)



![png](output_9_7.png)



```python
data=submic.copy()
ires='0.8'
desc='mictype'
root='Homeostatic'
neighbor=15 #default value
rand=np.random.randint(low=0,high=229,size=1)[0]
sc.tl.pca(data, svd_solver='arpack')
sc.pp.neighbors(data)#,n_neighbors=neighbor)
#sc.tl.draw_graph(data)
#sc.pl.draw_graph(data, color=desc, legend_loc='on data')
sc.tl.diffmap(data)#,n_comps=10)
sc.pp.neighbors(data,n_neighbors=neighbor,use_rep='X_diffmap')
sc.tl.draw_graph(data)
sc.pl.draw_graph(data, color=desc, legend_loc='on data')

sc.tl.paga(data, groups=desc)
sc.pl.paga(data,color=[desc])
sc.pl.paga_compare(data,title='',save='paga.png')
#sc.pl.draw_graph(data, color=desc, edges=True, size=20, legend_loc='on data')

data.uns['iroot'] = np.flatnonzero((data.obs[desc] == root))[rand]
sc.tl.dpt(data,n_dcs=10)
sc.pl.diffmap(data, color=['dpt_pseudotime', desc],title=['DPT_pseudotime','Subpopulation'],save="diffmap.png")
#sc.tl.draw_graph(data, init_pos='paga', maxiter=50)
# sc.pl.draw_graph(data, color='dpt_pseudotime')
# sc.pl.draw_graph(data, color='dpt_pseudotime', edges=True, size=20)
sc.pl.tsne(data, color=['dpt_pseudotime',desc],title=['DPT_pseudotime','Subpopulation'],save="tsne_pseudotime.png")
sc.tl.draw_graph(data)
sc.pl.draw_graph(data, color=desc, legend_loc='on data')
sc.pl.draw_graph(data, color='dpt_pseudotime')
sc.pl.draw_graph(data, color=['dpt_pseudotime','trem2'], edges=True, size=20,save='time_trem2.png',title=['DPT_pseudotime','Trem2'])
#sc.pl.draw_graph(data, color='trem2', edges=True, size=10)
sc.pl.draw_graph(data, color=['atscore','apoe'], edges=True, size=20,save='at_apoe.png',title=['AT Score','APOE'])
#sc.pl.draw_graph(data, color='apoe', edges=True, size=10)

```


![png](output_10_0.png)



![png](output_10_1.png)


    WARNING: saving figure to file figures\paga_comparepaga.png
    


![png](output_10_3.png)


    WARNING: saving figure to file figures\diffmapdiffmap.png
    


![png](output_10_5.png)


    WARNING: saving figure to file figures\tsnetsne_pseudotime.png
    


![png](output_10_7.png)



![png](output_10_8.png)



![png](output_10_9.png)


    WARNING: saving figure to file figures\draw_graph_fatime_trem2.png
    


![png](output_10_11.png)


    WARNING: saving figure to file figures\draw_graph_faat_apoe.png
    


![png](output_10_13.png)



```python
sc.pl.paga_compare(data,title='Subpopulations')
```


![png](output_11_0.png)

