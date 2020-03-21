
# Microglia contamination analysis


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
adata=sc.read("combIm.h5ad")
```


```python
res=['0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0','1.2','1.4','1.6','1.8','2.0']
```


```python
from matplotlib.cm import get_cmap
cmap_tab10 = get_cmap('tab10')

import matplotlib as mpl
sc.set_figure_params(dpi_save=500)
cmap = mpl.colors.ListedColormap([cmap_tab10.colors[0],cmap_tab10.colors[1],cmap_tab10.colors[2],cmap_tab10.colors[3],cmap_tab10.colors[4],cmap_tab10.colors[5],cmap_tab10.colors[6],
                                 cmap_tab10.colors[7],cmap_tab10.colors[8]])

```


```python
current_palette = sns.color_palette()
ires='1.4'
adata.obsm['X_tsne']=adata.obsm['X_tsne'+ires]
fig, ax0 = plt.subplots(1,1,figsize=[5,5])
sc.pl.tsne(adata,color='celltype',title='Cell type',ax=ax0,palette=current_palette,save='combIm.png')
ax0.axes.get_xaxis().set_visible(False)
ax0.axes.get_yaxis().set_visible(False)

```

    WARNING: saving figure to file figures\tsnecombIm.png
    


![png](output_10_1.png)

