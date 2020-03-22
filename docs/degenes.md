
# Find DE genes



```python
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import confusion_matrix
from IPython.display import display
import matplotlib.pyplot as plt
import os
from adjustText import adjust_text
import seaborn as sns
```


```python
from matplotlib.cm import get_cmap
cmap_tab10 = get_cmap('tab10')

import matplotlib as mpl
```


```python
sc.set_figure_params(dpi_save=300)
cmap = mpl.colors.ListedColormap([cmap_tab10.colors[0],cmap_tab10.colors[3],cmap_tab10.colors[2],cmap_tab10.colors[4],cmap_tab10.colors[1],cmap_tab10.colors[8],cmap_tab10.colors[9]])
```


```python
microglia=sc.read("microglia.h5ad")
```


```python
ires='0.8'
microglia.obsm.X_tsne=microglia.obsm['X_tsne'+ires]
ncols=4
nrows=1
fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[ncols*5,nrows*5])
ax=sc.pl.scatter(microglia,basis='tsne',color='mictype',title='Cell type',ax=axes[0],show=False,legend_fontsize=10)
axes[0].get_legend().remove()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
ax=sc.pl.scatter(microglia,basis='tsne',color='trem2',title='TREM2',legend_loc='upper right',ax=axes[1],show=False)

ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
ax=sc.pl.scatter(microglia,basis='tsne',color='atscore',title='AT score',legend_loc='upper right',ax=axes[2],show=False)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
ax=sc.pl.scatter(microglia,basis='tsne',color='apoe',title='APOE',legend_loc='upper right',ax=axes[3],show=False)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)


```


![png](output_5_0.png)



```python
microglia.obs['mictype'].value_counts()
```




    Motile           1147
    ARM              1119
    Oli_microglia     690
    Ex_microglia      392
    Dystrophic        278
    Homeostatic       229
    Ast_microglia     127
    Name: mictype, dtype: int64




```python
submic=microglia[(microglia.obs['mictype']=='ARM')|(microglia.obs['mictype']=='Motile')|
                 (microglia.obs['mictype']=='Dystrophic')|(microglia.obs['mictype']=="Homeostatic")]
                 
```

DE analysis by wilcoxon rank sum test applied in Scanpy. Two cases: one vs rest and other three vs HOMO

CASE 1: one vs rest


```python
# These script should run in LINUX not windows, sc.tl.rank_genes_groups perform difference in two platform
sc.tl.rank_genes_groups(submic, "mictype",
                       groups=['Homeostatic','ARM','Motile','Dystrophic'],
                        method='wilcoxon',
                        n_genes=submic.raw.var_names.size,
                        rankby_abs=True)
```

    C:\Users\wangk\AppData\Roaming\Python\Python37\site-packages\scanpy\tools\_rank_genes_groups.py:370: RuntimeWarning: overflow encountered in long_scalars
      (ns[imask] * (n_cells - ns[imask]) * (n_cells + 1) / 12))
    


```python
def rank_genes_groups_df(adata, group, pval_cutoff : float =None, logfc_cutoff=None): 
    d = pd.DataFrame() 
    for k in ['scores', 'names', 'logfoldchanges', 'pvals', 'pvals_adj']: 
        d[k] = adata.uns["rank_genes_groups"][k][group] 
    if pval_cutoff is not None: 
        d = d[d["pvals_adj"] < pval_cutoff] 
    if logfc_cutoff is not None: 
        d = d[d["logfoldchanges"].abs() > logfc_cutoff] 
    return d
```


```python
log_fc_threshold=0.1
pval_adj_threshold=0.05
```


```python

de=[]
de_homo = rank_genes_groups_df(submic, 'Homeostatic')
de_arm = rank_genes_groups_df(submic, 'ARM')
de_motile = rank_genes_groups_df(submic, 'Motile')
de_dys = rank_genes_groups_df(submic, 'Dystrophic')

de.append(de_homo)
de.append(de_arm)
de.append(de_motile)
de.append(de_dys)
savefile=["w_homo",'w_arm','w_mot','w_dys']


N=4

for i in range(N):

    de_pos=de[i][(de[i].logfoldchanges>log_fc_threshold)&(de[i].pvals_adj<pval_adj_threshold)]
    de_pos.sort_values('pvals_adj',inplace=True)
    
    de_neg=de[i][(de[i].logfoldchanges<-log_fc_threshold)&(de[i].pvals_adj<pval_adj_threshold)]
    de_neg.sort_values('pvals_adj',inplace=True)

    de_pos.to_csv(savefile[i]+"_pos.csv")
    de_neg.to_csv(savefile[i]+"_neg.csv")


```


```python

de_pos_homo=de_homo[(de_homo.logfoldchanges>log_fc_threshold)&(de_homo.pvals_adj<pval_adj_threshold)]
de_pos_homo.sort_values('pvals_adj',inplace=True)
de_pos_homo.head()
```

    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:2: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>scores</th>
      <th>names</th>
      <th>logfoldchanges</th>
      <th>pvals</th>
      <th>pvals_adj</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>3</th>
      <td>9.575580</td>
      <td>FRMD4A</td>
      <td>1.244736</td>
      <td>1.012879e-21</td>
      <td>6.690823e-18</td>
    </tr>
    <tr>
      <th>7</th>
      <td>8.031395</td>
      <td>SFMBT2</td>
      <td>1.433935</td>
      <td>9.637045e-16</td>
      <td>3.182995e-12</td>
    </tr>
    <tr>
      <th>8</th>
      <td>7.499418</td>
      <td>AC120193.1</td>
      <td>1.504748</td>
      <td>6.410157e-14</td>
      <td>1.881951e-10</td>
    </tr>
    <tr>
      <th>9</th>
      <td>7.431903</td>
      <td>AC068992.1</td>
      <td>4.262547</td>
      <td>1.070459e-13</td>
      <td>2.828473e-10</td>
    </tr>
    <tr>
      <th>10</th>
      <td>7.252969</td>
      <td>AC008691.1</td>
      <td>3.208752</td>
      <td>4.077321e-13</td>
      <td>9.794096e-10</td>
    </tr>
  </tbody>
</table>
</div>




```python
de_neg_homo=de_homo[(de_homo.logfoldchanges<-log_fc_threshold)&(de_homo.pvals_adj<pval_adj_threshold)]
de_neg_homo.sort_values('pvals_adj',inplace=True)
de_neg_homo.head()
```

    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:2: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>scores</th>
      <th>names</th>
      <th>logfoldchanges</th>
      <th>pvals</th>
      <th>pvals_adj</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>-12.647748</td>
      <td>NEAT1</td>
      <td>-2.312541</td>
      <td>1.151283e-36</td>
      <td>3.042035e-32</td>
    </tr>
    <tr>
      <th>1</th>
      <td>-11.580046</td>
      <td>SLC26A3</td>
      <td>-2.328851</td>
      <td>5.201821e-31</td>
      <td>6.872386e-27</td>
    </tr>
    <tr>
      <th>2</th>
      <td>-10.786842</td>
      <td>FKBP5</td>
      <td>-3.545546</td>
      <td>3.972054e-27</td>
      <td>3.498452e-23</td>
    </tr>
    <tr>
      <th>4</th>
      <td>-9.420601</td>
      <td>DPYD</td>
      <td>-3.074943</td>
      <td>4.485143e-21</td>
      <td>2.370219e-17</td>
    </tr>
    <tr>
      <th>5</th>
      <td>-8.901549</td>
      <td>SLC1A3</td>
      <td>-1.497096</td>
      <td>5.507224e-19</td>
      <td>2.425290e-15</td>
    </tr>
  </tbody>
</table>
</div>



CASE 2: Arm, Motile, Dysgtrophic vs HOMOE


```python
sc.tl.rank_genes_groups(submic, "mictype",
                        groups=['ARM','Motile','Dystrophic'],
                        reference="Homeostatic",
                        method='wilcoxon',
                        n_genes=submic.raw.var_names.size,
                        rankby_abs=True)
```


```python
de_arm_h = rank_genes_groups_df(submic, 'ARM')
de_mot_h = rank_genes_groups_df(submic, 'Motile')
de_dys_h = rank_genes_groups_df(submic, 'Dystrophic')
```


```python
de=[]

de.append(de_arm_h)
de.append(de_mot_h)
de.append(de_dys_h)
savefile=['w_arm_h','w_mot_h','w_dys_h']


N=3

for i in range(N):

    de_pos=de[i][(de[i].logfoldchanges>log_fc_threshold)&(de[i].pvals_adj<pval_adj_threshold)]
    de_pos.sort_values('pvals_adj',inplace=True)
    
    de_neg=de[i][(de[i].logfoldchanges<-log_fc_threshold)&(de[i].pvals_adj<pval_adj_threshold)]
    de_neg.sort_values('pvals_adj',inplace=True)

    de_pos.to_csv(savefile[i]+"_pos.csv")
    de_neg.to_csv(savefile[i]+"_neg.csv")
```

    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:14: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      
    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:17: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:14: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      
    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:17: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:14: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      
    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:17: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
    


```python
i=0 # arm_vs_homo
de_pos=de[i][(de[i].logfoldchanges>log_fc_threshold)&(de[i].pvals_adj<pval_adj_threshold)]
de_pos.sort_values('pvals_adj',inplace=True)
de_pos.head()
```

    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:3: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      This is separate from the ipykernel package so we can avoid doing imports until
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>scores</th>
      <th>names</th>
      <th>logfoldchanges</th>
      <th>pvals</th>
      <th>pvals_adj</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>12.081801</td>
      <td>SLC26A3</td>
      <td>2.513461</td>
      <td>1.317995e-33</td>
      <td>1.741269e-29</td>
    </tr>
    <tr>
      <th>2</th>
      <td>11.397553</td>
      <td>NEAT1</td>
      <td>2.180482</td>
      <td>4.300338e-30</td>
      <td>3.787594e-26</td>
    </tr>
    <tr>
      <th>3</th>
      <td>10.581952</td>
      <td>DPYD</td>
      <td>3.488856</td>
      <td>3.613500e-26</td>
      <td>2.386988e-22</td>
    </tr>
    <tr>
      <th>6</th>
      <td>9.791317</td>
      <td>FKBP5</td>
      <td>3.413069</td>
      <td>1.226873e-22</td>
      <td>4.631093e-19</td>
    </tr>
    <tr>
      <th>7</th>
      <td>9.542397</td>
      <td>SLC1A3</td>
      <td>1.691859</td>
      <td>1.395673e-21</td>
      <td>4.609732e-18</td>
    </tr>
  </tbody>
</table>
</div>




```python
i=1 # motile_vs_homo
de_pos=de[i][(de[i].logfoldchanges>log_fc_threshold)&(de[i].pvals_adj<pval_adj_threshold)]
de_pos.sort_values('pvals_adj',inplace=True)
de_pos.head()
```

    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:3: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      This is separate from the ipykernel package so we can avoid doing imports until
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>scores</th>
      <th>names</th>
      <th>logfoldchanges</th>
      <th>pvals</th>
      <th>pvals_adj</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>13.623844</td>
      <td>NEAT1</td>
      <td>2.520373</td>
      <td>2.889481e-42</td>
      <td>7.634876e-38</td>
    </tr>
    <tr>
      <th>1</th>
      <td>10.718767</td>
      <td>FKBP5</td>
      <td>3.671437</td>
      <td>8.310355e-27</td>
      <td>1.097923e-22</td>
    </tr>
    <tr>
      <th>2</th>
      <td>8.828892</td>
      <td>SLC26A3</td>
      <td>1.938071</td>
      <td>1.057182e-18</td>
      <td>9.311303e-15</td>
    </tr>
    <tr>
      <th>3</th>
      <td>7.947751</td>
      <td>SLC1A3</td>
      <td>1.368051</td>
      <td>1.899279e-15</td>
      <td>1.254616e-11</td>
    </tr>
    <tr>
      <th>5</th>
      <td>7.667789</td>
      <td>APOE</td>
      <td>1.923241</td>
      <td>1.749855e-14</td>
      <td>7.706070e-11</td>
    </tr>
  </tbody>
</table>
</div>




```python
i=2 # dystrohpic_vs_homo
de_pos=de[i][(de[i].logfoldchanges>log_fc_threshold)&(de[i].pvals_adj<pval_adj_threshold)]
de_pos.sort_values('pvals_adj',inplace=True)
de_pos.head()
```

    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:3: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      This is separate from the ipykernel package so we can avoid doing imports until
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>scores</th>
      <th>names</th>
      <th>logfoldchanges</th>
      <th>pvals</th>
      <th>pvals_adj</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>15.575459</td>
      <td>FTL</td>
      <td>3.922354</td>
      <td>1.068825e-54</td>
      <td>2.824156e-50</td>
    </tr>
    <tr>
      <th>1</th>
      <td>13.127926</td>
      <td>SLC26A3</td>
      <td>3.167222</td>
      <td>2.278097e-39</td>
      <td>3.009707e-35</td>
    </tr>
    <tr>
      <th>2</th>
      <td>11.868525</td>
      <td>RPS6</td>
      <td>4.330101</td>
      <td>1.724818e-32</td>
      <td>1.519162e-28</td>
    </tr>
    <tr>
      <th>3</th>
      <td>11.729336</td>
      <td>FTH1</td>
      <td>2.991559</td>
      <td>9.016316e-32</td>
      <td>5.955953e-28</td>
    </tr>
    <tr>
      <th>4</th>
      <td>10.799176</td>
      <td>RPL19</td>
      <td>4.029462</td>
      <td>3.473063e-27</td>
      <td>1.835375e-23</td>
    </tr>
  </tbody>
</table>
</div>


