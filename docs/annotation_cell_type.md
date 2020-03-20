
# Cell type annotation


```python
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import confusion_matrix
from IPython.display import display
import matplotlib.pyplot as plt
import seaborn as sns
```


```python
mic=["CD40","CD68","CEBPA","CEBPB","CX3CR1","ITGAM","PTPRC","TYROBP","APOE"]
ast=["GFAP","ALDH1L1","AQP4","SLC1A3","SLC1A2","GLUL"]
glu=["SLC17A7","GRIN1","GRIN2B"]
gab=["GAD1","GAD2","CCK"]
oli=["SOX10","MBP","MOG","OLIG1","OLIG2","PDGFRA","VCAN"]
end=['FLT1']
mks=mic+ast+glu+gab+oli+end
```


```python
res=['0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8']
```


```python
adata=sc.read("brain.h5ad")
```


```python
sc.set_figure_params(dpi_save=300)
```


```python
fig = plt.figure(figsize=(24,12))

ax1 = plt.subplot2grid((2, 8), (0, 0),colspan=1)
ax2 = plt.subplot2grid((2, 8), (0, 1), colspan=7)
ax3 = plt.subplot2grid((2, 8), (1, 0),colspan=1)
ax4 = plt.subplot2grid((2, 8), (1, 1), colspan=7)

ax1 = sns.violinplot(y='n_counts',data=adata.obs,ax=ax1)
ax1.set_xticklabels(labels=['ALL'])
ax1.set_ylabel("Counts")
ax2 = sns.violinplot('sampleID', y='n_counts', data=adata.obs,ax=ax2,palette=sns.color_palette("muted"))
ax2.set_ylabel("Counts")
ax3 = sns.violinplot(y='n_genes',data=adata.obs,ax=ax3)
ax3.set_xticklabels(labels=['ALL'])
ax3.set_ylabel("Counts")
ax4 = sns.violinplot('sampleID', y='n_genes', data=adata.obs,ax=ax4,palette=sns.color_palette("muted"))
ax4.set_ylabel("Counts")
plt.tight_layout()

```


![png](output_6_0.png)



```python
fig,axes=plt.subplots(nrows=4,ncols=4,figsize=(24,24))
axes[0][0]=adata.obs.plot.scatter(x='n_genes',y='n_counts',title='ALL',ax=axes[0][0])
axes[0][0].set_xlim(0,4100)
axes[0][0].set_ylim(0,12500)
axes[0][0].set_xlabel("Genes")
axes[0][0].set_ylabel("Counts")
i=0
j=1
str_cat=adata.obs['sampleID'].cat.categories
for i in range(0,4):
    for j in range(0,4):
        if i==0 and j==0:
            continue
        data=adata[adata.obs['sampleID']==str_cat[i*4+j-1]]
        axes[i][j]=data.obs.plot.scatter(x='n_genes',y='n_counts',title=str_cat[i*4+j-1],ax=axes[i][j])
        axes[i][j].set_xlim(0,4100)
        axes[i][j].set_ylim(0,12500)
        axes[i][j].set_xlabel("Genes")
        axes[i][j].set_ylabel("Counts")
plt.tight_layout()

```


![png](output_7_0.png)



```python
res=['0.2','0.4','0.6','0.8']
i=0
for ires in res:
    adata.obsm.X_tsne=adata.obsm['X_tsne'+ires]
    adata.obs['maxprob']=(np.max(adata.uns['prob_matrix'+ires],axis=1)*1000).astype(int)
    adata.obsm.X_tsne=adata.obsm['X_tsne'+ires]
    sc.pl.tsne(adata, color=['desc_'+ires,'n_counts','maxprob','sampleID'],title=['desc_'+ires,'Counts','maxprob','sampleID'])#,ax=axes[i])
    i=i+1

```


![png](output_8_0.png)



![png](output_8_1.png)



![png](output_8_2.png)



![png](output_8_3.png)



```python
ires='0.1'
adata.obsm.X_tsne=adata.obsm['X_tsne'+ires]
sc.pl.tsne(adata, color=['desc_'+ires],palette='tab20')
for cell in ['ast','mic','oli','glu','gab','end']:
    print(cell)
    sc.pl.tsne(adata, color=['desc_'+ires]+vars()[cell],color_map="YlGn",size=20)
```


![png](output_9_0.png)


    ast
    


![png](output_9_2.png)


    mic
    


![png](output_9_4.png)


    oli
    


![png](output_9_6.png)


    glu
    


![png](output_9_8.png)


    gab
    


![png](output_9_10.png)


    end
    


![png](output_9_12.png)



```python
for i in range(0,len(res)-1):
    print("desc_"+res[i]+" vs desc_"+res[i+1])
    tab=confusion_matrix(adata.obs['desc_'+res[i]],adata.obs['desc_'+res[i+1]])
    tab=tab[np.sum(tab,axis=1)!=0,]
    display(pd.DataFrame(tab))
```

    desc_0.2 vs desc_0.4
    


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
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
      <th>8</th>
      <th>9</th>
      <th>10</th>
      <th>11</th>
      <th>12</th>
      <th>13</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>44303</td>
      <td>356</td>
      <td>8</td>
      <td>68</td>
      <td>2</td>
      <td>28</td>
      <td>15</td>
      <td>139</td>
      <td>9</td>
      <td>26</td>
      <td>39</td>
      <td>455</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3</td>
      <td>52</td>
      <td>1</td>
      <td>12751</td>
      <td>1</td>
      <td>3</td>
      <td>1</td>
      <td>4071</td>
      <td>10</td>
      <td>0</td>
      <td>2</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>52</td>
      <td>14</td>
      <td>15823</td>
      <td>45</td>
      <td>4</td>
      <td>9</td>
      <td>3</td>
      <td>47</td>
      <td>4</td>
      <td>3</td>
      <td>3</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>20</td>
      <td>15793</td>
      <td>5</td>
      <td>612</td>
      <td>3</td>
      <td>4</td>
      <td>6</td>
      <td>1921</td>
      <td>5</td>
      <td>0</td>
      <td>4</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>38</td>
      <td>13</td>
      <td>0</td>
      <td>9</td>
      <td>6905</td>
      <td>6</td>
      <td>1</td>
      <td>15</td>
      <td>1</td>
      <td>0</td>
      <td>3</td>
      <td>19</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>5</th>
      <td>8</td>
      <td>5</td>
      <td>2</td>
      <td>6</td>
      <td>1</td>
      <td>6417</td>
      <td>10</td>
      <td>48</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>12</td>
      <td>9</td>
      <td>0</td>
      <td>7</td>
      <td>2</td>
      <td>11</td>
      <td>6189</td>
      <td>16</td>
      <td>1</td>
      <td>0</td>
      <td>11</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>4</td>
      <td>94</td>
      <td>4</td>
      <td>30</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>666</td>
      <td>3632</td>
      <td>1</td>
      <td>1051</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>26</td>
      <td>7</td>
      <td>2</td>
      <td>10</td>
      <td>1</td>
      <td>1</td>
      <td>7</td>
      <td>22</td>
      <td>2</td>
      <td>4081</td>
      <td>1</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>1965</td>
      <td>43</td>
      <td>17</td>
      <td>29</td>
      <td>3</td>
      <td>38</td>
      <td>8</td>
      <td>51</td>
      <td>7</td>
      <td>33</td>
      <td>38</td>
      <td>1720</td>
      <td>4</td>
      <td>1</td>
    </tr>
    <tr>
      <th>10</th>
      <td>3</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>651</td>
      <td>0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>456</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.4 vs desc_0.6
    


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
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
      <th>8</th>
      <th>9</th>
      <th>10</th>
      <th>11</th>
      <th>12</th>
      <th>13</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>46390</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>35</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3</td>
      <td>16366</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>18</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0</td>
      <td>0</td>
      <td>15859</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>13559</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6925</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>6517</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6239</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>0</td>
      <td>11</td>
      <td>1</td>
      <td>10</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>6968</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3671</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4142</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>10</th>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1148</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2217</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>657</td>
      <td>0</td>
    </tr>
    <tr>
      <th>13</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.6 vs desc_0.8
    


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
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
      <th>8</th>
      <th>9</th>
      <th>10</th>
      <th>11</th>
      <th>12</th>
      <th>13</th>
      <th>14</th>
      <th>15</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>46070</td>
      <td>114</td>
      <td>11</td>
      <td>51</td>
      <td>3</td>
      <td>16</td>
      <td>37</td>
      <td>13</td>
      <td>1</td>
      <td>9</td>
      <td>6</td>
      <td>13</td>
      <td>3</td>
      <td>46</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>32</td>
      <td>15160</td>
      <td>3</td>
      <td>11</td>
      <td>4</td>
      <td>4</td>
      <td>568</td>
      <td>5</td>
      <td>6</td>
      <td>0</td>
      <td>4</td>
      <td>584</td>
      <td>0</td>
      <td>2</td>
      <td>2</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>8</td>
      <td>1</td>
      <td>15767</td>
      <td>8</td>
      <td>3</td>
      <td>8</td>
      <td>26</td>
      <td>10</td>
      <td>1</td>
      <td>3</td>
      <td>13</td>
      <td>11</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>2</td>
      <td>0</td>
      <td>9968</td>
      <td>2</td>
      <td>1</td>
      <td>5</td>
      <td>1</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>3584</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>14</td>
      <td>1</td>
      <td>1</td>
      <td>4</td>
      <td>6888</td>
      <td>2</td>
      <td>2</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>6</td>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>7</td>
      <td>13</td>
      <td>1</td>
      <td>3</td>
      <td>0</td>
      <td>11</td>
      <td>18</td>
      <td>4183</td>
      <td>0</td>
      <td>0</td>
      <td>2273</td>
      <td>6</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>6133</td>
      <td>10</td>
      <td>4</td>
      <td>1</td>
      <td>0</td>
      <td>84</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>6</td>
      <td>98</td>
      <td>4</td>
      <td>41</td>
      <td>0</td>
      <td>5</td>
      <td>5719</td>
      <td>9</td>
      <td>29</td>
      <td>4</td>
      <td>5</td>
      <td>1071</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>3</td>
      <td>3</td>
      <td>0</td>
      <td>7</td>
      <td>0</td>
      <td>2</td>
      <td>21</td>
      <td>2</td>
      <td>3620</td>
      <td>0</td>
      <td>2</td>
      <td>8</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>11</td>
      <td>7</td>
      <td>4</td>
      <td>6</td>
      <td>2</td>
      <td>5</td>
      <td>6</td>
      <td>7</td>
      <td>1</td>
      <td>4082</td>
      <td>7</td>
      <td>5</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>10</th>
      <td>18</td>
      <td>2</td>
      <td>1</td>
      <td>4</td>
      <td>2</td>
      <td>12</td>
      <td>17</td>
      <td>4</td>
      <td>13</td>
      <td>0</td>
      <td>2</td>
      <td>4</td>
      <td>1067</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>505</td>
      <td>57</td>
      <td>18</td>
      <td>14</td>
      <td>1</td>
      <td>17</td>
      <td>105</td>
      <td>32</td>
      <td>0</td>
      <td>8</td>
      <td>9</td>
      <td>18</td>
      <td>3</td>
      <td>1477</td>
      <td>3</td>
      <td>0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>654</td>
      <td>0</td>
    </tr>
    <tr>
      <th>13</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>454</td>
    </tr>
  </tbody>
</table>
</div>



```python
for ires in res:
    ax = sc.pl.dotplot(adata, mks, groupby='desc_'+ires,var_group_positions=[(0,8),(9,14),(15,17),(18,20),(21,26),(27,27),(28,28)],var_group_labels=['mic','ast','glu','gab','oli','opc','end'])
```


![png](output_11_0.png)



![png](output_11_1.png)



![png](output_11_2.png)



![png](output_11_3.png)



```python
#annotation manually
res=['0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8']
celltype=[None]*15

celltype[0]=['oli','glu','ast','gab','oli','glu','mic','oli']#0.1
celltype[1]=['oli','glu','ast','opc','gab','gab','glu','mic','oli','end']
celltype[2]=['oli','glu','ast','glu','opc','gab','gab','glu','mic','oli','oli','end']
celltype[3]=['oli','glu','ast','glu','opc','gab','gab','glu','mic','oli','oli','end']#0.25
celltype[4]=['oli','ast','glu','glu','opc','gab','gab','glu','glu','mic','oli','oli','end']
celltype[5]=['oli','ast','glu','glu','opc','gab','gab','glu','glu','mic','glu','oli','oli','end']
celltype[6]=['oli','glu','ast','glu','opc','gab','gab','glu','glu','mic','glu','oli','oli','end']#0.4
celltype[7]=['oli','glu','ast','glu','opc','gab','gab','glu','glu','mic','glu','oli','oli','end']
celltype[8]=['oli','glu','ast','glu','opc','gab','gab','glu','glu','mic','glu','oli','oli','end']
celltype[9]=['oli','glu','ast','glu','opc','gab','gab','glu','glu','mic','glu','oli','oli','end']
celltype[10]=['oli','glu','ast','glu','opc','gab','gab','glu','glu','mic','glu','oli','oli','end']
celltype[11]=['oli','glu','ast','glu','opc','gab','gab','glu','glu','mic','glu','glu','oli','oli','end']#0.65
celltype[12]=['oli','glu','ast','glu','opc','gab','gab','glu','glu','mic','glu','glu','oli','oli','end']
celltype[13]=['oli','glu','ast','glu','opc','gab','gab','glu','glu','mic','glu','glu','oli','oli','end']
celltype[14]=['oli','glu','ast','glu','opc','gab','glu','gab','glu','mic','gab','glu','glu','oli','oli','end']
for i in range(0,len(res)):
   adata.obs['celltype_'+res[i]] =adata.obs['desc_'+res[i]]
   adata.obs['celltype_'+res[i]] .replace(range(0,len(celltype[i])), celltype[i], inplace = True)  
```


```python
for i in range(0,len(res)-1):
   print("desc_"+res[i]+" vs desc_"+res[i+1])
   tab=confusion_matrix(adata.obs['desc_'+res[i]],adata.obs['desc_'+res[i+1]])
   tab=tab[np.sum(tab,axis=1)!=0,]
 #  display(pd.DataFrame(tab))
   display(pd.DataFrame(tab,index=celltype[i],columns=celltype[i+1]))
```

    desc_0.1 vs desc_0.15
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>mic</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>37452</td>
      <td>16</td>
      <td>11</td>
      <td>88</td>
      <td>54</td>
      <td>46</td>
      <td>54</td>
      <td>79</td>
      <td>651</td>
      <td>3</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>12</td>
      <td>23891</td>
      <td>7</td>
      <td>3</td>
      <td>89</td>
      <td>25</td>
      <td>880</td>
      <td>8</td>
      <td>11</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>9</td>
      <td>12</td>
      <td>16219</td>
      <td>14</td>
      <td>25</td>
      <td>42</td>
      <td>4</td>
      <td>2</td>
      <td>20</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>4</td>
      <td>5</td>
      <td>1</td>
      <td>2</td>
      <td>6465</td>
      <td>6280</td>
      <td>85</td>
      <td>1</td>
      <td>5</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>536</td>
      <td>5</td>
      <td>3</td>
      <td>8720</td>
      <td>11</td>
      <td>10</td>
      <td>8</td>
      <td>7</td>
      <td>107</td>
      <td>1</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>11</td>
      <td>177</td>
      <td>12</td>
      <td>2</td>
      <td>108</td>
      <td>126</td>
      <td>14350</td>
      <td>5</td>
      <td>8</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>764</td>
      <td>18</td>
      <td>9</td>
      <td>10</td>
      <td>24</td>
      <td>19</td>
      <td>12</td>
      <td>5754</td>
      <td>113</td>
      <td>2</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>1689</td>
      <td>7</td>
      <td>9</td>
      <td>46</td>
      <td>12</td>
      <td>18</td>
      <td>9</td>
      <td>99</td>
      <td>5385</td>
      <td>458</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.15 vs desc_0.2
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>mic</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>38712</td>
      <td>29</td>
      <td>6</td>
      <td>320</td>
      <td>7</td>
      <td>12</td>
      <td>12</td>
      <td>16</td>
      <td>4</td>
      <td>1359</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>16</td>
      <td>16498</td>
      <td>4</td>
      <td>7501</td>
      <td>7</td>
      <td>7</td>
      <td>3</td>
      <td>65</td>
      <td>12</td>
      <td>16</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>26</td>
      <td>27</td>
      <td>15956</td>
      <td>93</td>
      <td>5</td>
      <td>17</td>
      <td>11</td>
      <td>35</td>
      <td>23</td>
      <td>77</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>1556</td>
      <td>13</td>
      <td>14</td>
      <td>49</td>
      <td>6978</td>
      <td>11</td>
      <td>12</td>
      <td>10</td>
      <td>10</td>
      <td>229</td>
      <td>2</td>
      <td>1</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>16</td>
      <td>27</td>
      <td>0</td>
      <td>276</td>
      <td>1</td>
      <td>6409</td>
      <td>4</td>
      <td>18</td>
      <td>2</td>
      <td>35</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>17</td>
      <td>3</td>
      <td>1</td>
      <td>304</td>
      <td>0</td>
      <td>24</td>
      <td>6187</td>
      <td>14</td>
      <td>4</td>
      <td>10</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>18</td>
      <td>271</td>
      <td>3</td>
      <td>9743</td>
      <td>3</td>
      <td>10</td>
      <td>14</td>
      <td>5313</td>
      <td>10</td>
      <td>17</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>700</td>
      <td>13</td>
      <td>4</td>
      <td>40</td>
      <td>0</td>
      <td>3</td>
      <td>3</td>
      <td>6</td>
      <td>4080</td>
      <td>1105</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>4382</td>
      <td>18</td>
      <td>26</td>
      <td>47</td>
      <td>10</td>
      <td>9</td>
      <td>13</td>
      <td>11</td>
      <td>18</td>
      <td>1109</td>
      <td>651</td>
      <td>6</td>
    </tr>
    <tr>
      <th>end</th>
      <td>7</td>
      <td>1</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>453</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.2 vs desc_0.25
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>mic</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>45286</td>
      <td>1</td>
      <td>1</td>
      <td>9</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>151</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>16887</td>
      <td>0</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>2</td>
      <td>0</td>
      <td>16004</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>32</td>
      <td>1</td>
      <td>18332</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>4</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>7007</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>6497</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6257</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>4</td>
      <td>1</td>
      <td>10</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>5472</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4159</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>24</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>3929</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>657</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>462</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.25 vs desc_0.3
    


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
      <th>oli</th>
      <th>ast</th>
      <th>glu</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>44466</td>
      <td>1</td>
      <td>32</td>
      <td>75</td>
      <td>0</td>
      <td>4</td>
      <td>5</td>
      <td>46</td>
      <td>5</td>
      <td>2</td>
      <td>681</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>4</td>
      <td>2</td>
      <td>12641</td>
      <td>15</td>
      <td>3</td>
      <td>3</td>
      <td>1</td>
      <td>4234</td>
      <td>15</td>
      <td>1</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>20</td>
      <td>15867</td>
      <td>41</td>
      <td>11</td>
      <td>2</td>
      <td>10</td>
      <td>3</td>
      <td>39</td>
      <td>6</td>
      <td>1</td>
      <td>10</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>26</td>
      <td>2</td>
      <td>276</td>
      <td>15451</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>2576</td>
      <td>3</td>
      <td>0</td>
      <td>21</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>25</td>
      <td>1</td>
      <td>3</td>
      <td>3</td>
      <td>6928</td>
      <td>1</td>
      <td>0</td>
      <td>8</td>
      <td>2</td>
      <td>0</td>
      <td>37</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>3</td>
      <td>2</td>
      <td>6</td>
      <td>7</td>
      <td>0</td>
      <td>6415</td>
      <td>4</td>
      <td>49</td>
      <td>6</td>
      <td>0</td>
      <td>5</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>9</td>
      <td>0</td>
      <td>7</td>
      <td>8</td>
      <td>1</td>
      <td>7</td>
      <td>6197</td>
      <td>18</td>
      <td>5</td>
      <td>0</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>0</td>
      <td>40</td>
      <td>20</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>740</td>
      <td>4670</td>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>10</td>
      <td>3</td>
      <td>6</td>
      <td>6</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>11</td>
      <td>8</td>
      <td>4109</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>776</td>
      <td>6</td>
      <td>18</td>
      <td>14</td>
      <td>1</td>
      <td>10</td>
      <td>2</td>
      <td>32</td>
      <td>8</td>
      <td>3</td>
      <td>3228</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>656</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>458</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.3 vs desc_0.35
    


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
      <th>oli</th>
      <th>ast</th>
      <th>glu</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>44932</td>
      <td>4</td>
      <td>65</td>
      <td>28</td>
      <td>2</td>
      <td>0</td>
      <td>3</td>
      <td>15</td>
      <td>7</td>
      <td>5</td>
      <td>5</td>
      <td>273</td>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>11</td>
      <td>15845</td>
      <td>3</td>
      <td>3</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>10</td>
      <td>0</td>
      <td>3</td>
      <td>1</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>4</td>
      <td>2</td>
      <td>4</td>
      <td>12932</td>
      <td>0</td>
      <td>2</td>
      <td>2</td>
      <td>119</td>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>15</td>
      <td>4</td>
      <td>15437</td>
      <td>19</td>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>122</td>
      <td>5</td>
      <td>2</td>
      <td>1</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>10</td>
      <td>0</td>
      <td>5</td>
      <td>9</td>
      <td>6904</td>
      <td>0</td>
      <td>1</td>
      <td>6</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>2</td>
      <td>6425</td>
      <td>7</td>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>8</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>2</td>
      <td>5</td>
      <td>6189</td>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>6</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>6</td>
      <td>2</td>
      <td>452</td>
      <td>252</td>
      <td>0</td>
      <td>8</td>
      <td>6</td>
      <td>6974</td>
      <td>28</td>
      <td>0</td>
      <td>8</td>
      <td>17</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>3</td>
      <td>3</td>
      <td>3</td>
      <td>5</td>
      <td>2</td>
      <td>5</td>
      <td>2</td>
      <td>13</td>
      <td>3633</td>
      <td>6</td>
      <td>1052</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>3</td>
      <td>1</td>
      <td>3</td>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>6</td>
      <td>0</td>
      <td>4095</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>963</td>
      <td>16</td>
      <td>35</td>
      <td>27</td>
      <td>1</td>
      <td>21</td>
      <td>13</td>
      <td>49</td>
      <td>7</td>
      <td>24</td>
      <td>14</td>
      <td>2838</td>
      <td>3</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>657</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>457</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.35 vs desc_0.4
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>45480</td>
      <td>201</td>
      <td>7</td>
      <td>16</td>
      <td>6</td>
      <td>12</td>
      <td>7</td>
      <td>43</td>
      <td>2</td>
      <td>11</td>
      <td>15</td>
      <td>150</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>24</td>
      <td>4</td>
      <td>15835</td>
      <td>3</td>
      <td>2</td>
      <td>2</td>
      <td>1</td>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>15</td>
      <td>15404</td>
      <td>4</td>
      <td>153</td>
      <td>2</td>
      <td>12</td>
      <td>5</td>
      <td>411</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>13</td>
      <td>24</td>
      <td>1</td>
      <td>13171</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>55</td>
      <td>4</td>
      <td>2</td>
      <td>4</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>3</td>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>6903</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>5</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>1</td>
      <td>6439</td>
      <td>6</td>
      <td>8</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>9</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>7</td>
      <td>6206</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>4</td>
      <td>680</td>
      <td>5</td>
      <td>208</td>
      <td>1</td>
      <td>5</td>
      <td>1</td>
      <td>6403</td>
      <td>14</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>21</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>13</td>
      <td>3638</td>
      <td>0</td>
      <td>6</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>15</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>5</td>
      <td>9</td>
      <td>0</td>
      <td>4103</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>6</td>
      <td>0</td>
      <td>1083</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>861</td>
      <td>45</td>
      <td>8</td>
      <td>9</td>
      <td>6</td>
      <td>36</td>
      <td>9</td>
      <td>50</td>
      <td>5</td>
      <td>24</td>
      <td>41</td>
      <td>2057</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>654</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>456</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.4 vs desc_0.45
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>46423</td>
      <td>5</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>16374</td>
      <td>0</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>13</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>15862</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>13555</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6925</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>6517</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6240</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>13</td>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>6970</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3671</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4144</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1152</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>2210</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>657</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.45 vs desc_0.5
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>46422</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>16389</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>15863</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>13564</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6926</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6517</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6240</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6990</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3673</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4146</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1151</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>4</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>2210</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>658</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.5 vs desc_0.55
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>46413</td>
      <td>3</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>5</td>
      <td>16373</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>16</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>15862</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>13551</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>12</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6926</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6516</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>6239</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>15</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6977</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3674</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4146</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>1150</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2205</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>657</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.55 vs desc_0.6
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>46388</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>30</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>3</td>
      <td>16368</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>15859</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>13558</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>10</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6926</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>6516</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6239</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>13</td>
      <td>1</td>
      <td>10</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>6976</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3673</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4143</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1146</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2216</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>657</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.6 vs desc_0.65
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>45842</td>
      <td>72</td>
      <td>5</td>
      <td>52</td>
      <td>2</td>
      <td>3</td>
      <td>3</td>
      <td>32</td>
      <td>0</td>
      <td>6</td>
      <td>11</td>
      <td>1</td>
      <td>366</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>32</td>
      <td>15074</td>
      <td>2</td>
      <td>19</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>591</td>
      <td>5</td>
      <td>0</td>
      <td>630</td>
      <td>0</td>
      <td>26</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>12</td>
      <td>3</td>
      <td>15777</td>
      <td>12</td>
      <td>2</td>
      <td>4</td>
      <td>6</td>
      <td>22</td>
      <td>2</td>
      <td>2</td>
      <td>15</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>4</td>
      <td>2</td>
      <td>0</td>
      <td>10019</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>2</td>
      <td>1</td>
      <td>3534</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>11</td>
      <td>4</td>
      <td>4</td>
      <td>3</td>
      <td>6893</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>3</td>
      <td>8</td>
      <td>0</td>
      <td>10</td>
      <td>0</td>
      <td>6455</td>
      <td>1</td>
      <td>12</td>
      <td>1</td>
      <td>0</td>
      <td>8</td>
      <td>0</td>
      <td>19</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>9</td>
      <td>2</td>
      <td>0</td>
      <td>6</td>
      <td>0</td>
      <td>24</td>
      <td>6174</td>
      <td>14</td>
      <td>3</td>
      <td>1</td>
      <td>5</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>7</td>
      <td>81</td>
      <td>3</td>
      <td>46</td>
      <td>0</td>
      <td>1</td>
      <td>2</td>
      <td>5752</td>
      <td>34</td>
      <td>7</td>
      <td>1052</td>
      <td>0</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>16</td>
      <td>3630</td>
      <td>0</td>
      <td>4</td>
      <td>8</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>8</td>
      <td>5</td>
      <td>2</td>
      <td>6</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>6</td>
      <td>0</td>
      <td>4089</td>
      <td>4</td>
      <td>0</td>
      <td>20</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>10</td>
      <td>2</td>
      <td>0</td>
      <td>6</td>
      <td>1</td>
      <td>0</td>
      <td>7</td>
      <td>13</td>
      <td>20</td>
      <td>0</td>
      <td>4</td>
      <td>1060</td>
      <td>26</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>118</td>
      <td>8</td>
      <td>4</td>
      <td>10</td>
      <td>0</td>
      <td>3</td>
      <td>1</td>
      <td>30</td>
      <td>1</td>
      <td>2</td>
      <td>8</td>
      <td>0</td>
      <td>2082</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>653</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>457</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.65 vs desc_0.7
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>46044</td>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>15259</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>15796</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>10193</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6904</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6498</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6198</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6480</td>
      <td>1</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3698</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4107</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>5274</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1070</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2561</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>654</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>457</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.7 vs desc_0.75
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>mic</th>
      <th>glu</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>46043</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>15255</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>11</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>15796</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>10189</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>6</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6904</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6498</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6198</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6483</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3699</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4107</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>6</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
      <td>5266</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1070</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>12</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>2553</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>654</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>457</td>
    </tr>
  </tbody>
</table>
</div>


    desc_0.75 vs desc_0.8
    


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
      <th>oli</th>
      <th>glu</th>
      <th>ast</th>
      <th>glu</th>
      <th>opc</th>
      <th>gab</th>
      <th>glu</th>
      <th>gab</th>
      <th>glu</th>
      <th>mic</th>
      <th>gab</th>
      <th>glu</th>
      <th>glu</th>
      <th>oli</th>
      <th>oli</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>oli</th>
      <td>45929</td>
      <td>52</td>
      <td>11</td>
      <td>11</td>
      <td>2</td>
      <td>20</td>
      <td>11</td>
      <td>5</td>
      <td>1</td>
      <td>4</td>
      <td>4</td>
      <td>3</td>
      <td>1</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>28</td>
      <td>15158</td>
      <td>3</td>
      <td>0</td>
      <td>5</td>
      <td>3</td>
      <td>49</td>
      <td>3</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>8</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>3</td>
      <td>0</td>
      <td>15755</td>
      <td>2</td>
      <td>4</td>
      <td>2</td>
      <td>8</td>
      <td>8</td>
      <td>0</td>
      <td>1</td>
      <td>12</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>11</td>
      <td>2</td>
      <td>4</td>
      <td>10047</td>
      <td>0</td>
      <td>5</td>
      <td>3</td>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>5</td>
      <td>112</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>6891</td>
      <td>1</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>3</td>
      <td>6</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>13</td>
      <td>13</td>
      <td>4169</td>
      <td>0</td>
      <td>1</td>
      <td>2286</td>
      <td>4</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>6120</td>
      <td>7</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>66</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>11</td>
      <td>134</td>
      <td>6</td>
      <td>2</td>
      <td>0</td>
      <td>15</td>
      <td>6267</td>
      <td>19</td>
      <td>5</td>
      <td>1</td>
      <td>6</td>
      <td>39</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>4</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>4</td>
      <td>18</td>
      <td>4</td>
      <td>3658</td>
      <td>0</td>
      <td>2</td>
      <td>1</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>5</td>
      <td>0</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>5</td>
      <td>3</td>
      <td>7</td>
      <td>1</td>
      <td>4076</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>2</td>
      <td>22</td>
      <td>5</td>
      <td>46</td>
      <td>0</td>
      <td>2</td>
      <td>54</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>5135</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>glu</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>1061</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>679</td>
      <td>84</td>
      <td>20</td>
      <td>8</td>
      <td>1</td>
      <td>23</td>
      <td>101</td>
      <td>46</td>
      <td>4</td>
      <td>24</td>
      <td>14</td>
      <td>9</td>
      <td>11</td>
      <td>1527</td>
      <td>4</td>
      <td>1</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>654</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>454</td>
    </tr>
  </tbody>
</table>
</div>



```python
celltype_idx=['glu','oli','gab','ast','mic','opc','end']
cells=adata.obs.keys()[(adata.obs.keys().str.startswith('celltype'))]
nn=len(res)
adata.obs['final_celltype']='none'
celltype=celltype_idx
for i in range(0,len(celltype)):
    if (celltype[i]!='opc')&(celltype[i]!='end'):
        adata.obs['final_celltype'][np.sum(adata.obs[cells]==celltype[i],axis=1)==nn]=celltype[i]
    else:
        adata.obs['final_celltype'][np.sum(adata.obs[cells]==celltype[i],axis=1)==nn-1]=celltype[i]
display(adata.obs['final_celltype'].value_counts())
```

    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:10: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      # Remove the CWD from sys.path while we load stuff.
    C:\Users\wangk\Anaconda3\lib\site-packages\ipykernel_launcher.py:14: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      
    


    oli     44182
    glu     39176
    ast     15683
    gab     12286
    none     8633
    opc      6853
    mic      3982
    end       444
    Name: final_celltype, dtype: int64



```python
tab=pd.DataFrame(adata.obs['final_celltype'].value_counts())
tab.final_celltype.sum(axis=0)-tab.loc['none']
```




    final_celltype    122606
    Name: none, dtype: int64




```python
for ires in [0.2,0.4,0.6,0.8]:
    print(ires)
    adata.obs['maxprob']=(np.max(adata.uns['prob_matrix'+str(ires)],axis=1)*1000).astype(int)
    adata.obsm.X_tsne=adata.obsm['X_tsne'+str(ires)]
    sc.pl.tsne(adata, color=['final_celltype','n_counts','maxprob'])
```

    0.2
    


![png](output_16_1.png)


    0.4
    


![png](output_16_3.png)


    0.6
    


![png](output_16_5.png)


    0.8
    


![png](output_16_7.png)



```python
for i in range(0,len(res)-1):
    print("celltype_"+res[i]+" vs celltype_"+res[i+1])
    tab=confusion_matrix(adata.obs['celltype_'+res[i]],adata.obs['celltype_'+res[i+1]],labels=celltype_idx)
    display(pd.DataFrame(tab,index=celltype_idx,columns=celltype_idx))
```

    celltype_0.1 vs celltype_0.15
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>39298</td>
      <td>42</td>
      <td>348</td>
      <td>19</td>
      <td>13</td>
      <td>5</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>99</td>
      <td>45820</td>
      <td>151</td>
      <td>23</td>
      <td>185</td>
      <td>8854</td>
      <td>462</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>90</td>
      <td>9</td>
      <td>12745</td>
      <td>1</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>16</td>
      <td>29</td>
      <td>67</td>
      <td>16219</td>
      <td>2</td>
      <td>14</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>30</td>
      <td>877</td>
      <td>43</td>
      <td>9</td>
      <td>5754</td>
      <td>10</td>
      <td>2</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.15 vs celltype_0.2
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>39391</td>
      <td>69</td>
      <td>34</td>
      <td>7</td>
      <td>22</td>
      <td>10</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>441</td>
      <td>46213</td>
      <td>46</td>
      <td>32</td>
      <td>22</td>
      <td>17</td>
      <td>6</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>642</td>
      <td>79</td>
      <td>12624</td>
      <td>1</td>
      <td>6</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>155</td>
      <td>104</td>
      <td>28</td>
      <td>15956</td>
      <td>23</td>
      <td>5</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>59</td>
      <td>1805</td>
      <td>6</td>
      <td>4</td>
      <td>4080</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>72</td>
      <td>1787</td>
      <td>23</td>
      <td>14</td>
      <td>10</td>
      <td>6978</td>
      <td>1</td>
    </tr>
    <tr>
      <th>end</th>
      <td>4</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>453</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.2 vs celltype_0.25
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>40750</td>
      <td>11</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>12</td>
      <td>50047</td>
      <td>1</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>2</td>
      <td>3</td>
      <td>12755</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>2</td>
      <td>7</td>
      <td>0</td>
      <td>16004</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>4159</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>7007</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>462</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.25 vs celltype_0.3
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>40681</td>
      <td>68</td>
      <td>7</td>
      <td>4</td>
      <td>1</td>
      <td>6</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>231</td>
      <td>49808</td>
      <td>21</td>
      <td>7</td>
      <td>5</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>106</td>
      <td>25</td>
      <td>12623</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>97</td>
      <td>30</td>
      <td>13</td>
      <td>15867</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>31</td>
      <td>15</td>
      <td>2</td>
      <td>3</td>
      <td>4109</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>16</td>
      <td>62</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>6928</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>458</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.3 vs celltype_0.35
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41063</td>
      <td>51</td>
      <td>29</td>
      <td>11</td>
      <td>8</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>252</td>
      <td>49667</td>
      <td>37</td>
      <td>20</td>
      <td>29</td>
      <td>3</td>
      <td>2</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>21</td>
      <td>13</td>
      <td>12626</td>
      <td>1</td>
      <td>2</td>
      <td>4</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>17</td>
      <td>16</td>
      <td>2</td>
      <td>15845</td>
      <td>3</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>12</td>
      <td>4</td>
      <td>4</td>
      <td>1</td>
      <td>4095</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>20</td>
      <td>14</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>6904</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>457</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.35 vs celltype_0.4
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41303</td>
      <td>38</td>
      <td>25</td>
      <td>10</td>
      <td>3</td>
      <td>6</td>
      <td>1</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>429</td>
      <td>49207</td>
      <td>65</td>
      <td>15</td>
      <td>35</td>
      <td>13</td>
      <td>1</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>18</td>
      <td>17</td>
      <td>12658</td>
      <td>2</td>
      <td>3</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>11</td>
      <td>26</td>
      <td>3</td>
      <td>15835</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>12</td>
      <td>16</td>
      <td>6</td>
      <td>0</td>
      <td>4103</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>6</td>
      <td>4</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>6903</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>456</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.4 vs celltype_0.45
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41773</td>
      <td>3</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>10</td>
      <td>49298</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>12757</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>15862</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4144</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6925</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.45 vs celltype_0.5
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41781</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>4</td>
      <td>49295</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>1</td>
      <td>1</td>
      <td>12757</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>15863</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4146</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6926</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.5 vs celltype_0.55
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41780</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>10</td>
      <td>49285</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>1</td>
      <td>1</td>
      <td>12756</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>15862</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>4146</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6926</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.55 vs celltype_0.6
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41766</td>
      <td>22</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>3</td>
      <td>49291</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>1</td>
      <td>12756</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>15859</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>4143</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6926</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>459</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.6 vs celltype_0.65
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41618</td>
      <td>119</td>
      <td>16</td>
      <td>5</td>
      <td>8</td>
      <td>5</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>226</td>
      <td>49063</td>
      <td>12</td>
      <td>9</td>
      <td>8</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>69</td>
      <td>34</td>
      <td>12654</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>54</td>
      <td>16</td>
      <td>10</td>
      <td>15777</td>
      <td>2</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>21</td>
      <td>28</td>
      <td>2</td>
      <td>2</td>
      <td>4089</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>11</td>
      <td>16</td>
      <td>2</td>
      <td>4</td>
      <td>0</td>
      <td>6893</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>457</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.65 vs celltype_0.7
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41996</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>5</td>
      <td>49271</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>1</td>
      <td>0</td>
      <td>12696</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>15796</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>4107</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6904</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>457</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.7 vs celltype_0.75
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>42000</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>7</td>
      <td>49265</td>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>0</td>
      <td>0</td>
      <td>12696</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>15796</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4107</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>6904</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>457</td>
    </tr>
  </tbody>
</table>
</div>


    celltype_0.75 vs celltype_0.8
    


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
      <th>glu</th>
      <th>oli</th>
      <th>gab</th>
      <th>ast</th>
      <th>mic</th>
      <th>opc</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>glu</th>
      <td>41835</td>
      <td>58</td>
      <td>85</td>
      <td>20</td>
      <td>2</td>
      <td>6</td>
      <td>1</td>
    </tr>
    <tr>
      <th>oli</th>
      <td>296</td>
      <td>48796</td>
      <td>112</td>
      <td>31</td>
      <td>28</td>
      <td>3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>gab</th>
      <td>33</td>
      <td>6</td>
      <td>12656</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ast</th>
      <td>10</td>
      <td>5</td>
      <td>22</td>
      <td>15755</td>
      <td>1</td>
      <td>4</td>
      <td>0</td>
    </tr>
    <tr>
      <th>mic</th>
      <td>5</td>
      <td>6</td>
      <td>19</td>
      <td>2</td>
      <td>4076</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>opc</th>
      <td>2</td>
      <td>2</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
      <td>6891</td>
      <td>0</td>
    </tr>
    <tr>
      <th>end</th>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>454</td>
    </tr>
  </tbody>
</table>
</div>



```python
adata.write("brain.h5ad")
```


```python
from matplotlib.cm import get_cmap
import matplotlib as mpl
cmap_tab10 = get_cmap('tab10')
cmap = mpl.colors.ListedColormap([cmap_tab10.colors[0],cmap_tab10.colors[3],cmap_tab10.colors[2],cmap_tab10.colors[4],cmap_tab10.colors[1],cmap_tab10.colors[8],cmap_tab10.colors[9]])
```


```python
ires='0.15'
adata.obsm['X_tsne']=adata.obsm['X_tsne'+ires+'_600_120']
ncols=1
nrows=1
fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[ncols*10,nrows*10])
ax=sc.pl.scatter(adata,basis='tsne',
                 color='final_celltype',
                 legend_loc='on data',
                 ax=axes,
                 size=10,
                 color_map=cmap,
                 show=False)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
```


![png](output_20_0.png)

