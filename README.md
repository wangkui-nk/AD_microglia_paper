# APOE and TREM2 regulate amyloid responsive microglia (ARM) in Alzheimer’s disease


Analysis | tools
------------ | -------------
Cluster analysis | [desc](https://github.com/eleozzr/desc)
Trajectory analysis  | [paga](https://github.com/theislab/paga)
Trajectory analysis  | [monocle3](https://cole-trapnell-lab.github.io/monocle3/)

**Clustering:**
The nuclei were clustered based on the 2,000 highly variable genes using [desc](https://github.com/eleozzr/desc).We used two hidden layers for encoder with 256 nodes in the first layer 
and 128 nodes in the second layer.To identify an optimal resolution parameter, we performed the DESC analysis 15 times 
using different resolutions ranging from 0.1 to 0.8 with an incremental 0.05 step. 

[**Batch effect analysis**](docs/batch_effect_analysis.md)

[**Cell type annotation**](docs/annotation_cell_type.md]

**Microglia recluster:**
To identify microglia subpopulations, microglia cells was re-clustered by [desc](https://github.com/eleozzr/desc). We 
selected 5,000 highly variable genes in order to detect subtle cellular differences among the 
microglia. A single 128-node layer was applied for the encoder in DESC and the resolutions 
ranged from 0.1 to 0.8 with an incremental 0.05 increasing step. 

[**Microglia contamination analysis**](docs/contamination.md)

**Trajectory analysis:**
[**Monocle3 analysis**](docs/monocle3.md)
[**Paga analysis**](docs/PAGA_analysis.md)


# Reference:
1. Xiangjie Li, Kui Wang, Yafei Lyu, Jingxiao Zhang, Dwight Stambolian, Katalin Susztak, Muredach P. Reilly, Gang
Hu, Mingyao Li. Deep learning enables accurate clustering with batch effect removal in single-cell RNA-seq analysis. Nature Communications.

