ml  anaconda3/2023.09
source activate base
conda activate /sc/arion/projects/adineto/py312/scanpy_env_apoe
conda list | grep 
mkdir results

python

# the python enviroment
import numpy as np
import pandas as pd
import scanpy as sc
#import pooch
import anndata as ad

sc.settings.verbosity = 3  
# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

adata = sc.read_h5ad("GSE254205_ad_raw.h5ad")
adata
adata.obs_names_make_unique()

# Quality Control
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,stripplot=False
)

from matplotlib import pyplot as plt

with plt.rc_context({"figure.figsize": (5, 5)}):  # Use this to set figure params like size and dpi
    sc.pl.violin(adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,stripplot=False, show=False)
    plt.savefig("plots/violinplots.pdf", bbox_inches="tight")

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# output cell meta data
#meta = adata.obs
#pd.DataFrame(meta).to_csv("sample_meta.csv")
#or 
#adata.obs.to_csv("results/sample_meta.csv")


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Doublet detection
#from scipy.sparse import csr_matrix # do not need to use
#adata.X = csr_matrix(adata.X)      # do not to use
sc.pp.scrublet(adata, batch_key="sample")


# Saving count data
adata.layers["counts"] = adata.X.copy()
adata.obs.to_csv("results/sample_meta_dbl.csv")
adata.write_h5ad("results/GSE254205_v1.h5ad")

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(adata)

#Dimensionality Reduction
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)
# Nearest neighbor graph constuction and visualization
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)

# Clustering
# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"])

sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
)


sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

#Manual cell-type annotation
for res in [0.02,0.25, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
    
sc.pl.umap(
    adata,
    color=["leiden_res_0.02","leiden_res_0.25", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)


sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=0.25)
with plt.rc_context({"figure.figsize": (3.5, 3.5)}):  # Use this to set figure params like size and dpi
    sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=1, show = False,
    )
    plt.savefig("plots/doubletplots.pdf", bbox_inches="tight")


#Marker gene set
marker_genes = {
    "Ast": ["AQP4","GJA1","GFAP"],
    "Mic": ["CSF1R","CD74","P2RY12"],
    
    "Endo": ["FLT1","APOLD1", "CD34"],
    "Ex": ["NRGN","CCK","GRIN1","GRIN2B","SLC17A7","SLC17A6","GLUL"],
    
    "In": ["GAD1","GAD2","GABBR1","GABBR2","SLC6A1","SLC32A1"],
    "Olig": ["MAG","MOG","MOBP","UGT8","ERMN","MBP"],
    "Opc": ["PDGFRA","VCAN","CA10"],
    "Neu": ["STMN2","VGF","SYNPR"]}
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.25", standard_scale="var")
adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.25"].map(
    {
        "0": "Olig",
        "1": "Olig",
        "2": "Ast",
        "3": "In",
        "4": "Lymphocytes",
        "5": "Ex",
        "6": "Ex",
        "7": "Ex",
        "8": "Ex",
        "9": "Endo",
        "10": "Opc",
        "11": "In",
        "12": "Ex",
        "13": "Ex",
        "14": "Ex",
        "15": "Un",
    }
)

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.25", standard_scale="var")
adata.write_h5ad("results/GSE254205_v1.h5ad")

