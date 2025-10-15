import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from skbio.diversity.alpha import shannon, simpson, chao1
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

def make_pca(abundance_df: pd.DataFrame, metadata_df: pd.DataFrame, hue_col: str, fig_path: str, interecept_dfs: bool = False, show_plot=False, ):
    
    if interecept_dfs:
        shared_index = abundance_df.index.intersection(metadata_df.index)
        abundance_df = abundance_df.loc[shared_index]
        metadata_df = metadata_df.loc[shared_index]

    # Compute Bray–Curtis distance matrix
    bc_dm = beta_diversity("braycurtis", abundance_df.values, ids=abundance_df.index)

    # Perform PCoA
    pcoa_results = pcoa(bc_dm)

    from scipy.stats import spearmanr

    # Correlate each taxon's abundance with PC1
    pcoa1_scores = pcoa_results.samples['PC1']
    pcoa2_scores = pcoa_results.samples['PC2']

    # Save coordinates to file
    pcoa_results.samples.to_csv("pcoa_coordinates.tsv", sep="\t")

    # Merge coordinates with metadata
    coords = pcoa_results.samples
    coords_with_meta = coords.join(metadata_df)

    if show_plot:
        # Plot with seaborn
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            data=coords_with_meta,
            x="PC1",
            y="PC2",
            hue=hue_col,
            palette="Set2",
            s=100,
            edgecolor="black"
        )
        plt.xlabel(f"PC1 ({pcoa_results.proportion_explained[0]*100:.1f}%)")
        plt.ylabel(f"PC2 ({pcoa_results.proportion_explained[1]*100:.1f}%)")
        plt.title(f"PCoA colored by {hue_col}")
        plt.grid(True)
        plt.legend(title=hue_col)
        plt.tight_layout()
        plt.savefig(fig_path)
        plt.show()

    return pcoa1_scores, pcoa2_scores

def plot_sample_richness(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    hue_col: str,
    fig_path: str,
    interecept_dfs: bool = False,
    metric: str = "observed",
    show_plot = False
):
    """
    Plot richness or alpha diversity per sample group.

    Parameters:
    - abundance_df: DataFrame with samples as rows and taxa/features as columns.
    - metadata_df: DataFrame with samples as index and metadata columns.
    - hue_col: Column in metadata_df to group samples.
    - interecept_dfs: If True, intersect the index of both dataframes.
    - metric: One of {"observed", "shannon", "simpson", "chao1"}.
    """
    if interecept_dfs:
        shared_index = abundance_df.index.intersection(metadata_df.index)
        abundance_df = abundance_df.loc[shared_index]
        metadata_df = metadata_df.loc[shared_index]

    # Select metric
    if metric == "observed":
        richness = (abundance_df > 0).sum(axis=1)
    elif metric == "shannon":
        richness = abundance_df.apply(shannon, axis=1)
    elif metric == "simpson":
        richness = abundance_df.apply(simpson, axis=1)
    elif metric == "chao1":
        richness = abundance_df.apply(chao1, axis=1)
    else:
        raise ValueError(f"Unsupported metric '{metric}'. Choose from: observed, shannon, simpson, chao1.")

    richness.name = metric

    # Merge with metadata
    richness_df = richness.to_frame().join(metadata_df)

    # Plot
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=richness_df, x=hue_col, y=metric, palette="Set2")
    sns.stripplot(data=richness_df, x=hue_col, y=metric, color="black", alpha=0.5, jitter=True)
    plt.title(f"{metric.capitalize()} by {hue_col}")
    plt.ylabel(f"{metric.capitalize()} index")
    plt.xlabel(hue_col)
    plt.tight_layout()
    plt.ylim(bottom=0)
    plt.grid(True, axis="y")
    plt.savefig(fig_path)
    if show_plot:
        plt.show()
    else:
        plt.close()

    return richness_df

def _plot_pcoa_variance_topn(pcoa_res,
                             n: int = 5,
                             fig_path: str = "pcoa_variance_topn.png",
                             show_plot: bool = True,
                             title: str | None = None):
    """
    Barplot of the PCoA variance explained: first n PCs individually,
    remaining PCs grouped into a single bar.
    """
    # proportions as percentages
    var_exp = np.asarray(pcoa_res.proportion_explained) * 100.0

    # guardrails
    if n < 1:
        n = 1
    n = min(n, len(var_exp))

    # split top-n and remainder
    top = var_exp[:n]
    rest = var_exp[n:]
    labels = [f"PC{i+1}" for i in range(n)]
    values = list(top)
    if rest.size > 0:
        labels.append(f"PC{n+1}+")
        values.append(rest.sum())

    # barplot
    plt.figure(figsize=(8, 5))
    bars = plt.bar(labels, values, edgecolor="black")
    plt.ylabel("Variance Explained (%)")
    plt.xlabel("Principal Coordinate")
    if title:
        plt.title(title)
    else:
        plt.title(f"PCoA variance explained (first {n} PCs + remainder)")

    # add data labels
    for rect, v in zip(bars, values):
        plt.text(rect.get_x() + rect.get_width()/2.0, v + 0.6, f"{v:.1f}%",
                 ha="center", va="bottom")

    plt.ylim((0,100))
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    if show_plot:
        plt.show()
    plt.close()

    # return a tidy table of what was plotted (useful in reports)
    out = pd.DataFrame({"component": labels, "variance_pct": values})
    return out

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix
from skbio.diversity import beta_diversity

def make_pcoa(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    hue_col: str,
    fig_path: str,
    metric: str = "euclidean",
    intersect_dfs: bool = False,
    show_plot: bool = False,
    title: str = "",
    # NEW:
    varplot_n: int = 5,
    varplot_path: str = "pcoa_variance_topn.png",
    varplot_title: str = ""
):
    # 1) Optional: align samples
    if intersect_dfs:
        shared = abundance_df.index.intersection(metadata_df.index)
        abundance_df = abundance_df.loc[shared]
        metadata_df  = metadata_df.loc[shared]

    # 2) distance matrix
    METRICS_REQUIRE_NONNEG = {"braycurtis", "jaccard", "unweighted_unifrac", "weighted_unifrac"}
    metric = metric.lower()
    if metric in METRICS_REQUIRE_NONNEG:
        validate = True
    else:
        validate = False

    dm = beta_diversity(metric, abundance_df.values, ids=abundance_df.index.astype(str), validate=validate)

    # 3) PCoA
    pcoa_res = pcoa(DistanceMatrix(dm, ids=abundance_df.index.astype(str)))

    # 4) save coordinates
    coords = pcoa_res.samples
    coords.to_csv("pcoa_coordinates.tsv", sep="\t", header=True)

    # 5) merge with metadata (if you want the scatter)
    coords_meta = coords.join(metadata_df)

    # 6) scatter of PC1 vs PC2 (optional)
    if show_plot:
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            data=coords_meta,
            x="PC1", y="PC2",
            hue=hue_col,
            palette="Set2",
            s=100, edgecolor="black"
        )
        var1 = pcoa_res.proportion_explained[0]*100
        var2 = pcoa_res.proportion_explained[1]*100
        plt.xlabel(f"PC1 ({var1:.1f}%)")
        plt.ylabel(f"PC2 ({var2:.1f}%)")
        plt.title(title if title else f"PCoA ({metric}) colored by {hue_col}")
        plt.grid(True)
        plt.legend(title=hue_col)
        plt.tight_layout()
        plt.savefig(fig_path, dpi=300)
        plt.show()
        plt.close()

    # 7) variance barplot: first n + remainder
    var_table = _plot_pcoa_variance_topn(
        pcoa_res,
        n=varplot_n,
        fig_path=varplot_path,
        show_plot=show_plot,
        title=varplot_title
    )

    # Return what’s useful to reuse downstream
    return {
        "PC1": pcoa_res.samples["PC1"],
        "PC2": pcoa_res.samples["PC2"],
        "variance_table": var_table,
        "pcoa_res": pcoa_res,
    }
