from __future__ import annotations
import warnings
from typing import Dict, Tuple, Any
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import chi2
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, permdisp, DistanceMatrix
from skbio.diversity.alpha import shannon, simpson, chao1


# ─────────────────────────────────────────────────────────────
# Helper: variance barplot (top n + "Others")
# ─────────────────────────────────────────────────────────────
def _plot_pcoa_variance_topn(
    pcoa_res,
    n: int = 5,
    fig_path: str = "pcoa_variance_topn.png",
    show_plot: bool = False,
    title: str = "",
) -> pd.DataFrame:
    """
    Build a small table of variance explained per axis, aggregating
    axes after the first n as 'Others'. Optionally plot a barplot.

    Returns
    -------
    var_table : DataFrame
        Index = axis labels (e.g. 'PC1', 'PC2', ..., 'Others')
        Column = 'variance_explained' (%)
    """
    var = pcoa_res.proportion_explained.copy()  # Series
    var_pct = 100 * var

    if len(var_pct) > n:
        top = var_pct.iloc[:n]
        others = var_pct.iloc[n:].sum()
        var_plot = pd.concat([top, pd.Series({"Others": others})])
    else:
        var_plot = var_pct

    var_table = var_plot.to_frame(name="variance_explained")

    if show_plot:
        plt.figure(figsize=(6, 4))
        ax = sns.barplot(
            x=var_table.index,
            y="variance_explained",
            data=var_table.reset_index()
        )
        ax.set_ylabel("Variance explained (%)")
        ax.set_xlabel("Axis")
        ax.set_title(title if title else "PCoA variance explained")
        plt.tight_layout()
        plt.savefig(fig_path, dpi=300)
        plt.show()
        plt.close()
    else:
        # If you prefer to always save the figure, even without showing:
        # plt.figure(...)
        # ...
        # plt.savefig(fig_path, dpi=300)
        # plt.close()
        pass

    return var_table


# ─────────────────────────────────────────────────────────────
# Helper: centroids + confidence ellipses
# ─────────────────────────────────────────────────────────────
def _plot_centroids_and_ellipses(
    ax: plt.Axes,
    coords_meta: pd.DataFrame,
    hue_col: str,
    hue_order: list[str],
    palette: Dict[str, Any],
    level: float = 0.95,
) -> Dict[str, Tuple[float, float]]:
    """
    Add centroids + confidence ellipses for each group on an existing PCoA plot.

    Parameters
    ----------
    ax : matplotlib Axes
        Axes with the scatter already drawn.
    coords_meta : DataFrame
        Must contain 'PC1', 'PC2' and `hue_col`.
    hue_col : str
        Column name with group labels.
    hue_order : list of str
        Order of groups (for consistent colors).
    palette : dict
        Mapping group -> color.
    level : float
        Confidence level (e.g. 0.95).

    Returns
    -------
    centroids : dict
        {group: (mean_PC1, mean_PC2)}
    """
    centroids: Dict[str, Tuple[float, float]] = {}

    # Chi-square cutoff for the confidence level with 2 DoF
    chi2_val = chi2.ppf(level, df=2)

    for group in hue_order:
        subset = coords_meta.loc[
            coords_meta[hue_col] == group, ["PC1", "PC2"]
        ].dropna()

        if subset.shape[0] == 0:
            continue

        x = subset["PC1"].values
        y = subset["PC2"].values

        mean_x = x.mean()
        mean_y = y.mean()
        centroids[group] = (mean_x, mean_y)

        # Need at least 2 points for ellipse (covariance)
        if subset.shape[0] < 2:
            # Still plot the centroid
            ax.scatter(
                mean_x, mean_y,
                marker="X",
                s=120,
                color=palette[group],
                edgecolor="black",
                linewidth=1,
                zorder=3,
            )
            continue

        # 2x2 covariance matrix
        cov = np.cov(x, y)
        eigvals, eigvecs = np.linalg.eigh(cov)
        order = eigvals.argsort()[::-1]
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]

        # Ellipse width/height (2 * sqrt(lambda * chi2_val))
        width, height = 2 * np.sqrt(eigvals * chi2_val)

        # Angle of the largest eigenvector
        angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))

        color = palette[group]

        ellip = Ellipse(
            (mean_x, mean_y),
            width=width,
            height=height,
            angle=angle,
            edgecolor=color,
            facecolor=color,
            alpha=0.2,
            lw=1.5,
            zorder=2,
        )
        ax.add_patch(ellip)

        # Plot centroid as a cross
        ax.scatter(
            mean_x, mean_y,
            marker="X",
            s=120,
            color=color,
            edgecolor="black",
            linewidth=1,
            zorder=3,
        )

    return centroids


# ─────────────────────────────────────────────────────────────
# Main function: PCoA + centroids/ellipses + PERMANOVA/PERMDISP
# ─────────────────────────────────────────────────────────────
def make_pcoa3(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    hue_col: str,
    fig_path: str,
    metric: str = "euclidean",
    intersect_dfs: bool = False,
    show_plot: bool = False,
    title: str = "",
    # Variance barplot options
    varplot_n: int = 5,
    varplot_path: str = "pcoa_variance_topn.png",
    varplot_title: str = "",
    # Centroids + ellipses
    add_centroids: bool = False,
    ellipse_level: float = 0.95,
    # PERMANOVA / PERMDISP
    run_stats: bool = False,
    n_perm: int = 999,
    stats_seed: int | None = None,
) -> Dict[str, Any]:
    """
    Run PCoA from an abundance table, plot PC1 vs PC2, optionally add
    centroids + confidence ellipses and run PERMANOVA / PERMDISP.

    Returns
    -------
    result : dict
        {
          "PC1": Series,
          "PC2": Series,
          "variance_table": DataFrame,
          "pcoa_res": OrdinationResults,
          "centroids": dict or None,
          "permanova": Series or None,
          "permdisp": Series or None,
        }
    """

    # 1) Optional: align samples
    if intersect_dfs:
        shared = abundance_df.index.intersection(metadata_df.index)
        abundance_df = abundance_df.loc[shared].copy()
        metadata_df = metadata_df.loc[shared].copy()
    else:
        abundance_df = abundance_df.copy()
        metadata_df = metadata_df.copy()

    # Ensure string IDs (good practice for skbio)
    abundance_df.index = abundance_df.index.astype(str)
    metadata_df.index = metadata_df.index.astype(str)

    # 2) distance matrix
    METRICS_REQUIRE_NONNEG = {
        "braycurtis",
        "jaccard",
        "unweighted_unifrac",
        "weighted_unifrac",
    }
    metric = metric.lower()
    validate = metric in METRICS_REQUIRE_NONNEG

    dm: DistanceMatrix = beta_diversity(
        metric,
        abundance_df.values,
        ids=abundance_df.index,
        validate=validate,
    )

    # 3) PCoA
    pcoa_res = pcoa(dm)

    # 4) Save coordinates
    coords = pcoa_res.samples
    coords.to_csv("pcoa_coordinates.tsv", sep="\t", header=True)

    # 5) Merge with metadata
    coords_meta = coords.join(metadata_df)

    centroids = None  # will be filled if add_centroids & show_plot

    # 6) Scatter of PC1 vs PC2
    if show_plot:
        plt.figure(figsize=(8, 6))

        # consistent mapping between groups and colors
        hue_order = sorted(coords_meta[hue_col].dropna().unique())
        raw_palette = sns.color_palette("Set2", n_colors=len(hue_order))
        palette = dict(zip(hue_order, raw_palette))

        ax = sns.scatterplot(
            data=coords_meta,
            x="PC1",
            y="PC2",
            hue=hue_col,
            hue_order=hue_order,
            palette=palette,
            s=100,
            edgecolor="black",
        )

        var1 = pcoa_res.proportion_explained[0] * 100
        var2 = pcoa_res.proportion_explained[1] * 100
        ax.set_xlabel(f"PC1 ({var1:.1f}%)")
        ax.set_ylabel(f"PC2 ({var2:.1f}%)")
        ax.set_title(title if title else f"PCoA ({metric}) colored by {hue_col}")
        ax.grid(True)
        ax.legend(title=hue_col)

        if add_centroids:
            centroids = _plot_centroids_and_ellipses(
                ax=ax,
                coords_meta=coords_meta,
                hue_col=hue_col,
                hue_order=hue_order,
                palette=palette,
                level=ellipse_level,
            )

        plt.tight_layout()
        plt.savefig(fig_path, dpi=300)
        plt.show()
        plt.close()

    # 7) PERMANOVA + PERMDISP
    permanova_res = None
    permdisp_res = None

    if run_stats:
        # Use DataFrame interface so IDs are handled by skbio
        grouping_df = metadata_df[[hue_col]].copy()

        # Keep only IDs that are in the distance matrix
        grouping_df = grouping_df.loc[list(dm.ids)]

        # Drop missing groups
        mask = grouping_df[hue_col].notna()
        grouping_df = grouping_df.loc[mask]

        if grouping_df.shape[0] < 2:
            warnings.warn(
                "Not enough samples with non-missing group labels for PERMANOVA/PERMDISP."
            )
        else:
            # Filter distance matrix to valid IDs
            dm_stats = dm.filter(grouping_df.index)

            # Check we have at least two groups
            if grouping_df[hue_col].nunique() < 2:
                warnings.warn(
                    "PERMANOVA/PERMDISP require at least two groups; skipping tests."
                )
            else:
                try:
                    permanova_res = permanova(
                        dm_stats,
                        grouping=grouping_df,
                        column=hue_col,
                        permutations=n_perm,
                        seed=stats_seed,
                    )
                except Exception as e:
                    warnings.warn(f"PERMANOVA failed: {e}")
                    permanova_res = None

                try:
                    permdisp_res = permdisp(
                        dm_stats,
                        grouping=grouping_df,
                        column=hue_col,
                        permutations=n_perm,
                        seed=stats_seed,
                    )
                except Exception as e:
                    warnings.warn(f"PERMDISP failed: {e}")
                    permdisp_res = None

    # 8) Variance barplot: first n + remainder
    var_table = _plot_pcoa_variance_topn(
        pcoa_res,
        n=varplot_n,
        fig_path=varplot_path,
        show_plot=show_plot,
        title=varplot_title,
    )

    return {
        "PC1": pcoa_res.samples["PC1"],
        "PC2": pcoa_res.samples["PC2"],
        "variance_table": var_table,
        "pcoa_res": pcoa_res,
        "centroids": centroids,
        "permanova": permanova_res,
        "permdisp": permdisp_res,
    }


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


def make_pcoa2(
    abundance_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    hue_col: str,
    fig_path: str,
    metric: str = "euclidean",
    intersect_dfs: bool = False,
    show_plot: bool = False,
    title: str = "",
    # NEW (existing)
    varplot_n: int = 5,
    varplot_path: str = "pcoa_variance_topn.png",
    varplot_title: str = "",
    # NEW (centroids)
    add_centroids: bool = True,
    centroid_marker: str = "X",
    centroid_size: int = 250,
):
    # 1) Optional: align samples
    if intersect_dfs:
        shared = abundance_df.index.intersection(metadata_df.index)
        abundance_df = abundance_df.loc[shared]
        metadata_df  = metadata_df.loc[shared]

    # 2) distance matrix
    METRICS_REQUIRE_NONNEG = {
        "braycurtis",
        "jaccard",
        "unweighted_unifrac",
        "weighted_unifrac",
    }
    metric = metric.lower()
    validate = metric in METRICS_REQUIRE_NONNEG

    dm = beta_diversity(
        metric,
        abundance_df.values,
        ids=abundance_df.index.astype(str),
        validate=validate,
    )

    # 3) PCoA
    pcoa_res = pcoa(DistanceMatrix(dm, ids=abundance_df.index.astype(str)))

    # 4) save coordinates
    coords = pcoa_res.samples
    coords.to_csv("pcoa_coordinates.tsv", sep="\t", header=True)

    # 5) merge with metadata (if you want the scatter)
    coords_meta = coords.join(metadata_df)

    # Precompute variance explained for labels
    var1 = pcoa_res.proportion_explained[0] * 100
    var2 = pcoa_res.proportion_explained[1] * 100

    # 6) scatter of PC1 vs PC2 (optional)
    centroids = None
    if show_plot:
        plt.figure(figsize=(8, 6))
        ax = sns.scatterplot(
            data=coords_meta,
            x="PC1", y="PC2",
            hue=hue_col,
            palette="Set2",
            s=100,
            edgecolor="black",
        )

        # --- NEW: compute and plot centroids ---
        if add_centroids:
            centroids = (
                coords_meta
                .groupby(hue_col)[["PC1", "PC2"]]
                .mean()
                .reset_index()
            )

            # Plot centroids as larger markers
            sns.scatterplot(
                data=centroids,
                x="PC1", y="PC2",
                hue=hue_col,
                palette="Set2",
                s=centroid_size,
                marker=centroid_marker,
                edgecolor="black",
                legend=False,  # keep legend from first scatter
                ax=ax,
            )

            # Optional: label centroids with group name
            for _, row in centroids.iterrows():
                ax.text(
                    row["PC1"],
                    row["PC2"],
                    str(row[hue_col]),
                    ha="center",
                    va="center",
                    fontsize=9,
                    weight="bold",
                )

        ax.set_xlabel(f"PC1 ({var1:.1f}%)")
        ax.set_ylabel(f"PC2 ({var2:.1f}%)")
        ax.set_title(title if title else f"PCoA ({metric}) colored by {hue_col}")
        ax.grid(True)
        ax.legend(title=hue_col)
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
        title=varplot_title,
    )

    # Return what’s useful to reuse downstream
    return {
        "PC1": pcoa_res.samples["PC1"],
        "PC2": pcoa_res.samples["PC2"],
        "variance_table": var_table,
        "pcoa_res": pcoa_res,
        "centroids": centroids,  # NEW: None if add_centroids=False
    }


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
