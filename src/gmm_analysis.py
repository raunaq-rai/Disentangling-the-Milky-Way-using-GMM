# gmm_analysis.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from tqdm import tqdm
import pygmmis
import pickle


def fit_gmm_fixed_components(df_bin, n_components, n_init=50):
    """
    Fit a Gaussian Mixture Model (GMM) to 3D velocity data using Extreme Deconvolution (XD).

    Parameters
    ----------
    df_bin : pd.DataFrame
        DataFrame containing velocity components and their uncertainties ('v_R', 'v_phi', 'v_Z').
    n_components : int
        Number of Gaussian components to fit.
    n_init : int, optional
        Number of initializations to avoid local minima (default: 50).

    Returns
    -------
    best_gmm : pygmmis.GMM
        Best-fitted GMM object with highest log-likelihood.
    """
    X = df_bin[['v_R', 'v_phi', 'v_Z']].values
    uncertainties = np.stack([
        df_bin['v_R_uncertainty'].values**2,
        df_bin['v_phi_uncertainty'].values**2,
        df_bin['v_Z_uncertainty'].values**2
    ], axis=1)
    cov_matrices = np.array([np.diag(u) for u in uncertainties])

    best_gmm = None
    best_logL = -np.inf

    for _ in tqdm(range(n_init), desc=f"Fitting GMM with {n_components} components"):
        gmm = pygmmis.GMM(K=n_components, D=3)
        pygmmis.initFromKMeans(gmm, X)
        logL, _ = pygmmis.fit(gmm, X, covar=cov_matrices, w=0.1, tol=1e-6, init_method='kmeans')

        if logL > best_logL:
            best_logL = logL
            best_gmm = gmm

    return best_gmm


def plot_gmm_with_contributions(df_bin, gmm, bin_label, bins=100, x_limits=(-400, 400), y_limits=(-400, 400),
                                component_colors=None, metallicity_range="", label="gmm_plot"):
    """
    Plot the velocity distribution (v_phi vs v_R) as a 2D histogram with GMM components shown as ellipses.
    Also visualizes the fractional contribution of each component.

    Parameters
    ----------
    df_bin : pd.DataFrame
        DataFrame containing velocity data ('v_R', 'v_phi').
    gmm : pygmmis.GMM
        Fitted GMM object.
    bin_label : str
        Title or label for the metallicity bin.
    bins : int, optional
        Number of bins in histogram (default: 100).
    x_limits : tuple
        Limits for the x-axis (v_R).
    y_limits : tuple
        Limits for the y-axis (v_phi).
    component_colors : list of str, optional
        List of colors to use for each GMM component.
    metallicity_range : str, optional
        Text describing the metallicity range (e.g., "-1.6 < [M/H] < -1.3").
    label : str, optional
        Filename label for saving the plot.

    Saves
    -----
    A figure is saved in ../figures/gmm_{label}.png
    """
    default_colors = ['blue', 'red', 'aqua', 'gold', 'green', 'purple', 'orange', 'pink']
    if component_colors is None:
        component_colors = default_colors
    elif len(component_colors) < gmm.K:
        raise ValueError(f"Not enough colors provided! Expected {gmm.K}, got {len(component_colors)}")

    fig, axes = plt.subplots(2, 1, figsize=(4, 6), gridspec_kw={'height_ratios': [1, 3]}, constrained_layout=True)

    df_filtered = df_bin.dropna(subset=["v_R", "v_phi"])
    H, xedges, yedges = np.histogram2d(df_filtered["v_R"], df_filtered["v_phi"],
                                       bins=bins, range=[x_limits, y_limits], density=True)

    H_min = np.min(H[H > 0]) if np.any(H > 0) else 1e-5
    H_max = np.max(H) if np.max(H) > 0 else 1
    H_normalized = (H - H_min) / (H_max - H_min)
    H_normalized = np.clip(H_normalized, 0, 1)

    ax = axes[1]
    im = ax.pcolormesh(xedges, yedges, H_normalized.T, cmap="Greys", shading='auto')

    total_weight = np.sum(gmm.amp)
    weights = (gmm.amp / total_weight) * 100

    for i in range(gmm.K):
        mean = gmm.mean[i][:2]
        cov = gmm.covar[i][:2, :2]
        eigvals, eigvecs = np.linalg.eigh(cov)
        order = eigvals.argsort()[::-1]
        eigvals, eigvecs = eigvals[order], eigvecs[:, order]
        angle = np.arctan2(*eigvecs[:, 0][::-1])
        width, height = 2 * np.sqrt(eigvals) * 2

        ell = Ellipse(xy=mean, width=width, height=height, angle=np.degrees(angle),
                      edgecolor=component_colors[i], facecolor='none', lw=2, label=f"Component {i+1}")
        ax.add_patch(ell)

    ax.set_xlabel(r"$v_R$ (km/s)")
    ax.set_ylabel(r"$v_\phi$ (km/s)")
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.grid(False)
    ax.text(x_limits[1] - 50, y_limits[0] + 20, f"{bin_label}\n{metallicity_range}",
            fontsize=12, color='black', ha='right', va='bottom')

    ax_frac = axes[0]
    bars = ax_frac.bar(range(1, gmm.K + 1), weights, color=component_colors[:gmm.K], alpha=0.8)
    for bar, weight in zip(bars, weights):
        ax_frac.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                     f"{weight:.4g}%", ha='center', fontsize=8)

    ax_frac.set_xticks(range(1, gmm.K + 1))
    ax_frac.set_xticklabels([f"Comp {i+1}" for i in range(gmm.K)])
    ax_frac.set_ylim(0, max(weights) * 1.2)

    plt.savefig(f'../figures/gmm_{label}.png', dpi=300, bbox_inches='tight')
    plt.show()


def extract_gmm_parameters(gmm, df, label, component_assignments=None):
    """
    Extracts and formats GMM parameters for structured reporting.

    Parameters
    ----------
    gmm : pygmmis.GMM
        Fitted GMM model.
    df : pd.DataFrame
        DataFrame the GMM was fitted on (used to count stars).
    label : str
        Label for the bin (e.g., a metallicity label).
    component_assignments : dict, optional
        Optional mapping from component names to GMM indices.

    Returns
    -------
    output : str
        Formatted string with means, standard deviations, and weights for each GMM component.
    """
    n_components = gmm.K
    num_stars = len(df)

    all_component_names = ["Stationary halo", "Prograde halo", "GS/E 1", "GS/E 2", "Thick Disc"]
    component_names = all_component_names[:n_components]

    means = gmm.mean.round(2)
    std_devs = np.sqrt(np.diagonal(gmm.covar, axis1=1, axis2=2)).round(2)
    weights = (gmm.amp / np.sum(gmm.amp) * 100).round(1)

    if component_assignments:
        ordered_indices = [component_assignments[name] for name in component_names if name in component_assignments]
    else:
        ordered_indices = np.argsort(means[:, 1])  # Sort by v_phi

    means = means[ordered_indices]
    std_devs = std_devs[ordered_indices]
    weights = weights[ordered_indices]

    df_results = pd.DataFrame({
        "Component": component_names,
        "Weights (%)": weights,
        r"$v_R$": means[:, 0],
        r"$\sigma_R$": std_devs[:, 0],
        r"$v_\phi$": means[:, 1],
        r"$\sigma_\phi$": std_devs[:, 1],
        r"$v_Z$": means[:, 2],
        r"$\sigma_Z$": std_devs[:, 2]
    })

    output = f"\n{label} ({num_stars} stars)\n"
    output += df_results.to_string(index=False, justify="center")
    return output


if __name__ == "__main__":
    print("This module provides GMM fitting, plotting, and parameter extraction utilities.")

