# bic_vs_components.py

import numpy as np
import pygmmis
import matplotlib.pyplot as plt
from tqdm import tqdm

def compute_bic_vs_n_components(df_bin, max_components=8, n_init=50):
    """
    Compute BIC for different numbers of Gaussian components (1 to max_components),
    using Extreme Deconvolution (XD) to account for measurement uncertainties.

    Args:
    - df_bin (pd.DataFrame): DataFrame containing velocity and uncertainty columns.
    - max_components (int): Maximum number of GMM components to evaluate.
    - n_init (int): Number of initializations per component count.

    Returns:
    - BIC_values (dict): Mapping from component count to list of BIC values.
    """
    X = df_bin[['v_R', 'v_phi', 'v_Z']].values
    n = len(X)

    # Diagonal covariance matrices per star
    uncertainties = np.stack([
        df_bin['v_R_uncertainty'].values**2,
        df_bin['v_phi_uncertainty'].values**2,
        df_bin['v_Z_uncertainty'].values**2
    ], axis=1)
    cov_matrices = np.array([np.diag(u) for u in uncertainties])

    BIC_values = {N: [] for N in range(1, max_components + 1)}

    print("\nComputing BIC for different numbers of components...")
    rng = np.random.RandomState()

    for N in tqdm(range(1, max_components + 1), desc="Fitting GMMs"):
        for _ in range(n_init):
            gmm = pygmmis.GMM(K=N, D=3)
            pygmmis.initFromKMeans(gmm, X)
            logL, _ = pygmmis.fit(gmm, X, covar=cov_matrices, w=0.1, tol=1e-6, init_method='kmeans')

            # BIC = k*ln(n) - 2*logL
            k = (1 + 3 + 6) * N - 1  # weights, means, covariances
            BIC = k * np.log(n) - 2 * logL
            BIC_values[N].append(BIC)

    return BIC_values


def plot_bic_vs_n_components(BIC_values, fig_name):
    """
    Plots BIC vs. Number of Components and highlights the minimum BIC point.

    Args:
    - BIC_values: Dictionary with BIC values for each component count.
    - fig_name: Path/filename to save the plot (e.g. 'output.png').
    """
    
    num_components = sorted([int(k) for k in BIC_values.keys()])
    
    smallest_bic = [min(BIC_values[str(N)] if str(N) in BIC_values else BIC_values[N]) for N in num_components]
    median_bic   = [np.median(BIC_values[str(N)] if str(N) in BIC_values else BIC_values[N]) for N in num_components]
    largest_bic  = [max(BIC_values[str(N)] if str(N) in BIC_values else BIC_values[N]) for N in num_components]
    q25_bic      = [np.percentile(BIC_values[str(N)] if str(N) in BIC_values else BIC_values[N], 25) for N in num_components]
    q75_bic      = [np.percentile(BIC_values[str(N)] if str(N) in BIC_values else BIC_values[N], 75) for N in num_components]

    # Find the minimum BIC value and corresponding component count
    min_bic = min(smallest_bic)
    min_index = smallest_bic.index(min_bic)
    best_n_components = num_components[min_index]

    plt.figure(figsize=(8, 6))
    plt.plot(num_components, smallest_bic, 'k-', label='Smallest BIC')
    plt.plot(num_components, q25_bic, 'g-.', label='25th Percentile BIC')
    plt.plot(num_components, median_bic, 'b--', label='Median BIC')
    plt.plot(num_components, q75_bic, 'm-.', label='75th Percentile BIC')
    plt.plot(num_components, largest_bic, 'r:', label='Largest BIC')

    # Highlight the minimum BIC value
    plt.plot(best_n_components, min_bic, 'ro', markersize=8, label='Minimum BIC')

    plt.xlabel("Number of GMM Components", fontsize=14)
    plt.ylabel("BIC Value", fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(fig_name, dpi=300)
    plt.show()
