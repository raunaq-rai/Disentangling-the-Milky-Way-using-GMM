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

