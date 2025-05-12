import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.spatial import cKDTree
from tqdm import tqdm

# Thick disc parameters (3D Gaussian)
disc_mean_3d = np.array([0, 180, 0])  # Mean velocity of thick disc
disc_disp_3d = np.array([70, 50, 60])  # Velocity dispersion of thick disc
disc_2sigma = 2 * disc_disp_3d[:2]  # Ellipse size for v_R and v_phi

def inside_ellipse(x, y, center, width, height):
    """
    Check if points (x, y) lie within a 2D ellipse.

    Parameters:
    - x, y: Arrays of coordinates.
    - center: (x0, y0) center of the ellipse.
    - width, height: Half-widths of the ellipse (2σ).

    Returns:
    - Boolean array indicating points inside the ellipse.
    """
    return ((x - center[0]) / width) ** 2 + ((y - center[1]) / height) ** 2 <= 1

def assign_uncertainties(mock_stars, obs_stars, obs_errors):
    """
    Assigns uncertainties to mock stars using nearest-neighbor errors
    from the observed stars.

    Parameters:
    - mock_stars: Generated stars (N, 3)
    - obs_stars: Observed stars (N, 3)
    - obs_errors: Observed errors (N, 3)

    Returns:
    - Noisy mock stars with applied uncertainties
    """
    tree = cKDTree(obs_stars)
    _, idx = tree.query(mock_stars)
    assigned_errors = obs_errors[idx]
    noisy_mock = np.random.normal(mock_stars, assigned_errors)
    return noisy_mock

def generate_mock_with_errors(gmm, obs_stars, obs_errors):
    """
    Generate mock stars from GMM and perturb with observed errors.

    Returns:
    - noisy_mock: Mock stars with errors applied
    """
    mock_stars = gmm.draw(len(obs_stars))[:, :3]
    return assign_uncertainties(mock_stars, obs_stars, obs_errors)

def compute_residual_map(df_bin, gmm, bins=100, vR_lim=(-400, 400), vphi_lim=(-400, 400)):
    """
    Computes normalized residual map between observed and mock star distributions.

    Returns:
    - H_residual: Residual map
    - xedges, yedges: Bin edges for plotting
    """
    obs_stars = df_bin[['v_R', 'v_phi', 'v_Z']].values
    obs_errors = df_bin[['v_R_uncertainty', 'v_phi_uncertainty', 'v_Z_uncertainty']].values
    mock_stars = generate_mock_with_errors(gmm, obs_stars, obs_errors)

    H_obs, xedges, yedges = np.histogram2d(obs_stars[:, 0], obs_stars[:, 1], bins=bins, range=[vR_lim, vphi_lim])
    H_mock, _, _ = np.histogram2d(mock_stars[:, 0], mock_stars[:, 1], bins=bins, range=[vR_lim, vphi_lim])

    H_residual = (H_obs - H_mock) / (H_obs + H_mock + 1e-5)
    return H_residual, xedges, yedges

def run_residual_analysis(df_bin, gmm, disc_fractions, n_realizations=200):
    """
    Runs a Monte Carlo analysis to test detectability of a disc component.

    Parameters:
    - df_bin: DataFrame of the velocity data
    - gmm: Fitted GMM object
    - disc_fractions: List of fractions to inject thick disc stars
    - n_realizations: Number of MC realizations

    Returns:
    - obs_mean, obs_std: Mean and std of residuals without injection
    - frac_means, frac_stds: Residuals with injected disc fractions
    """
    obs_stars = df_bin[['v_R', 'v_phi', 'v_Z']].values
    obs_errors = df_bin[['v_R_uncertainty', 'v_phi_uncertainty', 'v_Z_uncertainty']].values

    # Observed residual
    obs_residuals = []
    for _ in tqdm(range(n_realizations), desc="Observed residual MC"):
        baseline_mock = gmm.draw(len(obs_stars))[:, :3]
        noisy_mock = assign_uncertainties(baseline_mock, obs_stars, obs_errors)
        obs_in_disc = np.sum(inside_ellipse(obs_stars[:, 0], obs_stars[:, 1], disc_mean_3d[:2], *disc_2sigma))
        mock_in_disc = np.sum(inside_ellipse(noisy_mock[:, 0], noisy_mock[:, 1], disc_mean_3d[:2], *disc_2sigma))
        obs_residuals.append(obs_in_disc - mock_in_disc)

    obs_mean = np.mean(obs_residuals)
    obs_std = np.std(obs_residuals)

    # Injected residuals
    frac_means, frac_stds = [], []
    for frac in tqdm(disc_fractions, desc="Injected disc fractions"):
        mock_disc_residuals = []
        for _ in range(n_realizations):
            baseline_mock = gmm.draw(len(obs_stars))[:, :3]
            noisy_baseline = assign_uncertainties(baseline_mock, obs_stars, obs_errors)
            mock_in_disc = np.sum(inside_ellipse(noisy_baseline[:, 0], noisy_baseline[:, 1], disc_mean_3d[:2], *disc_2sigma))

            injected_mock = gmm.draw(len(obs_stars))[:, :3]
            n_inject = int(frac * len(obs_stars))
            if n_inject > 0:
                disc_stars = np.random.normal(disc_mean_3d, disc_disp_3d, size=(n_inject, 3))
                injected_mock = np.vstack([injected_mock, disc_stars])
            noisy_injected = assign_uncertainties(injected_mock, obs_stars, obs_errors)
            injected_in_disc = np.sum(inside_ellipse(noisy_injected[:, 0], noisy_injected[:, 1], disc_mean_3d[:2], *disc_2sigma))

            mock_disc_residuals.append(injected_in_disc - mock_in_disc)

        frac_means.append(np.mean(mock_disc_residuals))
        frac_stds.append(np.std(mock_disc_residuals))

    return obs_mean, obs_std, frac_means, frac_stds

def plot_residual_map(H_residual, xedges, yedges, name):
    """
    Plot the residual map with an overlaid disc ellipse.

    Parameters:
    - H_residual: Residual matrix
    - xedges, yedges: Axis bins
    - name: Title for the plot
    """
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.pcolormesh(xedges, yedges, H_residual.T, cmap='coolwarm', shading='auto', vmin=-1, vmax=1)
    ellipse = Ellipse(xy=(0, 180), width=2*70*2, height=2*50*2,
                      edgecolor='black', facecolor='none', linestyle='--', linewidth=2)
    ax.add_patch(ellipse)
    ax.set_xlabel('v_R (km/s)')
    ax.set_ylabel('v_phi (km/s)')
    ax.set_title(f"{name} Normalized Residual Map")
    ax.grid(True)
    fig.colorbar(im, ax=ax, orientation='vertical', label='Normalized Residual')
    plt.tight_layout()
    plt.show()

def plot_disc_fraction(disc_fractions, results, name):
    """
    Plot the residual vs. injected disc fraction curve.

    Parameters:
    - disc_fractions: List of fractions
    - results: Output from `run_residual_analysis`
    - name: Plot title
    """
    obs_mean, obs_std, frac_means, frac_stds = results
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.axhline(obs_mean, color='blue', linestyle='--', label='Observation result')
    ax.fill_between(disc_fractions, obs_mean - obs_std, obs_mean + obs_std, color='blue', alpha=0.3)
    ax.plot(disc_fractions, frac_means, color='black', marker='x', label='Mock star test')
    ax.fill_between(disc_fractions,
                    np.array(frac_means) - np.array(frac_stds),
                    np.array(frac_means) + np.array(frac_stds),
                    color='red', alpha=0.3)
    ax.set_xlabel('Disc fraction')
    ax.set_ylabel('Disc star residual')
    ax.set_ylim(-50, 250)
    ax.set_title(f"{name} Disc Residual vs. Disc Fraction")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()

