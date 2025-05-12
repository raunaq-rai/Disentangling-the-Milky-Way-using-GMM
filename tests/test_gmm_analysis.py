# test_gmm_analysis.py

import sys
import os

# Ensure the src directory is in the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt
from pygmmis import GMM
from gmm_analysis import fit_gmm_fixed_components, extract_gmm_parameters, plot_gmm_with_contributions

@pytest.fixture
def mock_df():
    """Generate a small DataFrame with mock velocity data and uncertainties."""
    np.random.seed(42)
    size = 100
    return pd.DataFrame({
        'v_R': np.random.normal(0, 50, size),
        'v_phi': np.random.normal(150, 30, size),
        'v_Z': np.random.normal(0, 20, size),
        'v_R_uncertainty': np.random.uniform(2, 5, size),
        'v_phi_uncertainty': np.random.uniform(2, 5, size),
        'v_Z_uncertainty': np.random.uniform(2, 5, size),
    })

def test_fit_gmm_fixed_components(mock_df):
    """Test that GMM fitting returns a valid GMM object with correct number of components."""
    gmm = fit_gmm_fixed_components(mock_df, n_components=3, n_init=5)
    assert isinstance(gmm, GMM)
    assert gmm.K == 3
    assert gmm.mean.shape == (3, 3)
    assert gmm.covar.shape == (3, 3, 3)

def test_extract_gmm_parameters(mock_df):
    """Test that GMM parameters can be extracted into a string."""
    gmm = fit_gmm_fixed_components(mock_df, n_components=3, n_init=2)
    output = extract_gmm_parameters(gmm, mock_df, label="[M/H] ~ -1.0")
    assert isinstance(output, str)
    assert "[M/H] ~ -1.0" in output
    assert "Component" in output

def test_plot_gmm_with_contributions_runs(mock_df, monkeypatch):
    """Test that plotting function runs without error (without saving file)."""

    # Prevent file saving during test
    monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: None)

    gmm = fit_gmm_fixed_components(mock_df, n_components=2, n_init=2)
    try:
        plot_gmm_with_contributions(
            mock_df, gmm,
            bin_label="[M/H] Bin",
            label="test_plot",
            metallicity_range="[-1.0, -0.8]"
        )
    except Exception as e:
        pytest.fail(f"Plotting failed with error: {e}")
