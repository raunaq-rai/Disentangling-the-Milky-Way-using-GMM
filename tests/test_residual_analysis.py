import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

import pytest
import numpy as np
import pandas as pd
from pygmmis import GMM
from residual_analysis import (
    inside_ellipse,
    assign_uncertainties,
    compute_residual_map,
    run_residual_analysis,
    plot_residual_map,
    plot_disc_fraction
)

@pytest.fixture
def mock_df():
    np.random.seed(42)
    n = 200
    df = pd.DataFrame({
        'v_R': np.random.normal(0, 100, size=n),
        'v_phi': np.random.normal(180, 60, size=n),
        'v_Z': np.random.normal(0, 80, size=n),
        'v_R_uncertainty': np.random.uniform(2, 5, size=n),
        'v_phi_uncertainty': np.random.uniform(2, 5, size=n),
        'v_Z_uncertainty': np.random.uniform(2, 5, size=n)
    })
    return df

@pytest.fixture
def simple_gmm():
    gmm = GMM(K=2, D=3)
    gmm.mean = np.array([[0, 180, 0], [50, 150, 0]])
    gmm.covar = np.array([np.eye(3) * 100 for _ in range(2)])
    gmm.amp = np.array([0.5, 0.5])
    return gmm

def test_inside_ellipse():
    assert inside_ellipse(0, 180, [0, 180], 140, 100) is True
    assert inside_ellipse(200, 180, [0, 180], 140, 100) is False

def test_assign_uncertainties_runs(mock_df):
    obs = mock_df[['v_R', 'v_phi', 'v_Z']].values
    err = mock_df[['v_R_uncertainty', 'v_phi_uncertainty', 'v_Z_uncertainty']].values
    noisy = assign_uncertainties(obs, obs, err)
    assert noisy.shape == obs.shape

def test_compute_residual_map_shape(mock_df, simple_gmm):
    residual, xedges, yedges = compute_residual_map(mock_df, simple_gmm, bins=50)
    assert residual.shape == (50, 50)

def test_run_residual_analysis_runs(mock_df, simple_gmm):
    disc_fractions = [0.0, 0.2, 0.5]
    result = run_residual_analysis(mock_df, simple_gmm, disc_fractions, n_realizations=5)
    assert len(result) == 4
    assert len(result[2]) == len(disc_fractions)

def test_plotting_functions(mock_df, simple_gmm):
    # Just run them to check they don't error
    residual, xedges, yedges = compute_residual_map(mock_df, simple_gmm)
    plot_residual_map(residual, xedges, yedges, name="Test")

    results = run_residual_analysis(mock_df, simple_gmm, [0.0, 0.3], n_realizations=5)
    plot_disc_fraction([0.0, 0.3], results, name="Test")
