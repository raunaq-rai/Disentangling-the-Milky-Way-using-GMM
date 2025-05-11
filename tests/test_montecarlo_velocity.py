import sys
import os

# Ensure the src directory is in the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

import numpy as np
import pandas as pd
import astropy.units as u
import pytest
from astropy.coordinates import Galactocentric

from montecarlo_velocity import (
    generate_monte_carlo_samples,
    compute_velocity_components_with_uncertainty,
    compute_velocity_components_for_samples,
    process_data_monte_carlo
)

@pytest.fixture
def dummy_data():
    """Creates a small dummy dataframe for testing."""
    return pd.DataFrame({
        'ra': [10.0, 20.0],
        'dec': [30.0, -10.0],
        'rpgeo': [1000.0, 1500.0],  # pc
        'pmra': [5.0, -3.0],        # mas/yr
        'pmdec': [2.0, -2.5],       # mas/yr
        'radial_velocity': [20.0, -30.0],  # km/s
        'parallax_error': [0.01, 0.01],
        'pmra_error': [0.1, 0.1],
        'pmdec_error': [0.1, 0.1],
        'rpgeo_error': [10.0, 15.0],
        'radial_velocity_error': [1.0, 1.5]
    })

@pytest.fixture
def galactocentric_frame():
    """Returns a standard Galactocentric frame."""
    return Galactocentric()

def test_generate_monte_carlo_samples_shape(dummy_data):
    n_samples = 10
    output = generate_monte_carlo_samples(dummy_data, n_samples, correlation_pmra_pmdec=0.0)
    assert all(arr.shape == (2, n_samples) for arr in output)

def test_compute_velocity_components_with_uncertainty_shape(dummy_data, galactocentric_frame):
    df = dummy_data
    v_R, v_phi, v_Z = compute_velocity_components_with_uncertainty(
        df['ra'].values * u.deg,
        df['dec'].values * u.deg,
        df['rpgeo'].values * u.pc,
        df['pmra'].values * u.mas/u.yr,
        df['pmdec'].values * u.mas/u.yr,
        df['radial_velocity'].values * u.km/u.s,
        galactocentric_frame
    )
    assert len(v_R) == len(df)
    assert len(v_phi) == len(df)
    assert len(v_Z) == len(df)
    assert np.all(np.isfinite(v_R) | np.isnan(v_R))

def test_process_data_monte_carlo_output(dummy_data, galactocentric_frame):
    df_out = process_data_monte_carlo(dummy_data, galactocentric_frame, chunk_size=1, n_samples=5)
    assert 'v_R_uncertainty' in df_out.columns
    assert 'v_phi_uncertainty' in df_out.columns
    assert 'v_Z_uncertainty' in df_out.columns
    assert len(df_out) == len(dummy_data)

