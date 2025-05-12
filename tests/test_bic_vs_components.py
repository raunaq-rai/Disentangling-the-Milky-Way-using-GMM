import sys
import os

# Ensure the src directory is in the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

# tests/test_bic_vs_components.py

import numpy as np
import pandas as pd
from bic_vs_components import compute_bic_vs_n_components

def test_compute_bic_output_structure():
    # Create a small synthetic DataFrame with mock values
    n_samples = 20
    data = {
        'v_R': np.random.normal(0, 50, n_samples),
        'v_phi': np.random.normal(200, 30, n_samples),
        'v_Z': np.random.normal(0, 20, n_samples),
        'v_R_uncertainty': np.full(n_samples, 5.0),
        'v_phi_uncertainty': np.full(n_samples, 5.0),
        'v_Z_uncertainty': np.full(n_samples, 5.0),
    }
    df_mock = pd.DataFrame(data)

    # Run the function
    bic_results = compute_bic_vs_n_components(df_mock, max_components=2, n_init=3)

    # Check output is a dictionary with correct keys and values
    assert isinstance(bic_results, dict)
    assert all(isinstance(k, int) for k in bic_results.keys())
    assert set(bic_results.keys()) == {1, 2}
    assert all(isinstance(v, list) for v in bic_results.values())
    assert all(len(v) == 3 for v in bic_results.values())  # n_init = 3

    print("[✓] test_compute_bic_output_structure passed.")

