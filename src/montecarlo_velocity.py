"""
montecarlo_velocity.py

This module provides utilities for computing velocity uncertainties of stars using
Monte Carlo sampling, based on astrometric measurements and the Galactocentric frame.
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
import astropy.units as u
from astropy.coordinates import ICRS, Galactocentric, CylindricalRepresentation, CylindricalDifferential


def generate_monte_carlo_samples(df_chunk, n_samples, correlation_pmra_pmdec):
    """
    Generate Monte Carlo samples for each star using its uncertainties.

    Parameters:
    - df_chunk (pd.DataFrame): Subset of the dataframe containing stellar data.
    - n_samples (int): Number of samples per star.
    - correlation_pmra_pmdec (float): Correlation coefficient between PMRA and PMDEC.

    Returns:
    - Tuple of arrays: (ra_samples, dec_samples, distance_samples,
                        pmra_samples, pmdec_samples, vlos_samples)
    """
    num_stars = len(df_chunk)

    ra, dec, distance = df_chunk['ra'].values, df_chunk['dec'].values, df_chunk['rpgeo'].values
    pmra, pmdec = df_chunk['pmra'].values, df_chunk['pmdec'].values
    vlos = df_chunk['radial_velocity'].values

    pmra_err = df_chunk['pmra_error'].values
    pmdec_err = df_chunk['pmdec_error'].values
    dist_err = df_chunk['rpgeo_error'].values
    vlos_err = df_chunk['radial_velocity_error'].values
    parallax_err = df_chunk['parallax_error'].values

    ra_samples = np.random.normal(ra[:, None], parallax_err[:, None], (num_stars, n_samples))
    dec_samples = np.random.normal(dec[:, None], parallax_err[:, None], (num_stars, n_samples))
    distance_samples = np.random.normal(distance[:, None], dist_err[:, None], (num_stars, n_samples))
    vlos_samples = np.random.normal(vlos[:, None], vlos_err[:, None], (num_stars, n_samples))

    pmra_samples = np.zeros((num_stars, n_samples))
    pmdec_samples = np.zeros((num_stars, n_samples))

    for i in tqdm(range(num_stars), desc="Generating Monte Carlo Samples"):
        cov = correlation_pmra_pmdec * pmra_err[i] * pmdec_err[i]
        cov_matrix = np.array([[pmra_err[i]**2, cov], [cov, pmdec_err[i]**2]])
        samples = np.random.multivariate_normal([pmra[i], pmdec[i]], cov_matrix, n_samples)
        pmra_samples[i, :], pmdec_samples[i, :] = samples[:, 0], samples[:, 1]

    dec_samples = np.clip(dec_samples, -90, 90)
    ra_samples = np.mod(ra_samples, 360)
    distance_samples = np.clip(distance_samples, 1e-5, None)

    return ra_samples, dec_samples, distance_samples, pmra_samples, pmdec_samples, vlos_samples


def compute_velocity_components_with_uncertainty(ra, dec, distance, pmra, pmdec, vlos, gc_frame):
    """
    Compute cylindrical velocity components (v_R, v_phi, v_Z) in km/s for stars, 
    transforming input astrometric data to a Galactocentric frame.

    Parameters
    ----------
    ra, dec : astropy.units.Quantity
        Sky positions (e.g., in degrees).
    distance : astropy.units.Quantity
        Distances (e.g., in parsecs).
    pmra, pmdec : astropy.units.Quantity
        Proper motions (e.g., in mas/yr).
    vlos : astropy.units.Quantity
        Radial velocities (e.g., in km/s).
    gc_frame : Galactocentric
        Target Galactocentric frame.

    Returns
    -------
    v_R, v_phi, v_Z : np.ndarray
        Cylindrical velocity components in km/s, with NaNs for invalid inputs.
    """
    # Identify valid rows
    mask = (
        np.isfinite(ra) & np.isfinite(dec) & np.isfinite(distance) &
        np.isfinite(pmra) & np.isfinite(pmdec) & np.isfinite(vlos)
    )

    # Preallocate output arrays
    v_R_full = np.full_like(ra.value, np.nan)
    v_phi_full = np.full_like(ra.value, np.nan)
    v_Z_full = np.full_like(ra.value, np.nan)

    if not np.any(mask):
        return v_R_full, v_phi_full, v_Z_full

    coords = ICRS(
        ra=ra[mask], dec=dec[mask], distance=distance[mask],
        pm_ra_cosdec=pmra[mask], pm_dec=pmdec[mask], radial_velocity=vlos[mask]
    )

    cg = coords.transform_to(gc_frame)

    try:
        cg_cyl = cg.represent_as(CylindricalRepresentation, CylindricalDifferential)
        v_R = cg_cyl.differentials['s'].d_rho.to(u.km/u.s).value
        v_phi = -(cg_cyl.differentials['s'].d_phi.to(u.rad/u.s) * cg_cyl.rho.to(u.km)).value
        v_Z = cg_cyl.differentials['s'].d_z.to(u.km/u.s).value

        v_R_full[mask] = v_R
        v_phi_full[mask] = v_phi
        v_Z_full[mask] = v_Z

    except (AttributeError, TypeError) as e:
        # Catch cases where differentials are missing entirely (e.g. empty cg)
        pass

    return v_R_full, v_phi_full, v_Z_full


def compute_velocity_components_for_samples(ra_samples, dec_samples, distance_samples,
                                            pmra_samples, pmdec_samples, vlos_samples, gc_frame):
    """
    Compute (v_R, v_phi, v_Z) for all Monte Carlo samples.

    Returns:
    - v_R_samples, v_phi_samples, v_Z_samples: np.ndarray
    """
    num_stars, n_samples = ra_samples.shape

    v_R_samples = np.zeros((num_stars, n_samples))
    v_phi_samples = np.zeros((num_stars, n_samples))
    v_Z_samples = np.zeros((num_stars, n_samples))

    for i in tqdm(range(n_samples), desc="Computing Velocities for Monte Carlo Samples"):
        v_R, v_phi, v_Z = compute_velocity_components_with_uncertainty(
            ra_samples[:, i] * u.deg, dec_samples[:, i] * u.deg,
            distance_samples[:, i] * u.pc,
            pmra_samples[:, i] * u.mas/u.yr, pmdec_samples[:, i] * u.mas/u.yr,
            vlos_samples[:, i] * u.km/u.s, gc_frame
        )
        v_R_samples[:, i], v_phi_samples[:, i], v_Z_samples[:, i] = v_R, v_phi, v_Z

    return v_R_samples, v_phi_samples, v_Z_samples


def process_data_monte_carlo(df, gc_frame, chunk_size=100000, n_samples=100, correlation_pmra_pmdec=0.0):
    """
    Process dataset in chunks, compute Monte Carlo velocity uncertainties.

    Returns:
    - pd.DataFrame: With added v_R_uncertainty, v_phi_uncertainty, v_Z_uncertainty
    """
    num_chunks = len(df) // chunk_size + 1
    df_combined = []

    for chunk_num in range(num_chunks):
        start_idx, end_idx = chunk_num * chunk_size, min((chunk_num + 1) * chunk_size, len(df))
        df_chunk = df.iloc[start_idx:end_idx].reset_index(drop=True)

        ra_samp, dec_samp, dist_samp, pmra_samp, pmdec_samp, vlos_samp = generate_monte_carlo_samples(
            df_chunk, n_samples, correlation_pmra_pmdec
        )

        v_R_samp, v_phi_samp, v_Z_samp = compute_velocity_components_for_samples(
            ra_samp, dec_samp, dist_samp, pmra_samp, pmdec_samp, vlos_samp, gc_frame
        )

        df_chunk['v_R_uncertainty'] = np.std(v_R_samp, axis=1)
        df_chunk['v_phi_uncertainty'] = np.std(v_phi_samp, axis=1)
        df_chunk['v_Z_uncertainty'] = np.std(v_Z_samp, axis=1)

        df_combined.append(df_chunk)
        print(f"Processed chunk {chunk_num + 1}/{num_chunks}")

    return pd.concat(df_combined, ignore_index=True)
