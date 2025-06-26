.. Disentagling the Components of the Milky Way documentation master file, created by
   sphinx-quickstart on Thu Jun 26 15:27:43 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Disentangling the Components of the Milky Way
=============================================

Overview
--------

This project investigates the formation history of the Milky Way by analysing its
stellar components—specifically old, low-metallicity stars.  Using Gaia DR3 data
we apply Gaussian Mixture Modelling (GMM) to disentangle Galactic structures in
velocity space.  We first replicate the analysis of Zhang *et al.* (2024) and
extend it by including alpha-element abundances, following Viswanathan
*et al.* (2024).

Goals
-----

- **Identify and characterise** Milky Way components across metallicity regimes  
- **Test** for a disc component at very low metallicities (:math:`[\mathrm{M/H}]<-2`)  
- **Examine** how chemical composition correlates with Galactic kinematic structure  

Scientific Motivation
---------------------

The Milky Way comprises several structural components including the thin and
thick discs and the stellar halo. While the halo is usually regarded as the home
of the oldest stars, recent studies hint that a disc may exist even at very low
metallicities. Confirming this would challenge traditional pictures of Galactic
evolution and shed new light on the Galaxy’s earliest assembly.

Resources
---------

*Data*

- Gaia DR3 RGB sample with metallicity — `Andrae et al. 2023 <https://zenodo.org/records/7945154>`_
- Distances — `Bailer-Jones et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021AJ....161..147B/abstract>`_
- Alpha abundances — `Li et al. 2024 <https://arxiv.org/abs/2309.14294>`_

*Key publications*

- Zhang et al. 2024 — `On the existence of a very metal-poor disc <https://arxiv.org/pdf/2311.09294>`_
- Belokurov & Deason 2024 — `Galactic Archaeology with Gaia <https://arxiv.org/pdf/2402.12443>`_
- Bovy et al. 2011 — `Extreme Deconvolution <https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-2B/Extreme-deconvolution--Inferring-complete-distribution-functions-from-noisy-heterogeneous/10.1214/10-AOAS439.full>`_

Methodology
-----------

1. **Gaussian Mixture Modelling**  
   The 3-D velocity distribution :math:`(v_R,\,v_\phi,\,v_Z)` is modelled as a
   sum of multivariate Gaussians, each representing a kinematic component
   (disc, halo, GS/E, etc.).

2. **Extreme Deconvolution**  
   We use *pyGMMis* to perform Extreme Deconvolution, incorporating the
   per-star velocity uncertainties provided by Gaia directly in the fit.

3. **Metallicity binning & component selection**  
   Stars are split into metallicity bins (:math:`[\mathrm{M/H}]<-1`).  For each
   bin we run 50 random initialisations and choose the number of Gaussians that
   minimises the Bayesian Information Criterion (BIC).

Interactive Visualisation
-------------------------

Interactive Plotly figures are listed in the project *README* and main report.


What's in this documentation?
-----------------------------

This documentation provides:

- **API documentation** for all modules in the `src/` directory, including:
  - Data loading and preprocessing
  - Model fitting and evaluation
  - Plotting and visualisation functions

Each module includes:
- Function and class docstrings
- Argument descriptions
- Return values


Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Project Modules:

   modules
   src

