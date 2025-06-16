# Disentangling the Components of the Milky Way with Gaussian Mixture Modelling

## Overview
This project aims to study the formation history of the Milky Way by analysing its stellar components, with a particular focus on old, low-metallicity stars. The goal is to disentangle the components of the Galaxy (e.g., thin disc, thick disc, and halo) using data from the Gaia DR3 release and apply Gaussian Mixture Modelling (GMM) to identify structures in velocity space.

## Repository
This repo includes the workflow used and notebook used in this project - note that it is a working progress!!!

## Project Goals

The primary goal of this project is to leverage recent advancements in astrometric and spectroscopic data provided by the Gaia satellite's Data Release 3 (DR3) to disentangle the components of the Milky Way using Gaussian Mixture Modelling (GMM). Specifically, we aim to:

1. **Identify and characterise the components of the Milky Way** in different metallicity regimes, focusing on the oldest and most metal-poor stars.
2. **Test for the presence of a disc component at very low metallicities** ([M/H] < -2) and place quantitative limits on its contribution if present.
3. **Refine the methodology** for analysing stellar kinematics and chemical properties, paving the way for future investigations into Galactic formation scenarios.

This project introduces novel analysis techniques to address open questions about the formation of the earliest structures in the Milky Way, leveraging state-of-the-art data and statistical methods.

## Why This Matters

### Galactic Formation History
The Milky Way's structure is composed of several key components, including the young thin disc, the older thick disc, and the ancient spherical halo. By studying the properties and dynamics of the stars within these components, we can infer their formation history and evolutionary processes. Of particular interest are the stars formed in the early Universe, which are identifiable by their extremely low metallicities. 

### Early Disc Formation
While the halo is commonly associated with the oldest stars, recent studies have hinted at the existence of a disc component even at very low metallicities. The confirmation of such a disc would challenge traditional views of Galactic formation, suggesting that the Milky Way's disc structure began to form earlier than previously thought. Understanding the properties of this ancient disc, if it exists, would provide crucial insights into the assembly of the Milky Way and the broader context of galaxy formation in the Universe.

## Resources

- **Data:**
  - Gaia DR3 vetted RGB sample with metalicity ([Andrae et al. 2023](https://zenodo.org/records/7945154)).
  - Distances from [Bailer-Jones et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021AJ....161..147B/abstract).
  - Alpha Abundance samples [Li et al. (2024)](https://arxiv.org/abs/2309.14294)
- **Key publications:**
  - [Zhang, Ardern-Arentsen & Belokurov (2024)](https://arxiv.org/pdf/2311.09294): "On the existence of a very metal-poor disc in the Milky Way".
  - [Belokurov & Deason (2024)](https://arxiv.org/pdf/2402.12443): "Galactic Archaeology with Gaia".
  - [Bovy, Hogg & Roweis (2011)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-2B/Extreme-deconvolution--Inferring-complete-distribution-functions-from-noisy-heterogeneous/10.1214/10-AOAS439.full): "Extreme Deconvolution: Inferring complete distribution functions from noisy, heterogeneous and incomplete observations".

## Interactive 3D GMM Plots with Weights

Click any of the links below to open the full interactive 3D Plotly widget (with component‐weight bar on the right), served via RawGitHack:

- [VMP : −3 < [M/H] < −2](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/VMP__-3%5BM_H%5D-2.html)  
- [IMP : −2 < [M/H] < −1.6](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/IMP__-2%5BM_H%5D-1.6.html)  
- [MP1 : −1.6 < [M/H] < −1.3](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/MP1__-1.6%5BM_H%5D-1.3.html)  
- [MP2 : −1.3 < [M/H] < −1.0](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/MP2__-1.3%5BM_H%5D-1.0.html)  
