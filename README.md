# Disentangling the Components of the Milky Way with Gaussian Mixture Modelling

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)


## Overview
This project aims to study the formation history of the Milky Way by analysing its stellar components, with a particular focus on old, low-metallicity stars. The goal is to disentangle the components of the Galaxy (e.g., thin disc, thick disc, and halo) using data from the Gaia DR3 release and apply Gaussian Mixture Modelling (GMM) to identify structures in velocity space. This report replicates and extends the analysis done by [Zhang et al. (2024)](https://arxiv.org/pdf/2311.09294).

We then extend this analysis by incorporating alpha-element abundance to analysis as done by [Viswanathan et al . (2024)](https://arxiv.org/abs/2411.12165). We fit a number of gaussian components for each bin governed by BIC score and observe the difference between Milky Way components in the high and low alpha regimes.


## Table of Contents
- [Overview](#overview)
- [Project&nbsp;Goals](#project-goals)
- [Why&nbsp;This&nbsp;Matters](#why-this-matters)
- [Resources](#resources)
- [Usage](#usage)
- [Methodology](#methodology)
- [Documentation](#documentation)
- [License](#license)
- [Replication&nbsp;Results](#replication-results-interactive-3d-gmm-plots-with-weights)
- [Expansion&nbsp;Results](#expansion-results-interactive-3d-gmm-plots-with-weights)
- [Authors and Acknowledgment](#authors-and-acknowledgment)

## Project Goals

The primary goal of this project is to leverage recent advancements in astrometric and spectroscopic data provided by the Gaia satellite's Data Release 3 (DR3) to disentangle the components of the Milky Way using Gaussian Mixture Modelling (GMM). Specifically, we aim to:

1. **Identify and characterise the components of the Milky Way** in different metallicity regimes, focusing on the oldest and most metal-poor stars.
2. **Test for the presence of a disc component at very low metallicities** ([M/H] < -2) and place quantitative limits on its contribution if present.
3. **Expand on analysis** for further understanding the relationship between chemical composition and milkyway components.

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

## Usage

create environment:

```bash
git clone 'https://gitlab.developers.cam.ac.uk/phy/data-intensive-science-mphil/assessments/projects/rsr45.git'
cd rsr45
conda env create -f environment.yaml
conda activate research_project_MW_GMM
python -m ipykernel install --user --name research_project_MW_GMM --display-name "research_project_MW_GMM"
```

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).


## Methodology

1. **Gaussian Mixture Modeling**  
   - We model the 3D velocity distribution \((v_R,\,v_\phi,\,v_Z)\) of stars as a sum of multivariate Gaussians.  
   - Each Gaussian represents a kinematic sub-population (e.g.\ halo, disc, GS/E) and yields a probabilistic (“soft”) assignment of stars to components.  

2. **Extreme Deconvolution with Uncertainties**  
   - Fitting is done via the pyGMMis package, which extends EM with “Extreme Deconvolution” to incorporate per-star velocity uncertainties (from Gaia) directly into the likelihood.  

3. **Metallicity Bins & Component Selection**  
   - We split stars into four metallicity bins ([M/H] < –1) and fit a separate GMM to each.  
   - The optimal number of Gaussians is chosen by minimizing the Bayesian Information Criterion (BIC) over multiple (50) random initializations, with k-means seeding to ensure stable convergence.  

## Documentation

- Documentation was created using [Sphinx](https://www.sphinx-doc.org/).
- To build it locally:

  ```bash
  cd docs
  make html
  open build/html/index.html 
  ```

- This documentation provides API documentation for all modules in the src directory.

## Replication Results: Interactive 3D GMM Plots with Weights

2D Visualisations are shown in the report. Below we use plotly to include an interactive visualisation for the 3-D velocity distribution. Please click on each link!!!

Click any of the links below to open the full interactive 3D Plotly widget (with component‐weight bar on the right), served via RawGitHack:

- [VMP : −3 < [M/H] < −2](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/VMP__-3%5BM_H%5D-2.html)  
- [IMP : −2 < [M/H] < −1.6](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/IMP__-2%5BM_H%5D-1.6.html)  
- [MP1 : −1.6 < [M/H] < −1.3](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/MP1__-1.6%5BM_H%5D-1.3.html)  
- [MP2 : −1.3 < [M/H] < −1.0](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/MP2__-1.3%5BM_H%5D-1.0.html)  

## Expansion Results: Interactive 3D GMM Plots with Weights

2D Visualisations are shown in the report. Below we use plotly to include an interactive visualisation for the 3-D velocity distribution. Please click on each link!!!

**High Alpha**

- **High–α VMP** (`-3 < [M/H] < -2`):  
  [High–α VMP](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/VMP_high___-3[M_H]-2.html)

- **High–α IMP** (`-2 < [M/H] < -1.6`):  
  [High–α IMP](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/IMP_high___-2[M_H]-1.6.html)

- **High–α MP1** (`-1.6 < [M/H] < -1.3`):  
  [High–α MP1](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/MP1_high___-1.6[M_H]-1.3.html)

- **High–α MP2** (`-1.3 < [M/H] < -1.0`):  
  [High–α MP2](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/MP2_high___-1.3[M_H]-1.0.html)


**Low Alpha**

- **Low–α VMP** (`-3 < [M/H] < -2`):  
  [Low–α VMP](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/VMP_low____-3[M_H]-2.html)

- **Low–α IMP** (`-2 < [M/H] < -1.6`):  
  [Low–α IMP](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/IMP_low____-2[M_H]-1.6.html)

- **Low–α MP1** (`-1.6 < [M/H] < -1.3`):  
  [Low–α MP1](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/MP1_low____-1.6[M_H]-1.3.html)

- **Low–α MP2** (`-1.3 < [M/H] < -1.0`):  
  [Low–α MP2](https://raw.githack.com/raunaq-rai/Disentangling-the-Milky-Way-using-GMM/main/figures/MP2_low____-1.3[M_H]-1.0.html)


## Authors and Acknowledgment
This project is maintained by [Raunaq Rai](https://www.linkedin.com/in/raunaq-rai/) at the University of Cambridge.

Special thanks to my supervisor, [Dr. Anke Arentsen](https://www.ast.cam.ac.uk/people/anke.ardern-arentsen), for her guidance, expertise, and encouragement throughout this project.


30th June 2025
