# EHT Sgr A* Polarization Imaging Research
## Raw Research Notes
**Date:** 2026-03-01

---

## Sources

### EHT Data Releases
1. [EHT Data Products Page](https://eventhorizontelescope.org/for-astronomers/data)
2. [GitHub: 2024-D02-01 - EHT Sgr A* Polarized Data](https://github.com/eventhorizontelescope/2024-D02-01)
3. [CyVerse: EHTC FirstSgrAPol Mar2024](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstSgrAPol_Mar2024)
4. [GitHub: 2022-D02-01 - First Sgr A* EHT Results: Calibrated Data (Stokes I only)](https://github.com/eventhorizontelescope/2022-D02-01)

### EHT Sgr A* Polarization Papers
5. [EHT Paper VII: Polarization of the Ring (ApJL 964 L25)](https://iopscience.iop.org/article/10.3847/2041-8213/ad2df0)
6. [EHT Paper VIII: Physical Interpretation of the Polarized Ring (ApJL 964 L26)](https://iopscience.iop.org/article/10.3847/2041-8213/ad2df1)

### Software
7. [eht-imaging (ehtim) GitHub](https://github.com/achael/eht-imaging)
8. [ehtim Documentation](https://achael.github.io/eht-imaging/)
9. [ehtim Imager API](https://achael.github.io/eht-imaging/imager.html)
10. [ehtim Image API](https://achael.github.io/eht-imaging/image.html)
11. [ehtim Obsdata API](https://achael.github.io/eht-imaging/obsdata.html)
12. [SMILI GitHub](https://github.com/astrosmili/smili)
13. [eht-dmc (D-term Modeling Code) GitHub](https://github.com/dpesce/eht-dmc)
14. [eht-dmc PyPI](https://pypi.org/project/eht-dmc/)
15. [Comrade.jl Paper](https://www.theoj.org/joss-papers/joss.04457/10.21105.joss.04457.pdf)
16. [StarWarps example](https://github.com/achael/eht-imaging/blob/main/examples/example_starwarps.py)

### Key Method Papers
17. [Chael et al. 2016 - High Resolution Linear Polarimetric Imaging for the EHT](https://arxiv.org/abs/1605.06156)
18. [Pesce 2021 - DMC for Simultaneous Calibration and Full-Stokes Imaging](https://arxiv.org/abs/2102.03328)
19. [Bouman et al. 2017 - StarWarps: Reconstructing Video from Interferometric Measurements](https://arxiv.org/abs/1711.01357)
20. [Broderick et al. 2020 - THEMIS: Parameter Estimation Framework for EHT](https://iopscience.iop.org/article/10.3847/1538-4357/ab91a4)

### QED/Axion Physics Papers
21. [Banks & Fischler 2026 - Proposal to Search for the CP Violating Electromagnetic Vacuum Angle at the EHT (arXiv:2601.20965)](https://arxiv.org/abs/2601.20965)
22. [Chen et al. 2022 - Probing Axions via Light Circular Polarization and EHT (arXiv:2209.13572)](https://arxiv.org/abs/2209.13572)
23. [Axion constraints from M87 EHT polarimetry (Nature Astronomy 2022)](https://www.nature.com/articles/s41550-022-01620-3)
24. [Probing Black Hole Magnetic Fields with QED](https://www.mdpi.com/2075-4434/6/2/57)
25. [Vacuum birefringence and X-ray polarization from BH accretion disks](https://arxiv.org/abs/1803.03798)
26. [Circular Polarization of Simulated Images of Black Holes (2024)](https://arxiv.org/abs/2406.15653)

### Scattering and Variability
27. [Johnson et al. 2018 - Model for Anisotropic Interstellar Scattering](https://ui.adsabs.harvard.edu/abs/2018arXiv180501242P)
28. [EHT Expanding Sgr A* dynamical imaging with African extension (A&A 2023)](https://www.aanda.org/articles/aa/full_html/2023/04/aa45344-22/aa45344-22.html)

---

## Key Findings

### 1. EHT Data Access

**Sgr A* Polarization Data (2024 Release):**
- **Repository:** `eventhorizontelescope/2024-D02-01` on GitHub
- **Also hosted on CyVerse:** `EHTC_FirstSgrAPol_Mar2024`
- **Observation dates:** April 6-7, 2017
- **Format:** UVFITS files
- **Two pipelines:** EHT-HOPS (`hops_data/`) and rPICARD/CASA (`casa_data/`)
- **Frequency bands:** Two 2-GHz-wide bands centered at **227.1 GHz (low band)** and **229.1 GHz (high band)**
- **Time averaging:** 10 seconds
- **Frequency averaging:** Averaged across 32 IFs (but lo and hi band available separately)
- **Calibration includes:** Absolute EVPA calibration from ALMA, D-term calibration using M87-derived leakage solutions, R/L gain corrections
- **Note:** R/L gains calibrated **assuming zero intrinsic Stokes V** -- this is important for any Stokes V science
- **File naming:** `[Data-Release-Tag]_[Source]_[year]_[day-of-year]_[band]_[pipeline]_[stages].[format]`

**Stokes I only data (2022 Release):**
- `eventhorizontelescope/2022-D02-01` -- ALL polarization information removed (RL, LR set to zero)
- Not suitable for polarimetric work

**Key point:** The 2024-D02-01 data **does** include full polarization (RR, LL, RL, LR), but frequency has been averaged within each band. For single-frequency imaging, use only the lo-band or hi-band UVFITS separately.

### 2. EHT Imaging Software

**eht-imaging (ehtim):**
- Most practical and well-documented tool
- Python library: `pip install ehtim`
- Full Stokes (I, Q, U, V) imaging capability
- Key Imager methods:
  - `make_image_I()` -- Stokes I only
  - `make_image_P()` -- Stokes P (linear polarization) only
  - `make_image_IP()` -- Stokes I and P simultaneously
  - `make_image_IV()` -- Stokes I and V simultaneously
- Polarimetric data terms: `'m'` (polarimetric ratio), `'pvis'` (polarimetric visibilities)
- Polarimetric regularizers: `'hw'` (Hadamard-Walsh), `'ptv'` (polarimetric total variation)
- Transform: `'mcv'` (modified convolutional vector) for polarimetric imaging
- Image class has `evpa()`, `lin_polfrac()`, `circ_polfrac()`, `mvec`, `chivec`, `pvec` properties
- StarWarps dynamic imaging built-in
- Scattering deblurring built-in

**SMILI (Sparse Modeling Imaging Library for Interferometry):**
- Python-interfaced library using sparse sampling
- Handles dual polarization (R/L gains separately)
- GitHub: `astrosmili/smili`
- Less well-documented than ehtim

**DMC (D-term Modeling Code) by Pesce:**
- Bayesian full-Stokes imaging
- Simultaneously fits D-terms + image + gains
- Uses HMC sampling via PyMC3
- `pip install eht-dmc`
- GitHub: `dpesce/eht-dmc`
- Best for uncertainty quantification

**THEMIS:**
- Full Bayesian parameter estimation framework
- Used for EHT analysis of M87 and Sgr A*
- C++ based, computationally intensive
- Handles geometric model fitting + imaging

**Comrade.jl:**
- Julia Bayesian VLBI modeling package
- Hybrid imaging methodology
- Composable modeling framework

### 3. Polarimetric Imaging Method

**Chael et al. 2016 method (implemented in ehtim):**
1. First image Stokes I using bispectrum (immune to atmospheric phases)
2. Then jointly image Stokes Q, U using **polarimetric ratios** (also robust to gain errors)
3. The polarimetric ratio m = (Q + iU) / I is used as the primary data product
4. This ratio cancels most gain errors, making it more robust

**Practical workflow in ehtim:**
1. Load UVFITS data
2. Flag bad data, average
3. First reconstruct Stokes I (using closure quantities: bispectrum, closure amplitudes)
4. Then reconstruct linear polarization using polarimetric ratios
5. Extract EVPA = 0.5 * arctan2(U, Q) and P = sqrt(Q^2 + U^2)

**D-terms (polarization leakage):**
- Already calibrated in the 2024 release using M87-derived solutions
- DMC can further refine D-terms simultaneously with imaging
- Typical D-term accuracy ~2%

### 4. Dynamic Imaging / Time Series

**Sgr A* variability:** Variable on timescales of **minutes** (much shorter than observation duration of hours)

**Approaches:**
- **Snapshot imaging:** Segment data into short time windows, image each independently. Time resolution limited by u,v coverage (~minutes to tens of minutes achievable)
- **StarWarps (Bouman et al. 2017):** Gaussian Markov Model to propagate information across time frames; reconstructs video rather than single image
- **ehtim dynamic imaging:** Uses dynamical regularizers between snapshot frames
- **100-minute snapshot intervals** used in some analyses

**For polarimetric time series:**
- Paper VII combined both days and both bands for imaging (data-limited)
- Time-resolved polarimetry is extremely challenging with 2017 data
- May require StarWarps-like approach or Bayesian methods
- ngEHT will dramatically improve time resolution

### 5. Sgr A* vs M87* Challenges

**Sgr A* specific challenges:**
- **Intraday variability** on minute timescales (M87 is essentially static during observation)
- **Interstellar scattering:** Plasma between Earth and Galactic Center blurs the image
  - Must deblur (divide visibilities by scattering kernel) -- amplifies noise on long baselines
  - ehtim has built-in `deblur` function using Sgr A* scattering model
- **Compact source:** Smaller angular size than M87
- **Faraday rotation:** Strong and time-variable in Sgr A*

**Current Sgr A* polarization results:**
- Linear polarization: highly polarized ring, spiral EVPA pattern, peak ~40% fractional polarization
- Circular polarization: modest ~5-10% dipole structure (negative west, positive east)
- Stokes V resolved structure could NOT be constrained -- only upper limit of 3.7% CP fraction
- Polarization structure is **time variable**

### 6. Physics: QED Birefringence and the Observable C

**Banks & Fischler (2026) -- arXiv:2601.20965:**
This is the key paper motivating this research question.

**The CP-violating electromagnetic vacuum angle:**
- The coupling: (e^2/32pi^2) * theta_EM * integral(E.B d^4x)
- Predicts a "Fischler-Kundu effect": universal Hall current on BH horizon
- The horizon acts with an impedance that differs for left/right circular polarization

**The observable C:**
- **Definition:** C(t) = (1/pi) * [chi(2pi, t) - chi(0, t)]
  - chi(phi, t) = EVPA as function of azimuthal angle around the ring
- **Equivalently:** C(t,omega) = integral_0^{2pi} [Im(P* d_phi P) / (2pi |P|^2)] dphi
  - where P = Q + iU is the complex polarization
- **Properties:** Parity-odd, measures winding number of EVPA around the ring
- **Key result:** Time-averaged <C> != 0 **if and only if** theta_QED != 0
  - The FK mechanism produces a time-independent azimuthal EVPA pattern
  - Plasma effects (Faraday rotation) average to zero if RM is azimuthally uniform
- **Current status:** "Current data does not appear to be sufficient" -- may need 345 GHz or more time averaging

**Chen et al. (2022):**
- Axion-like particles create circular polarization via photon-axion scattering in BH magnetic field
- Future circular polarization measurements can constrain ultralight axion parameter space

**M87 axion constraints (Nature Astronomy 2022):**
- Axion cloud around spinning BH rotates EVPA
- EHT M87 data used to constrain axion-photon coupling
- Mass range 10^-21 to 10^-20 eV constrained

---

## Follow-up Questions Investigated
- How exactly are lo/hi band files named in 2024-D02-01? -> File naming convention documented, band is part of filename
- Can you do single-frequency Stokes V imaging? -> Yes but R/L gains were calibrated assuming V=0, which complicates Stokes V science
- Is the data already frequency-averaged within each band? -> Yes, 32 IFs averaged, but lo and hi band kept separate
- What is the time resolution achievable? -> Minutes possible but very noisy; EHT Paper VII combined all data

---

## Gaps in Available Information
- Exact UVFITS filenames in 2024-D02-01 (would need to clone the repo)
- Whether Stokes V can be reliably extracted given the V=0 calibration assumption
- Detailed DMC usage examples for Sgr A* polarimetric imaging
- Whether anyone has published time-resolved polarimetric images of Sgr A*
- Specific code for computing the C observable from ehtim images
