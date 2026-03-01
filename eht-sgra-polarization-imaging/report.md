# Extracting Polarization Images of Sgr A* from EHT Data
## Comprehensive Research Report
**Date:** 2026-03-01

---

## Executive Summary

Generating time-series polarimetric images of Sgr A* at a single EHT frequency is feasible in principle using publicly available data and open-source tools, but the current 2017 data impose severe limitations on time resolution. The **eht-imaging (ehtim)** library is the most practical tool for this task: it supports full-Stokes imaging, handles single-frequency band selection, and includes dynamic imaging capabilities. The public polarized data from the 2024-D02-01 release provides calibrated UVFITS files with separate lo-band (227.1 GHz) and hi-band (229.1 GHz) products. However, Sgr A*'s rapid variability and sparse u,v coverage mean that time-resolved polarimetric imaging pushes against the limits of what the 2017 array can deliver. The CP-violating observable C proposed by Banks & Fischler (2026) is mathematically well-defined but, by the authors' own assessment, likely requires next-generation EHT data or 345 GHz observations.

---

## 1. EHT Data Access

### Available Data Releases

| Release | Repository | Content | Polarization? |
|---------|-----------|---------|---------------|
| 2022-D02-01 | [GitHub](https://github.com/eventhorizontelescope/2022-D02-01) | Sgr A* Stokes I calibrated data | **No** (RL, LR set to zero) |
| 2024-D02-01 | [GitHub](https://github.com/eventhorizontelescope/2024-D02-01) | Sgr A* polarized calibrated data | **Yes** (full RR, LL, RL, LR) |
| CyVerse mirror | [CyVerse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstSgrAPol_Mar2024) | Same as 2024-D02-01 | **Yes** |

### Data Specifications (2024-D02-01)

- **Source:** Sagittarius A*
- **Observation dates:** 2017 April 6 and April 7
- **Frequency bands:** Two 2-GHz-wide bands centered at **227.1 GHz (lo)** and **229.1 GHz (hi)**
- **Format:** UVFITS
- **Time averaging:** 10 seconds
- **Frequency averaging:** 32 IFs averaged within each band (lo and hi kept separate)
- **Pipelines:** Two independent calibration pipelines provide separate UVFITS files:
  - `hops_data/` -- EHT-HOPS pipeline
  - `casa_data/` -- rPICARD/CASA pipeline
- **Stations (8):** ALMA, APEX, SMT (Arizona), JCMT (Hawaii), LMT (Mexico), IRAM 30m (Spain), SMA (Hawaii), South Pole Telescope

### File Naming Convention

```
[Data-Release-Tag]_[Source]_[year]_[day-of-year]_[band]_[pipeline]_[stages].[format]
```

For example, a file might be named:
```
2024-D02-01_SgrA_2017_096_lo_hops_netcal.uvfits
2024-D02-01_SgrA_2017_096_hi_casa_netcal.uvfits
```

where `096` and `097` are day-of-year for April 6 and 7, and `lo`/`hi` denote the frequency band.

### Calibration Applied

The released data includes:
- A-priori amplitude calibration
- Fringe fitting
- **Absolute EVPA calibration** from ALMA
- **D-term (polarization leakage) calibration** using M87-derived leakage solutions
- **R/L gain corrections** for polarimetric gain offsets

**Critical caveat for Stokes V science:** The R/L complex gains were calibrated **assuming zero intrinsic Stokes V visibilities**. This means any intrinsic circular polarization signal may be partially absorbed into the gain calibration. For studies of Stokes V (e.g., the C observable), this assumption must be carefully revisited or circumvented.

### How to Download

```bash
# Clone the polarization data repository
git clone https://github.com/eventhorizontelescope/2024-D02-01.git

# Or download from CyVerse (browser-based):
# https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstSgrAPol_Mar2024
```

---

## 2. Imaging Software Landscape

### Comparison of EHT Imaging Tools

| Tool | Language | Method | Polarization | Dynamic Imaging | Best For |
|------|----------|--------|-------------|-----------------|----------|
| **eht-imaging** | Python | RML (regularized maximum likelihood) | Full Stokes | Yes (StarWarps) | General-purpose; most practical |
| **SMILI** | Python/C | Sparse modeling | Dual-pol gains | Limited | Sparse regularization |
| **DMC** | Python | Bayesian (HMC) | Full Stokes + D-terms | No | Uncertainty quantification |
| **THEMIS** | C++ | Bayesian | Yes | Limited | Model fitting, parameter estimation |
| **Comrade.jl** | Julia | Bayesian (hybrid) | Yes | Yes | Flexible Bayesian modeling |
| **CASA/tclean** | Python/C++ | CLEAN | Yes | No | Standard radio imaging |
| **Difmap** | C | CLEAN + self-cal | Limited | No | Traditional VLBI imaging |

### Recommendation for Time-Series Polarization

**eht-imaging is the clear first choice.** It is the most mature, best-documented Python library with native support for:
- Full-Stokes polarimetric imaging (I, Q, U, V)
- UVFITS loading with frequency/IF selection
- Observation splitting by time/scan
- StarWarps dynamic imaging
- Scattering deblurring for Sgr A*
- EVPA and polarization fraction extraction

For uncertainty quantification, **DMC** can complement eht-imaging by providing Bayesian posterior distributions over images and calibration parameters simultaneously.

### Installation

```bash
# eht-imaging
pip install ehtim

# Optional: faster FFTs (Python <= 3.11)
conda install -c conda-forge pynfft

# DMC (Bayesian full-Stokes)
pip install eht-dmc
```

---

## 3. Polarimetric Imaging: Theory and Practice

### The Polarimetric Ratio Method (Chael et al. 2016)

The standard EHT approach to polarimetric imaging, implemented in ehtim, proceeds in two stages:

1. **Stage 1 -- Stokes I reconstruction:** Use closure quantities (bispectrum, closure amplitudes) that are immune to station-based gain errors and atmospheric phase corruption.

2. **Stage 2 -- Linear polarization reconstruction:** Image Stokes Q and U jointly using the **polarimetric ratio**:

   m = (Q + iU) / I

   This ratio cancels most station-based gain errors, analogous to how closure phases cancel atmospheric phases. The key data product is the complex polarimetric ratio measured on each baseline.

### Handling Calibration

- **D-terms (polarization leakage):** Already calibrated in the 2024 release. If further refinement is needed, DMC can fit D-terms simultaneously with the image. Typical residual D-term errors are ~2%.
- **Gain calibration:** The polarimetric ratio method is inherently robust to gain errors.
- **EVPA calibration:** Absolute EVPA is set by ALMA, which has well-characterized feed rotation angles.
- **Faraday rotation:** Strong and time-variable in Sgr A*. At 230 GHz, the rotation measure (RM) of Sgr A* is ~(-5 x 10^5) rad/m^2, but the Faraday rotation at this frequency is modest (~few degrees).

### Extracting Polarization Products

From a reconstructed ehtim Image object with Stokes I, Q, U, V:

- **EVPA** (Electric Vector Position Angle): chi = 0.5 * arctan2(U, Q)
- **Linear polarization fraction:** m_L = sqrt(Q^2 + U^2) / I
- **Circular polarization fraction:** m_C = V / I
- **Complex polarization:** P = Q + iU
- **Polarization magnitude:** |P| = sqrt(Q^2 + U^2)

These are all available as built-in Image methods in ehtim (see code examples below).

---

## 4. Dynamic Imaging and Time Series

### The Variability Challenge

Sgr A* varies on timescales as short as **minutes**, much shorter than the ~8-hour observing window needed for full u,v coverage. This violates the static-source assumption of standard imaging.

### Strategies

1. **Snapshot imaging:** Divide the observation into short time segments (e.g., ~10-100 minutes), image each independently. The shorter the segment, the worse the u,v coverage and the noisier the image. For the 2017 EHT array with 8 stations, segments shorter than ~30 minutes yield very sparse u,v coverage.

2. **StarWarps (Bouman et al. 2017):** Reconstructs a video (sequence of images) by using a Gaussian Markov Model as a temporal prior. Information propagates between adjacent time frames, allowing recovery of dynamics even with incomplete snapshots. Implemented in ehtim.

3. **Dynamical regularizers (ehtim):** Three regularizers enforce temporal smoothness between consecutive snapshot frames, penalizing rapid changes in image structure.

4. **Stochastic optics imaging:** ehtim's `make_image_I_stochastic_optics` handles scattering mitigation for Sgr A* during imaging.

### Achievable Time Resolution

- With 2017 EHT data: **~10-30 minute snapshots** for Stokes I are feasible but noisy
- For polarimetric snapshots: much more challenging due to lower SNR of polarized signal
- The EHT Paper VII combined both days and both bands to produce a single polarimetric image per day
- Time-resolved polarimetric movies will likely require **ngEHT** (more stations, higher sensitivity)

### What the EHT Collaboration Did (Paper VII)

For the published Sgr A* polarization images:
- Data coherently averaged for **120 seconds**
- High and low bands **combined** into a single dataset
- One image reconstructed per observing day (April 6 and April 7)
- Multiple imaging methods: eht-imaging, THEMIS, snapshot m-ring modeling, DoG-HiT

---

## 5. Code Examples

### 5a. Loading EHT UVFITS Data and Selecting a Single Frequency Band

```python
import ehtim as eh
import numpy as np
import matplotlib.pyplot as plt

# Load a polarimetric UVFITS file (lo-band, April 6, HOPS pipeline)
# Adjust path to your actual data location
obs = eh.obsdata.load_uvfits(
    '2024-D02-01/hops_data/2024-D02-01_SgrA_2017_096_lo_hops_netcal.uvfits',
    polrep='stokes'  # Use Stokes representation (I, Q, U, V)
)

# If your UVFITS has multiple IFs and you want specific ones:
# obs = eh.obsdata.load_uvfits('file.uvfits', IF=[0,1,2,3])  # select specific IFs
# IF='all' averages all IFs (default)

# Inspect the data
print("Source:", obs.source)
print("RA:", obs.ra, "DEC:", obs.dec)
print("RF (Hz):", obs.rf)  # Reference frequency
print("BW (Hz):", obs.bw)  # Bandwidth
print("Stations:", obs.tarr['site'])
print("Number of data points:", len(obs.data))
print("Time range (hr):", obs.data['time'].min(), "-", obs.data['time'].max())

# The observation should contain all four correlation products
# In Stokes representation: vis (I), qvis (Q), uvis (U), vvis (V)
print("Data fields:", obs.data.dtype.names)
```

### 5b. Data Preparation

```python
# Add scan information
obs.add_scans()

# Apply scattering deblurring for Sgr A*
# This divides visibilities by the scattering kernel
obs_deblur = obs.deblur()

# Flag problematic stations if needed (JCMT has known polarimetric issues)
obs_flagged = obs_deblur.flag_sites(['JC'])

# Coherently average to improve SNR (e.g., 120 seconds as in Paper VII)
obs_avg = obs_flagged.avg_coherent(120.0, scan_avg=False)

# Optionally average over scans
# obs_scan = obs_flagged.avg_coherent(0, scan_avg=True)
```

### 5c. Stokes I Imaging

```python
# Set up image parameters
npix = 64              # Number of pixels
fov = 120 * eh.RADPERUAS  # Field of view in radians (120 microarcseconds)
zbl = 2.0              # Total compact flux density (Jy) -- approximate for Sgr A*

# Create a Gaussian prior/initial image
prior = eh.image.make_square(obs_avg, npix, fov)
prior = prior.add_gauss(zbl, (50*eh.RADPERUAS, 50*eh.RADPERUAS, 0, 0, 0))

# Set up the Imager for Stokes I
# Use bispectrum (closure phases + amplitudes) for robustness
imgr = eh.imager.Imager(
    obs_avg, prior, prior,
    data_term={'bs': 1},              # Bispectrum data term
    reg_term={'simple': 1, 'flux': 100, 'cm': 50},  # Regularizers
    maxit=200
)

# Round 1: coarse imaging
imgr.make_image_I(grads=True)
im1 = imgr.out_last().copy()

# Round 2: refine with blurred previous image as initialization
imgr.init_next = im1.blur_circ(30*eh.RADPERUAS)
imgr.dat_term_next = {'bs': 1}
imgr.reg_term_next = {'tv': 1, 'flux': 100, 'cm': 50}
imgr.make_image_I(grads=True)
im_I = imgr.out_last().copy()

# Display Stokes I image
im_I.display()
```

### 5d. Polarimetric Imaging (Stokes Q, U)

```python
# Now image linear polarization using the polarimetric ratio method
# Start from the Stokes I image

# Round 3: initial polarimetric imaging
imgr.init_next = im_I.blur_circ(15*eh.RADPERUAS)
imgr.transform_next = 'mcv'       # Modified convolutional vector transform
imgr.dat_term_next = {'m': 10}    # Polarimetric ratio data term
imgr.reg_term_next = {'hw': 1}    # Hadamard-Walsh regularizer for polarization
imgr.make_image_P(grads=True)     # Image linear polarization (P = Q + iU)
im_pol1 = imgr.out_last().copy()

# Round 4: refined polarimetric imaging with stronger constraints
imgr.init_next = im_pol1
imgr.dat_term_next = {'m': 100}
imgr.reg_term_next = {'hw': 1, 'ptv': 100}  # Add polarimetric total variation
imgr.make_image_P(grads=True)
im_pol = imgr.out_last().copy()

# Display with polarization ticks
im_pol.display(plotp=True)  # Shows EVPA ticks overlaid on Stokes I
```

### 5e. Extracting EVPA and Polarization Fraction Maps

```python
# Get 2D image arrays
imarr = im_pol.imarr()           # Stokes I (2D array, Jy/pixel)
qarr = im_pol.imarr(pol='Q')     # Stokes Q
uarr = im_pol.imarr(pol='U')    # Stokes U
varr = im_pol.imarr(pol='V')    # Stokes V (may be zero if not imaged)

# Compute linear polarization
P_arr = np.sqrt(qarr**2 + uarr**2)  # Linear polarization intensity
mL_arr = P_arr / imarr              # Fractional linear polarization

# EVPA map (radians, measured East of North)
evpa_arr = 0.5 * np.arctan2(uarr, qarr)  # Range: [-pi/2, pi/2]

# Or use built-in methods:
total_evpa = im_pol.evpa()           # Total (integrated) EVPA
total_lin_frac = im_pol.lin_polfrac() # Total fractional linear polarization
total_circ_frac = im_pol.circ_polfrac() # Total fractional circular polarization

# Per-pixel fractional polarization vector (complex: |m|*exp(2i*chi))
mvec = im_pol.mvec   # 1D array of fractional polarization per pixel
chivec = im_pol.chivec  # 1D array of EVPA per pixel (radians)

print(f"Total EVPA: {np.degrees(total_evpa):.1f} deg")
print(f"Linear polarization fraction: {total_lin_frac*100:.1f}%")
print(f"Circular polarization fraction: {total_circ_frac*100:.1f}%")

# Plot EVPA map
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Stokes I
axes[0].imshow(imarr, origin='lower', cmap='afmhot')
axes[0].set_title('Stokes I')

# Linear polarization fraction
im_mL = axes[1].imshow(mL_arr, origin='lower', cmap='magma', vmin=0, vmax=0.5)
axes[1].set_title('Linear Pol. Fraction')
plt.colorbar(im_mL, ax=axes[1])

# EVPA
im_evpa = axes[2].imshow(np.degrees(evpa_arr), origin='lower', cmap='hsv',
                          vmin=-90, vmax=90)
axes[2].set_title('EVPA (deg)')
plt.colorbar(im_evpa, ax=axes[2])

plt.tight_layout()
plt.savefig('sgra_polarization_maps.png', dpi=150)
plt.show()
```

### 5f. Time-Series Polarimetric Imaging

```python
# Split observation into time segments for a time series
obs_full = eh.obsdata.load_uvfits(
    '2024-D02-01/hops_data/2024-D02-01_SgrA_2017_096_lo_hops_netcal.uvfits',
    polrep='stokes'
)
obs_full.add_scans()

# Option 1: Split by scans (natural segmentation from telescope scheduling)
obs_list = obs_full.split_obs()  # One Obsdata per scan

# Option 2: Use StarWarps for dynamic imaging
from ehtim.imaging import starwarps as sw

# Split into time segments
obs_segments = sw.splitObs(obs_full)
print(f"Number of time segments: {len(obs_segments)}")

# Set up StarWarps parameters
npix = 30
fov = 100 * eh.RADPERUAS
flux = 2.0
fwhm = 50 * eh.RADPERUAS

# Initialize prior
emptyprior = eh.image.make_square(obs_full, npix, fov)
gaussprior = emptyprior.add_gauss(flux, (fwhm, fwhm, 0, 0, 0))

meanImg = [gaussprior.copy()]
imCov = [sw.gaussImgCovariance_2(meanImg[0], powerDropoff=2.0, frac=0.5)]
noiseCov_img = np.eye(npix**2) * 1e-7

# Initialize flow basis (no motion initially)
init_x, init_y, flowbasis_x, flowbasis_y, initTheta = \
    sw.affineMotionBasis_noTranslation(meanImg[0])

# Run StarWarps to get time-resolved images
expVal_t, expVal_t_t, _, loglikelihood, apxImgs = sw.computeSuffStatistics(
    meanImg, imCov, obs_segments, noiseCov_img,
    initTheta, init_x, init_y, flowbasis_x, flowbasis_y, initTheta,
    method='phase', measurement={'vis': 1},
    interiorPriors=True, numLinIters=5, compute_expVal_tm1_t=False
)

# expVal_t is a list of Image objects, one per time segment
# Save as movie
sw.movie(expVal_t, out='sgra_dynamic.mp4')

# NOTE: StarWarps as shown here works on Stokes I only.
# For polarimetric time series, you would need to:
# 1. Run StarWarps for Stokes I to get the dynamic I model
# 2. For each time segment, use the I model as a prior and image P
# This is an active area of research.
```

### 5g. Computing the C Observable (Banks-Fischler)

```python
def compute_C_observable(image, ring_radius_uas=25.0, ring_width_uas=10.0, nphi=360):
    """
    Compute the CP-violating observable C from a polarimetric image.

    C = (1/pi) * [chi(2pi) - chi(0)]

    where chi(phi) is the EVPA as a function of azimuthal angle around the ring.

    This is equivalent to the winding number of the EVPA around the photon ring.

    Parameters
    ----------
    image : ehtim.image.Image
        Polarimetric image with Stokes I, Q, U
    ring_radius_uas : float
        Radius of the ring in microarcseconds
    ring_width_uas : float
        Width of the annular region to average over
    nphi : int
        Number of azimuthal samples

    Returns
    -------
    C : float
        The C observable (winding number / pi)
    chi_phi : array
        EVPA as function of azimuthal angle
    phi : array
        Azimuthal angles (radians)
    """
    # Image parameters
    fov = image.fovx()
    npix = image.xdim
    psize = image.psize  # pixel size in radians

    # Get Stokes arrays
    imarr = image.imarr()
    qarr = image.imarr(pol='Q')
    uarr = image.imarr(pol='U')

    # Image center
    cx, cy = npix // 2, npix // 2

    # Ring parameters in pixels
    r_ring = ring_radius_uas * eh.RADPERUAS / psize
    dr = ring_width_uas * eh.RADPERUAS / psize

    # Sample azimuthal angles
    phi = np.linspace(0, 2*np.pi, nphi, endpoint=False)

    chi_phi = np.zeros(nphi)

    for i, p in enumerate(phi):
        # Sample points along this azimuthal angle within the ring width
        r_samples = np.linspace(r_ring - dr/2, r_ring + dr/2, 10)
        q_sum, u_sum = 0.0, 0.0
        w_sum = 0.0

        for r in r_samples:
            # Convert to pixel coordinates
            # Note: azimuthal angle measured East of North
            px = cx + r * np.sin(p)
            py = cy + r * np.cos(p)

            # Bilinear interpolation
            px_int, py_int = int(px), int(py)
            if 0 <= px_int < npix-1 and 0 <= py_int < npix-1:
                fx, fy = px - px_int, py - py_int

                q_val = (qarr[py_int, px_int] * (1-fx)*(1-fy) +
                         qarr[py_int, px_int+1] * fx*(1-fy) +
                         qarr[py_int+1, px_int] * (1-fx)*fy +
                         qarr[py_int+1, px_int+1] * fx*fy)
                u_val = (uarr[py_int, px_int] * (1-fx)*(1-fy) +
                         uarr[py_int, px_int+1] * fx*(1-fy) +
                         uarr[py_int+1, px_int] * (1-fx)*fy +
                         uarr[py_int+1, px_int+1] * fx*fy)
                i_val = (imarr[py_int, px_int] * (1-fx)*(1-fy) +
                         imarr[py_int, px_int+1] * fx*(1-fy) +
                         imarr[py_int+1, px_int] * (1-fx)*fy +
                         imarr[py_int+1, px_int+1] * fx*fy)

                weight = i_val  # Weight by intensity
                q_sum += q_val * weight
                u_sum += u_val * weight
                w_sum += weight

        if w_sum > 0:
            chi_phi[i] = 0.5 * np.arctan2(u_sum / w_sum, q_sum / w_sum)

    # Unwrap the EVPA to track continuous winding
    chi_unwrapped = np.unwrap(2 * chi_phi) / 2  # Unwrap the 2*chi to handle pi ambiguity

    # C = (1/pi) * [chi(2pi) - chi(0)]
    # Use the total accumulated phase change
    C = (chi_unwrapped[-1] - chi_unwrapped[0]) / np.pi

    # Alternative: integral form
    # C_alt = np.trapz(np.gradient(chi_unwrapped, phi), phi) / np.pi

    return C, chi_phi, phi


# Usage:
# C, chi_phi, phi = compute_C_observable(im_pol, ring_radius_uas=25.0)
# print(f"C observable = {C:.4f}")
```

---

## 6. Physics Context

### The Banks-Fischler Proposal (arXiv:2601.20965)

Banks and Fischler (2026) propose using EHT polarization data to search for a non-zero CP-violating electromagnetic vacuum angle theta_EM. The key physics is:

**The Fischler-Kundu (FK) Effect:** A black hole horizon acts as a membrane with impedance that depends on the topological theta-term in electromagnetism. The coupling:

```
L_theta = (e^2 / 32 pi^2) * theta_EM * integral(E . B d^4x)
```

induces a Hall conductivity on the horizon. Left and right circular polarizations experience different impedances:

```
J_{\pm} = (sigma \mp i * theta_QED) * E_{\pm}
```

This helicity asymmetry imprints a persistent, time-independent, parity-odd polarization structure in the emitted radiation.

**The Observable C:** The winding number of the EVPA around the photon ring:

```
C(t) = (1/pi) * [chi(2*pi, t) - chi(0, t)]
```

where chi(phi, t) is the EVPA at azimuthal position phi around the ring. Crucially:
- The FK contribution to C is **time-independent** (set by horizon physics)
- Plasma Faraday rotation contributions **average to zero** over time (unless there is a persistent parity-odd azimuthal RM gradient)
- Therefore: **time-averaged C is nonzero if and only if theta_QED is nonzero**

**Current feasibility:** The authors state "current data does not appear to be sufficient" for this measurement. They suggest 345 GHz observations or extended time averaging with future EHT/ngEHT data may be required.

### Axion Constraints from EHT Polarimetry

Several complementary approaches use EHT polarization to probe fundamental physics:

1. **Axion superradiance clouds (Nature Astronomy 2022):** An axion cloud around a spinning black hole rotates the EVPA. EHT M87* data constrains the axion-photon coupling g_{a gamma gamma} for masses ~10^{-21}--10^{-20} eV.

2. **Axion-induced circular polarization (Chen et al. 2022):** Axion-photon scattering in BH magnetic fields generates circular polarization. Future Stokes V measurements can probe ultralight axion parameter space.

3. **QED vacuum birefringence (Caiazzo & Heyl 2018):** The QED Euler-Heisenberg effective Lagrangian makes the vacuum birefringent in strong magnetic fields. This changes X-ray polarization from BH accretion disks, though the effect at mm wavelengths is very small.

### Circular Polarization of Sgr A* (EHT Paper VII)

The EHT detected a modest (~5-10%) circular polarization dipole in Sgr A*:
- Negative Stokes V in the western ring
- Positive Stokes V in the eastern ring
- The resolved Stokes V structure could **not** be confidently constrained
- Upper limit on resolved CP fraction: 3.7%
- CP in GRMHD simulations encodes the **handedness (vertical twist)** of the magnetic field

---

## 7. Practical Considerations

### Sgr A* vs. M87*: Key Differences

| Challenge | Sgr A* | M87* |
|-----------|--------|------|
| Variability | Minutes | Weeks to months |
| Interstellar scattering | Severe (Galactic Center) | Negligible |
| Faraday rotation | Strong, time-variable | Moderate |
| Angular size | ~52 uas | ~42 uas |
| Flux density (230 GHz) | ~2.4 Jy | ~0.5 Jy |
| Imaging difficulty | Very hard (variability) | Hard (low flux) |

### Interstellar Scattering

The interstellar medium toward the Galactic Center scatters Sgr A*'s emission:
- The scattering kernel is anisotropic and well-characterized
- ehtim's `obs.deblur()` divides visibilities by the scattering kernel
- Deblurring **amplifies noise** on long baselines
- The scattering model from Johnson et al. (2018) is the standard
- ehtim includes `make_image_I_stochastic_optics` for joint scattering + imaging

### Single-Frequency Considerations

Using only one band (lo or hi) means:
- **Half the sensitivity** compared to combining both bands
- **No rotation measure** estimation (requires two frequencies)
- **Cleaner interpretation** for frequency-dependent effects
- The FK effect is frequency-independent, so single-frequency is actually preferred for isolating it from frequency-dependent Faraday rotation
- The two EHT bands (227.1 and 229.1 GHz) are close enough that differences are small for most purposes

### Workflow Summary for the Complete Pipeline

1. Download `2024-D02-01` data
2. Load single-band (lo or hi) UVFITS with ehtim
3. Flag, deblur, and average
4. Reconstruct Stokes I using closure quantities
5. Reconstruct Stokes Q, U using polarimetric ratios
6. Extract EVPA and polarization fraction maps
7. For time series: split into scans, image each (or use StarWarps)
8. For the C observable: extract EVPA along the ring azimuth, compute winding number
9. Time-average C across all available snapshots

---

## Key References

### EHT Sgr A* Polarization
- [EHT Collaboration et al. 2024, Paper VII, ApJL 964, L25](https://iopscience.iop.org/article/10.3847/2041-8213/ad2df0)
- [EHT Collaboration et al. 2024, Paper VIII, ApJL 964, L26](https://iopscience.iop.org/article/10.3847/2041-8213/ad2df1)

### Data
- [2024-D02-01 GitHub Repository](https://github.com/eventhorizontelescope/2024-D02-01)
- [CyVerse Data Mirror](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstSgrAPol_Mar2024)
- [EHT Data Products Page](https://eventhorizontelescope.org/for-astronomers/data)

### Software
- [eht-imaging GitHub](https://github.com/achael/eht-imaging) | [Docs](https://achael.github.io/eht-imaging/) | [PyPI](https://pypi.org/project/ehtim/)
- [eht-dmc GitHub](https://github.com/dpesce/eht-dmc) | [PyPI](https://pypi.org/project/eht-dmc/)
- [SMILI GitHub](https://github.com/astrosmili/smili)

### Method Papers
- [Chael et al. 2016 -- Polarimetric Imaging for EHT](https://arxiv.org/abs/1605.06156)
- [Pesce 2021 -- DMC: Bayesian Full-Stokes Imaging](https://arxiv.org/abs/2102.03328)
- [Bouman et al. 2017 -- StarWarps Dynamic Imaging](https://arxiv.org/abs/1711.01357)
- [Broderick et al. 2020 -- THEMIS](https://iopscience.iop.org/article/10.3847/1538-4357/ab91a4)

### Fundamental Physics
- [Banks & Fischler 2026 -- CP Violating Vacuum Angle at EHT](https://arxiv.org/abs/2601.20965)
- [Chen et al. 2022 -- Probing Axions via Circular Polarization](https://arxiv.org/abs/2209.13572)
- [Axion constraints from M87 EHT polarimetry (Nature Astronomy)](https://www.nature.com/articles/s41550-022-01620-3)
- [Circular Polarization of Simulated BH Images (2024)](https://arxiv.org/abs/2406.15653)
- [Caiazzo & Heyl 2018 -- Vacuum Birefringence and BH Accretion Disks](https://arxiv.org/abs/1803.03798)
