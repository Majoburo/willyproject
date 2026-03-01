# EHT Banks-Fischler C Observable

Computation of the Banks-Fischler C observable (arXiv:2601.20965) from Event Horizon Telescope (EHT) M87 polarimetric VLBI data.

## Overview

This repository implements the topological observable **C** proposed by Banks & Fischler (2026) to probe CP-violating electromagnetic vacuum angles using polarimetric images from supermassive black holes. The observable measures the winding of the polarization field around the black hole shadow:

```
C = (1/2π) ∫ Im(P* dP/dφ) / |P|² dφ
```

where P(φ) = Q + iU is the complex linear polarization extracted along an azimuthal path around the ring.

## Script

### `compute_C.py`

Main analysis pipeline that processes EHT M87 data through the following steps:

1. **Load & prepare** - Load UVFITS, flag stations, coherently average (120s), rescale noise
2. **Image Stokes I** - Self-calibrated imaging using closure amplitudes and closure phases
3. **Image Q,U** - Joint I+P imaging followed by polarimetric ratio imaging
4. **Find ring center** - Optimize ring centroid from Stokes I brightness distribution
5. **Extract P(φ)** - Sample Q+iU in annulus around detected ring
6. **Compute C** - Integrate topological charge density

**Usage:**
```bash
# Run all observation days × both frequency bands
python compute_C.py

# Single observation (April 11, low band)
python compute_C.py --day April11 --band lo

# Multiple random seeds for uncertainty estimation
python compute_C.py --nseeds 10
```

**Output:**
- `output/{day}_{band}_pol.fits` - Polarimetric image cubes
- `output/{day}_{band}.png` - Diagnostic plots (Stokes I, |P(φ)|, EVPA, ticks)
- Console summary table with C values and statistics

## Installation

Requires Python 3.11+ and the EHT imaging library (`ehtim`).

```bash
# Using uv (recommended)
uv sync

# Or with pip
pip install ehtim numpy scipy matplotlib astropy h5py pandas networkx
```

## Data

Expected data directory structure:
```
data/m87/hops_data/
  ├── April05/
  ├── April06/
  ├── April10/
  └── April11/
      ├── *_101_lo_hops_*.uvfits
      └── *_101_hi_hops_*.uvfits
```

**Download EHT M87 data:**
- [EHT Data Portal](https://eventhorizontelescope.org/for-astronomers/data)
- [2017 M87 Polarimetric Data Release](https://github.com/eventhorizontelescope/2021-D01-02)

**Note:** The script expects data organized by observation day (April05, April06, April10, April11) with separate lo/hi band UVFITS files. You may need to reorganize the downloaded data to match the expected structure shown above.

## Background

See the `eht-sgra-polarization-imaging/` directory for detailed research notes on:
- EHT data releases and calibration pipelines
- Polarimetric imaging techniques
- QED/axion physics background
- Software tools (`ehtim`, SMILI, Comrade.jl)

## References

- **Banks & Fischler 2026** - [arXiv:2601.20965](https://arxiv.org/abs/2601.20965) - CP-violating vacuum angle observable
- **EHT Collaboration 2021** - [ApJL 910, L12](https://iopscience.iop.org/article/10.3847/2041-8213/abe71d) - M87 polarization
- **Chael et al. 2016** - [arXiv:1605.06156](https://arxiv.org/abs/1605.06156) - Polarimetric imaging methods

## License

MIT License
