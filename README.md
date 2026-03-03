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

Main analysis pipeline that processes EHT data through the following steps:

1. **Load & prepare** — Load UVFITS, flag stations, coherently average, rescale noise
2. **Image Stokes I** — Multi-round imaging with recentering (amp+cphase for pre-calibrated M87, logcamp+selfcal for Sgr A*)
3. **Image Q,U** — Polarimetric imaging using gain-robust `m` data term with gentle `pvis` stabilizer, initialized with explicit EVPA winding patterns
4. **Extract P(φ)** — Sample Q+iU in annulus at fixed ring radius (21 µas M87, 25 µas Sgr A*)
5. **Compute C** — Phase-increment winding number: C = (1/2π) Σ arg(P[i+1]·P[i]*)
6. **Bayesian posterior** — χ²-based posterior over integer C values across multiple winding initializations and random seeds

**Supported sources:** M87 (`--source m87`) and Sgr A* (`--source sgra`).

**Usage:**
```bash
# Run all M87 observation days × both frequency bands
python compute_C.py --source m87

# Single observation
python compute_C.py --source m87 --day April11 --band hi

# Sgr A* with 30-min time frames (default)
python compute_C.py --source sgra
```

**Output:**
- `output/{day}_{band}_C{n}_pol.fits` — Polarimetric image per winding basin
- `output/{day}_{band}_C{n}.png` — Diagnostic plots (Stokes I + annulus, |P(φ)|, EVPA unwrapped, EVPA ticks)
- Console summary with per-observation and combined posteriors

### `azimuthal_spectrum.py`

Modular diagnostic that decomposes P(φ) into azimuthal Fourier modes to characterize the nature of the polarization structure:

```
P(φ) = Σ_m  a_m · exp(i·m·φ)
```

Two decompositions are computed:
- **P(φ)** directly — includes brightness modulation from the asymmetric crescent
- **exp(2iχ(φ))** — EVPA phase only, with brightness removed

An FK signal (C = n) produces a peak at mode m = 2n in the exp(2iχ) spectrum. If m = +2 stands out above the noise floor, it is consistent with a coherent C = 1 winding. A flat or red spectrum indicates turbulent plasma or imaging artifacts.

**Usage:**
```bash
# Analyze one or more polarimetric FITS files
python azimuthal_spectrum.py output/April11_hi_C+0_pol.fits --ring-r 21

# Sgr A* (different ring radius)
python azimuthal_spectrum.py output/April06_lo_f0_C+0_pol.fits --ring-r 25

# Multiple files at once
python azimuthal_spectrum.py output/*_pol.fits --ring-r 21
```

**Library API:**
```python
from azimuthal_spectrum import azimuthal_power_spectrum, plot_spectrum
modes, power_P, power_chi, a_m_P, a_m_chi = azimuthal_power_spectrum(phi, P, good)
plot_spectrum(modes, power_P, power_chi, a_m_P, a_m_chi, "label")
```

**Output:**
- `output/{label}_azspec.png` — 2×2 panel: P(φ) and exp(2iχ) power spectra (linear + log)

### `test_synthetic.py`

Synthetic validation tests. Creates ring images with known EVPA winding patterns (C = 0, 1, −1, 2), runs the extraction pipeline, and verifies that the computed C matches the input.

**Usage:**
```bash
python test_synthetic.py
```

## Installation

Requires Python 3.11+ and the EHT imaging library (`ehtim`).

### Python packages

```bash
# Using uv (recommended)
uv sync

# Or with pip
pip install ehtim numpy scipy matplotlib astropy h5py pandas networkx
```

### NFFT (optional but recommended)

The NFFT library enables fast non-uniform FFTs (`ttype="nfft"`) in ehtim,
which is significantly faster than `ttype="direct"` and more accurate than
`ttype="fast"` for polarimetric imaging. Without it, polarimetric imaging
is limited to `ttype="direct"`.

#### macOS (Apple Silicon)

GCC has an internal compiler error with NFFT's OpenMP code on ARM macOS.
Use Clang with Homebrew's `libomp` instead.

**1. Install dependencies:**

```bash
brew install fftw libomp
```

**2. Build and install the NFFT C library:**

```bash
cd /tmp
curl -L https://github.com/NFFT/nfft/releases/download/3.5.3/nfft-3.5.3.tar.gz -o nfft-3.5.3.tar.gz
tar xzf nfft-3.5.3.tar.gz
cd nfft-3.5.3
./configure --enable-all --enable-openmp \
  --with-fftw3=$(brew --prefix fftw) \
  CC=clang \
  CFLAGS="-I$(brew --prefix fftw)/include -I$(brew --prefix libomp)/include -Xpreprocessor -fopenmp" \
  LDFLAGS="-L$(brew --prefix fftw)/lib -L$(brew --prefix libomp)/lib -lomp" \
  LIBS="-lfftw3_threads"
make -j$(sysctl -n hw.ncpu)
sudo make install
```

**3. Build and install pyNFFT (patched for Cython 3 / Python 3.11+):**

The upstream `pynfft` package is incompatible with modern Python/Cython.
Clone it and apply three fixes before building:

```bash
cd /tmp
git clone https://github.com/pyNFFT/pyNFFT.git
cd pyNFFT
```

**Fix 1 — Cython 3 relative imports.** In `pynfft/nfft.pxd`, change:
```
from cnfft3 cimport nfft_plan
```
to:
```
from pynfft.cnfft3 cimport nfft_plan
```

In `pynfft/nfft.pyx`, change:
```
from cnfft3 cimport *
```
to:
```
from pynfft.cnfft3 cimport *
```

**Fix 2 — nogil blocks.** In `pynfft/nfft.pyx`, Cython 3 disallows
`&self._plan` inside `nogil` blocks. For each of the five `cdef` methods
(`_precompute`, `_trafo`, `_trafo_direct`, `_adjoint`, `_adjoint_direct`),
change from:
```python
cdef void _trafo(self):
    with nogil:
        nfft_trafo(&self._plan)
```
to:
```python
cdef void _trafo(self):
    cdef nfft_plan *p = &self._plan
    with nogil:
        nfft_trafo(p)
```

**Fix 3 — Remove `util` extension.** In `setup.py`, delete the
`ext_modules.append(...)` block for `pynfft.util` in both
`get_extensions()` and `get_cython_extensions()` (it references symbols
not present in NFFT 3.5.3).

**Re-cythonize and install:**

```bash
pip install cython numpy  # build dependencies
python -m cython pynfft/nfft.pyx

CFLAGS="-I$(brew --prefix fftw)/include -I/usr/local/include" \
LDFLAGS="-L$(brew --prefix fftw)/lib -L/usr/local/lib -L$(brew --prefix libomp)/lib \
  -lnfft3_threads -lfftw3 -lfftw3_threads -lomp -lm" \
pip install --no-build-isolation .
```

**Verify:**

```bash
python -c "from pynfft.nfft import NFFT; print('pynfft OK')"
python -c "import ehtim; print('ehtim loaded')"  # should NOT print "No NFFT installed!"
```

## Data

Expected data directory structure:
```
data/
├── m87/hops_data/
│   ├── April05/
│   ├── April06/
│   ├── April10/
│   └── April11/
│       ├── *_101_lo_hops_zbl-dtcal+selfcal.uvfits
│       └── *_101_hi_hops_zbl-dtcal+selfcal.uvfits
└── sgra/hops_data/
    ├── April06/
    └── April07/
        ├── *_097_lo_hops_*dtermcal.uvfits
        └── *_097_hi_hops_*dtermcal.uvfits
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
