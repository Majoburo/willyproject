"""
Polarimetric imaging of Sgr A* from EHT 2017 data.

Produces Stokes I and linear polarization (Q, U) images at a single
frequency band, with EVPA ticks overlaid. Designed to later compute
the Banks-Fischler C observable on rings around the photon ring.
"""

import ehtim as eh
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────────
DATA_DIR = Path("data/hops_data")
BAND = "lo"  # "lo" (227.1 GHz) or "hi" (229.1 GHz)
DAY = "April06"  # "April06" (DOY 096) or "April07" (DOY 097)
DOY = "096" if DAY == "April06" else "097"

NPIX = 64
FOV_UAS = 120  # microarcseconds
ZBL = 2.0  # total compact flux (Jy)
AVG_TIME = 120.0  # coherent averaging time (seconds)

OUT_DIR = Path("output")
OUT_DIR.mkdir(exist_ok=True)


# ── Step 1: Load and prepare data ─────────────────────────────────
uvfits = list(DATA_DIR.glob(f"{DAY}/*_{DOY}_{BAND}_hops_*.uvfits"))[0]
print(f"Loading {uvfits}")

obs = eh.obsdata.load_uvfits(str(uvfits), polrep="stokes")
print(f"  {obs.source}, {obs.rf/1e9:.1f} GHz, {len(obs.data)} visibilities")

# Deblur interstellar scattering
obs = obs.deblur()

# Flag JCMT (known polarimetric issues)
obs = obs.flag_sites(["JC"])

# Coherent average to improve SNR
obs.add_scans()
obs = obs.avg_coherent(AVG_TIME, scan_avg=False)
print(f"  After flagging + averaging: {len(obs.data)} visibilities")


# ── Step 2: Image Stokes I ────────────────────────────────────────
fov = FOV_UAS * eh.RADPERUAS
prior = eh.image.make_square(obs, NPIX, fov)
prior = prior.add_gauss(ZBL, (50 * eh.RADPERUAS, 50 * eh.RADPERUAS, 0, 0, 0))

imgr = eh.imager.Imager(
    obs, prior, prior,
    data_term={"bs": 1},
    reg_term={"simple": 1, "flux": 100, "cm": 50},
    maxit=200,
    ttype="direct",  # DFT mode (no NFFT library needed)
)

# Round 1: coarse
imgr.make_image_I(grads=True, show_updates=False)
im = imgr.out_last().copy()

# Round 2: refine
imgr.init_next = im.blur_circ(30 * eh.RADPERUAS)
imgr.dat_term_next = {"bs": 1}
imgr.reg_term_next = {"tv": 1, "flux": 100, "cm": 50}
imgr.make_image_I(grads=True, show_updates=False)
im_I = imgr.out_last().copy()
print(f"  Stokes I imaging done, total flux = {im_I.total_flux():.2f} Jy")


# ── Step 3: Image linear polarization (Q, U) ──────────────────────
imgr.init_next = im_I.blur_circ(15 * eh.RADPERUAS)
imgr.transform_next = "mcv"
imgr.dat_term_next = {"m": 10}
imgr.reg_term_next = {"hw": 1}
imgr.make_image_P(grads=True, show_updates=False)
im_pol = imgr.out_last().copy()

# Refine
imgr.init_next = im_pol
imgr.dat_term_next = {"m": 100}
imgr.reg_term_next = {"hw": 1, "ptv": 100}
imgr.make_image_P(grads=True, show_updates=False)
im_pol = imgr.out_last().copy()
print(f"  Polarization imaging done")
print(f"  Linear pol fraction: {im_pol.lin_polfrac()*100:.1f}%")
print(f"  EVPA: {np.degrees(im_pol.evpa()):.1f} deg")


# ── Step 4: Save and plot ─────────────────────────────────────────
im_pol.save_fits(str(OUT_DIR / f"sgra_{DAY}_{BAND}_pol.fits"))

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Stokes I
imarr = im_pol.imarr()
axes[0].imshow(imarr, origin="lower", cmap="afmhot", extent=[-FOV_UAS/2, FOV_UAS/2]*2)
axes[0].set_title("Stokes I")
axes[0].set_xlabel("ΔRA (μas)")
axes[0].set_ylabel("ΔDEC (μas)")

# Linear polarization fraction
qarr = im_pol.imarr(pol="Q")
uarr = im_pol.imarr(pol="U")
mL = np.sqrt(qarr**2 + uarr**2) / np.where(imarr > 0, imarr, 1)
im1 = axes[1].imshow(mL, origin="lower", cmap="magma", vmin=0, vmax=0.5,
                      extent=[-FOV_UAS/2, FOV_UAS/2]*2)
axes[1].set_title("Linear pol. fraction")
plt.colorbar(im1, ax=axes[1])

# EVPA
evpa = 0.5 * np.arctan2(uarr, qarr)
im2 = axes[2].imshow(np.degrees(evpa), origin="lower", cmap="hsv", vmin=-90, vmax=90,
                      extent=[-FOV_UAS/2, FOV_UAS/2]*2)
axes[2].set_title("EVPA (deg)")
plt.colorbar(im2, ax=axes[2])

plt.tight_layout()
plt.savefig(str(OUT_DIR / f"sgra_{DAY}_{BAND}_polarization.png"), dpi=150)
print(f"  Saved to {OUT_DIR}/")
plt.close()


# ── Step 5: Image with EVPA ticks ─────────────────────────────────
im_pol.display(plotp=True,
               export_pdf=str(OUT_DIR / f"sgra_{DAY}_{BAND}_evpa_ticks.pdf"),
               show=False)
plt.savefig(str(OUT_DIR / f"sgra_{DAY}_{BAND}_evpa_ticks.png"), dpi=150)
plt.close()
print("Done.")
