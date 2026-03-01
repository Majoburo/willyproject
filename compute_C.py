"""
Compute the Banks-Fischler C observable (arXiv:2601.20965, eq. 7)
from EHT polarimetric VLBI data.

Pipeline per (day, band):
  1. Load single-band UVFITS, flag, average, rescale noise
  2. Image Stokes I (closure quantities + self-cal)
  3. Image Q,U (joint IP then polarimetric ratios)
  4. Find ring center from Stokes I
  5. Extract P(phi) = Q + iU in annulus
  6. Compute C = (1/2pi) * integral[ Im(P* dP/dphi) / |P|^2 ] dphi

Usage:
  python compute_C.py              # run all days x both bands
  python compute_C.py --day April11 --band lo   # single run
  python compute_C.py --nseeds 10  # 10 random seeds per (day, band)
"""

import argparse
import ehtim as eh
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize

# ── Configuration ──────────────────────────────────────────────────
DATA_DIR = Path("data/m87/hops_data")
DAY_TO_DOY = {"April05": "095", "April06": "096",
              "April10": "100", "April11": "101"}
ALL_DAYS = list(DAY_TO_DOY.keys())
ALL_BANDS = ["lo", "hi"]

NPIX = 64
FOV_UAS = 128.0
ZBL = 0.6
AVG_TIME = 120.0
NOISE_RESCALE = 6.0

RING_RMIN_UAS = 15.0
RING_RMAX_UAS = 35.0
ANNULUS_HW_UAS = 5.0
NPHI = 360
P_MASK_FRAC = 0.05
NSEEDS = 5

OUT_DIR = Path("output")
OUT_DIR.mkdir(exist_ok=True)


# ── 1. Load and prepare ───────────────────────────────────────────
def load_and_prepare(day, band):
    doy = DAY_TO_DOY[day]
    uvfits = list(DATA_DIR.glob(f"{day}/*_{doy}_{band}_hops_*.uvfits"))[0]
    obs = eh.obsdata.load_uvfits(str(uvfits), polrep="stokes")
    obs = obs.flag_sites(["JC"])
    obs.add_scans()
    obs = obs.avg_coherent(AVG_TIME, scan_avg=False)
    obs = obs.rescale_noise(NOISE_RESCALE)
    return obs


# ── 2. Image Stokes I with self-cal ───────────────────────────────
def image_stokes_I(obs):
    fov = FOV_UAS * eh.RADPERUAS
    prior = eh.image.make_square(obs, NPIX, fov)
    prior = prior.add_gauss(ZBL, (40 * eh.RADPERUAS, 40 * eh.RADPERUAS, 0, 0, 0))

    imgr = eh.imager.Imager(
        obs, prior, prior,
        data_term={"cphase": 1, "logcamp": 1},
        reg_term={"simple": 1, "flux": 100, "cm": 50},
        maxit=300, ttype="direct",
    )
    imgr.make_image_I(grads=True, show_updates=False)
    im = imgr.out_last().copy()

    imgr.init_next = im.blur_circ(20 * eh.RADPERUAS)
    imgr.reg_term_next = {"tv": 1, "flux": 100, "cm": 50}
    imgr.make_image_I(grads=True, show_updates=False)
    im = imgr.out_last().copy()

    # Self-cal + re-image with amplitudes
    obs_sc = eh.self_cal.self_cal(obs, im, method="both", ttype="direct", processes=1)
    imgr_sc = eh.imager.Imager(
        obs_sc, im.blur_circ(10 * eh.RADPERUAS), prior,
        data_term={"amp": 1, "cphase": 1},
        reg_term={"tv": 2, "flux": 100, "cm": 50},
        maxit=300, ttype="direct",
    )
    imgr_sc.make_image_I(grads=True, show_updates=False)
    im = imgr_sc.out_last().copy()

    # Second self-cal round
    obs_sc2 = eh.self_cal.self_cal(obs_sc, im, method="both", ttype="direct", processes=1)
    imgr_f = eh.imager.Imager(
        obs_sc2, im.blur_circ(5 * eh.RADPERUAS), prior,
        data_term={"amp": 1, "cphase": 1},
        reg_term={"tv": 2, "flux": 100, "cm": 50},
        maxit=300, ttype="direct",
    )
    imgr_f.make_image_I(grads=True, show_updates=False)
    return imgr_f, imgr_f.out_last().copy()


# ── 3. Image polarization ─────────────────────────────────────────
def image_polarization(imgr, seed=None):
    obs = imgr.obs_next
    prior = imgr.prior_next
    im_I = imgr.out_last().copy()

    imgr_ip = eh.imager.Imager(
        obs, im_I.blur_circ(15 * eh.RADPERUAS), prior,
        data_term={"amp": 1, "cphase": 1},
        reg_term={"tv": 1, "flux": 100, "cm": 50},
        maxit=500, ttype="direct",
    )
    if seed is not None:
        np.random.seed(seed)
    imgr_ip.make_image_IP(grads=True, show_updates=False)

    imgr_ip.init_next = imgr_ip.out_last().copy()
    imgr_ip.transform_next = "mcv"
    imgr_ip.dat_term_next = {"pvis": 5}
    imgr_ip.reg_term_next = {"hw": 1, "ptv": 10}
    imgr_ip.make_image_P(grads=True, show_updates=False)

    imgr_ip.init_next = imgr_ip.out_last().copy()
    imgr_ip.dat_term_next = {"pvis": 10, "m": 5}
    imgr_ip.reg_term_next = {"hw": 1, "ptv": 30}
    imgr_ip.make_image_P(grads=True, show_updates=False)
    return imgr_ip.out_last().copy()


# ── 4. Find ring center ───────────────────────────────────────────
def find_ring_center(im):
    I = im.imarr()
    ny, nx = I.shape
    psize = im.psize
    sig_pix = 3 * eh.RADPERUAS / psize
    I_s = gaussian_filter(I, sig_pix) if sig_pix > 0.5 else I

    x = (np.arange(nx) - (nx - 1) / 2.0) * psize
    y = (np.arange(ny) - (ny - 1) / 2.0) * psize
    X, Y = np.meshgrid(x, y)

    thresh = 0.3 * I_s.max()
    m = I_s > thresh
    x0 = np.average(X[m], weights=I_s[m]) if m.any() else 0.0
    y0 = np.average(Y[m], weights=I_s[m]) if m.any() else 0.0

    rmin = RING_RMIN_UAS * eh.RADPERUAS
    rmax = RING_RMAX_UAS * eh.RADPERUAS

    def neg_ring(xy):
        R = np.sqrt((X - xy[0])**2 + (Y - xy[1])**2)
        mask = (R >= rmin) & (R <= rmax)
        return -np.mean(I_s[mask]) if mask.any() else 1e30

    res = minimize(neg_ring, [x0, y0], method="Nelder-Mead")
    x0, y0 = res.x

    R = np.sqrt((X - x0)**2 + (Y - y0)**2)
    edges = np.linspace(rmin, rmax, 101)
    centers = 0.5 * (edges[:-1] + edges[1:])
    profile = [np.mean(I_s[(R >= edges[j]) & (R < edges[j+1])])
               if ((R >= edges[j]) & (R < edges[j+1])).any() else 0
               for j in range(100)]
    r0 = centers[np.argmax(profile)]
    return x0, y0, r0


# ── 5. Extract P(phi) ─────────────────────────────────────────────
def extract_P_of_phi(im, x0, y0, r0):
    I = im.imarr()
    Q = im.imarr(pol="Q")
    U = im.imarr(pol="U")
    psize = im.psize
    ny, nx = I.shape

    x = (np.arange(nx) - (nx - 1) / 2.0) * psize
    y = (np.arange(ny) - (ny - 1) / 2.0) * psize
    X, Y = np.meshgrid(x, y)

    R = np.sqrt((X - x0)**2 + (Y - y0)**2)
    hw = ANNULUS_HW_UAS * eh.RADPERUAS
    annulus = (R >= r0 - hw) & (R <= r0 + hw)
    phi_pix = np.arctan2(Y - y0, X - x0) % (2 * np.pi)

    edges = np.linspace(0, 2 * np.pi, NPHI + 1)
    phi = 0.5 * (edges[:-1] + edges[1:])
    P_phi = np.zeros(NPHI, dtype=complex)
    weights = np.zeros(NPHI)

    k = np.clip(np.digitize(phi_pix[annulus], edges) - 1, 0, NPHI - 1)
    w = np.clip(I[annulus], 0, None)
    Pvals = Q[annulus] + 1j * U[annulus]

    for j in range(NPHI):
        m = k == j
        if m.any():
            weights[j] = w[m].sum()
            P_phi[j] = np.average(Pvals[m], weights=w[m])

    absP = np.abs(P_phi)
    good = (absP > P_MASK_FRAC * np.nanmax(absP)) & (weights > 0)
    return phi, P_phi, good


# ── 6. Compute C ──────────────────────────────────────────────────
def compute_C(phi, P, good):
    phi_g, P_g = phi[good], P[good]
    if len(phi_g) < 10:
        return np.nan
    dphi = np.diff(phi_g)
    dP = np.diff(P_g)
    P_mid = 0.5 * (P_g[:-1] + P_g[1:])
    absP2 = np.abs(P_mid)**2
    integrand = np.imag(np.conj(P_mid) * dP / dphi) / np.where(absP2 > 0, absP2, 1)
    return np.sum(integrand * dphi) / (2 * np.pi)


# ── 7. Plot ───────────────────────────────────────────────────────
def plot_results(im, x0, y0, r0, phi, P, good, C_val, label):
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    fov = im.fovx() / eh.RADPERUAS
    ext = [-fov/2, fov/2, -fov/2, fov/2]

    axes[0].imshow(im.imarr(), origin="lower", cmap="afmhot", extent=ext)
    theta = np.linspace(0, 2*np.pi, 256)
    hw = ANNULUS_HW_UAS * eh.RADPERUAS
    for r in [r0 - hw, r0 + hw]:
        axes[0].plot((x0 + r*np.cos(theta))/eh.RADPERUAS,
                     (y0 + r*np.sin(theta))/eh.RADPERUAS, "c-", lw=0.8)
    axes[0].set_title("Stokes I + annulus")

    axes[1].plot(np.degrees(phi), np.abs(P))
    axes[1].axhline(np.nanmax(np.abs(P)) * P_MASK_FRAC, color="r", ls="--", lw=0.8)
    axes[1].set_title("|P(φ)|")
    axes[1].set_xlabel("φ (deg)")

    chi = 0.5 * np.angle(P)
    chi_uw = np.unwrap(2 * chi[good]) / 2
    axes[2].plot(np.degrees(phi[good]), np.degrees(chi_uw), ".-", ms=2)
    axes[2].set_title(f"EVPA(φ), C = {C_val:.3f}")
    axes[2].set_xlabel("φ (deg)")

    im.display(plotp=True, show=False, axis=axes[3])
    axes[3].set_title("EVPA ticks")

    plt.suptitle(label, fontsize=14)
    plt.tight_layout()
    plt.savefig(str(OUT_DIR / f"{label}.png"), dpi=150)
    plt.close()


# ── Run one (day, band) ──────────────────────────────────────────
def run_one(day, band, nseeds=1):
    label = f"{day}_{band}"
    print(f"\n{'='*60}")
    print(f"  {label}  ({nseeds} seed{'s' if nseeds > 1 else ''})")
    print(f"{'='*60}")

    obs = load_and_prepare(day, band)
    print(f"  {obs.source}, {obs.rf/1e9:.1f} GHz, {len(obs.data)} vis")

    # Stokes I imaging is deterministic — do it once
    imgr, im_I = image_stokes_I(obs)

    C_list = []
    for s in range(nseeds):
        seed = 42 + s
        im_pol = image_polarization(imgr, seed=seed)
        x0, y0, r0 = find_ring_center(im_pol)
        phi, P, good = extract_P_of_phi(im_pol, x0, y0, r0)
        C_val = compute_C(phi, P, good)
        C_list.append(C_val)

        suffix = f"{label}_s{s}" if nseeds > 1 else label
        print(f"  seed {s}: r={r0/eh.RADPERUAS:.1f} uas, "
              f"good={good.sum()}/{len(good)}, "
              f"mL={im_pol.lin_polfrac()*100:.1f}%, C={C_val:.4f}")

        im_pol.save_fits(str(OUT_DIR / f"{suffix}_pol.fits"))
        plot_results(im_pol, x0, y0, r0, phi, P, good, C_val, suffix)

    C_arr = np.array(C_list)
    C_mean = np.nanmean(C_arr)
    C_std = np.nanstd(C_arr)
    if nseeds > 1:
        print(f"  >> C = {C_mean:.4f} +/- {C_std:.4f} (std over {nseeds} seeds)")

    return {"day": day, "band": band, "C_mean": C_mean, "C_std": C_std,
            "C_all": C_list, "nseeds": nseeds}


# ── Main ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--day", default=None, help="Single day (e.g. April11)")
    parser.add_argument("--band", default=None, help="Single band (lo or hi)")
    parser.add_argument("--nseeds", type=int, default=NSEEDS,
                        help=f"Random seeds per (day, band) (default {NSEEDS})")
    args = parser.parse_args()

    days = [args.day] if args.day else ALL_DAYS
    bands = [args.band] if args.band else ALL_BANDS

    results = []
    for day in days:
        for band in bands:
            try:
                r = run_one(day, band, nseeds=args.nseeds)
                results.append(r)
            except Exception as e:
                print(f"  FAILED {day} {band}: {e}")
                results.append({"day": day, "band": band,
                                "C_mean": np.nan, "C_std": np.nan,
                                "C_all": [], "nseeds": args.nseeds})

    # Summary table
    print(f"\n{'='*60}")
    print(f"  SUMMARY  ({args.nseeds} seed{'s' if args.nseeds > 1 else ''} each)")
    print(f"{'='*60}")
    print(f"{'Day':<12} {'Band':<6} {'C_mean':>8} {'C_std':>8} {'N':>4}")
    print("-" * 42)
    for r in results:
        n = len([c for c in r.get("C_all", []) if np.isfinite(c)])
        print(f"{r['day']:<12} {r['band']:<6} {r['C_mean']:>8.3f} "
              f"{r['C_std']:>8.3f} {n:>4}")

    # Pool all C values across days, bands, and seeds
    all_C = [c for r in results for c in r.get("C_all", []) if np.isfinite(c)]
    if all_C:
        all_C = np.array(all_C)
        print(f"\n  Grand mean C = {np.mean(all_C):.4f} "
              f"+/- {np.std(all_C)/np.sqrt(len(all_C)):.4f} "
              f"(N={len(all_C)})")
        print(f"  Grand std C  = {np.std(all_C):.4f}")
        for band in ALL_BANDS:
            Cb = [c for r in results if r["band"] == band
                  for c in r.get("C_all", []) if np.isfinite(c)]
            if Cb:
                Cb = np.array(Cb)
                print(f"  Mean C ({band}) = {np.mean(Cb):.4f} "
                      f"+/- {np.std(Cb)/np.sqrt(len(Cb)):.4f} (N={len(Cb)})")
