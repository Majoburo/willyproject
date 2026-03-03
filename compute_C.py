"""
Compute the Banks-Fischler C observable (arXiv:2601.20965, Eq. 7)
from EHT polarimetric VLBI data.

C measures the net EVPA winding number around the photon ring.
A non-zero time-averaged ⟨C⟩ indicates the Fischler-Kundu horizon
Hall effect (non-zero θ_QED).

Usage:
  python compute_C.py                              # M87, all days × both bands
  python compute_C.py --day April11 --band lo      # single M87 observation
  python compute_C.py --source sgra                # Sgr A* (time-resolved)
  python compute_C.py --source sgra --frame-sec 1800  # 30-min frames
  python compute_C.py --use-pvis                   # add pvis stabilizer
"""

import argparse
import ehtim as eh
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Source configurations ──────────────────────────────────────────
SOURCES = {
    "m87": {
        "data_dir":    Path("data/m87/hops_data"),
        "day_to_doy":  {"April05": "095", "April06": "096",
                        "April10": "100", "April11": "101"},
        "glob_pattern": "{day}/*_{doy}_{band}_hops_zbl-dtcal+selfcal.uvfits",
        "glob_fallback": "{day}/*_{doy}_{band}_hops_ALMArot.uvfits",
        "zbl":          0.6,    # Jy
        "fov_uas":      128.0,
        "ring_r_uas":   21.0,   # photon ring radius (fixed by BH mass)
        "avg_time":     120.0,  # s — source is static
        "frame_sec":    0,      # no framing: one image per night
        "deblur":       False,
        "selfcal":      False,  # already self-calibrated to SMILI images
        "flag_sites":   ["JC"],
    },
    "sgra": {
        "data_dir":    Path("data/sgra/hops_data"),
        "day_to_doy":  {"April06": "096", "April07": "097"},
        "glob_pattern": "{day}/*_{doy}_{band}_hops_*dtermcal.uvfits",
        "glob_fallback": None,
        "zbl":          2.4,    # Jy
        "fov_uas":      200.0,
        "ring_r_uas":   25.0,   # photon ring radius (fixed by BH mass)
        "avg_time":     60.0,   # s — shorter; source varies on ~20s timescale
        "frame_sec":    1800,   # 30-min frames by default
        "deblur":       True,   # interstellar scattering must be removed
        "selfcal":      True,   # not pre-self-calibrated — do it ourselves
        "flag_sites":   ["JC"],
    },
}

# ── Shared configuration ───────────────────────────────────────────
OUT_DIR = Path("output"); OUT_DIR.mkdir(exist_ok=True)
ALL_BANDS = ["lo", "hi"]

NPIX           = 160
NOISE_RESCALE  = 2.0
SYSTEMATIC_NOISE = 0.02
ANNULUS_HW_UAS = 5.0
NPHI           = 360
P_MASK_FRAC    = 0.05
C_WINDINGS     = [-2, -1, 0, 1, 2]
C_RANGE        = range(-3, 4)


# ── Load data ─────────────────────────────────────────────────────
def load_and_prepare(day, band, cfg, frame_sec=None):
    doy = cfg["day_to_doy"][day]
    pattern = cfg["glob_pattern"].format(day=day, doy=doy, band=band)
    candidates = list(cfg["data_dir"].glob(pattern))
    if not candidates and cfg["glob_fallback"]:
        candidates = list(cfg["data_dir"].glob(
            cfg["glob_fallback"].format(day=day, doy=doy, band=band)))
    uvfits = candidates[0]
    print(f"  Using: {uvfits.name}")

    obs = eh.obsdata.load_uvfits(str(uvfits), polrep="stokes")
    obs = obs.flag_sites(cfg["flag_sites"])
    obs.add_scans()
    obs = obs.avg_coherent(cfg["avg_time"], scan_avg=False)

    if cfg["deblur"]:
        obs = obs.deblur()

    fsec = frame_sec if frame_sec is not None else cfg["frame_sec"]
    if fsec > 0:
        frames = [f for f in obs.split_obs(t_gather=fsec) if len(f.data) >= 30]
        print(f"  Split into {len(frames)} frames (~{fsec:.0f}s each)")
        return frames
    return [obs]


# ── Image Stokes I ────────────────────────────────────────────────
def image_stokes_I(obs, cfg):
    fov = cfg["fov_uas"] * eh.RADPERUAS
    prior = eh.image.make_square(obs, NPIX, fov)
    prior = prior.add_gauss(cfg["zbl"], (40*eh.RADPERUAS, 40*eh.RADPERUAS, 0, 0, 0))
    obs_I = obs.rescale_noise(NOISE_RESCALE)

    if cfg["selfcal"]:
        # Not pre-self-calibrated: start with closure amplitudes only
        # (skip cphase — bispectra computation fails on some short frames),
        # then self-cal before using amplitudes
        imgr = eh.imager.Imager(obs_I, prior, prior,
            data_term={"logcamp": 1},
            reg_term={"simple": 1, "flux": 100, "cm": 50},
            maxit=300, ttype="nfft")
        imgr.make_image_I(grads=True, show_updates=False)
        im = imgr.out_last().copy()

        imgr.init_next = im.blur_circ(20*eh.RADPERUAS)
        imgr.reg_term_next = {"tv": 1, "flux": 100, "cm": 50}
        imgr.make_image_I(grads=True, show_updates=False)
        im = imgr.out_last().copy()

        for blur in [10, 5]:
            obs_sc = eh.self_cal.self_cal(obs, im, method="both",
                                          ttype="nfft", processes=1)
            imgr_sc = eh.imager.Imager(obs_sc, im.blur_circ(blur*eh.RADPERUAS), prior,
                data_term={"amp": 1, "logcamp": 1},
                reg_term={"tv": 2, "flux": 100, "cm": 50},
                maxit=300, ttype="nfft")
            imgr_sc.make_image_I(grads=True, show_updates=False)
            im = imgr_sc.out_last().copy()
            obs = obs_sc
        return imgr_sc, im

    else:
        # Pre-self-calibrated: use amplitudes directly
        imgr = eh.imager.Imager(obs_I, prior, prior,
            data_term={"amp": 1, "cphase": 1},
            reg_term={"simple": 1, "flux": 100, "cm": 50},
            maxit=300, ttype="nfft")
        imgr.make_image_I(grads=True, show_updates=False)
        im = imgr.out_last().copy()

        imgr.init_next = im.blur_circ(10*eh.RADPERUAS)
        imgr.reg_term_next = {"tv": 2, "flux": 100, "cm": 50}
        imgr.make_image_I(grads=True, show_updates=False)
        im = imgr.out_last().copy()

        imgr.init_next = im.blur_circ(5*eh.RADPERUAS)
        imgr.make_image_I(grads=True, show_updates=False)
        return imgr, imgr.out_last().copy()


# ── Image polarization ────────────────────────────────────────────
def set_evpa_winding(im, n, pol_frac=0.2):
    """Set Q/U to EVPA = n·φ pattern as optimizer starting point."""
    ny, nx = im.ydim, im.xdim
    y, x = np.mgrid[0:ny, 0:nx]
    x = (x - (nx-1)/2.0) * im.psize
    y = (y - (ny-1)/2.0) * im.psize
    phi = np.arctan2(y, x)

    I_2d = im.imvec.reshape(ny, nx)
    P = pol_frac * np.clip(I_2d, 0, None)
    chi = n * phi

    im_out = im.copy()
    im_out.qvec = (P * np.cos(2*chi)).ravel()
    im_out.uvec = (P * np.sin(2*chi)).ravel()
    return im_out


def image_polarization(imgr, im_I, winding=None, seed=None, use_pvis=False):
    """Image Q,U using gain-robust 'm' data term with explicit winding init."""
    prior = imgr.prior_next
    res = imgr.obs_next.res()

    if winding is not None:
        init_im = set_evpa_winding(im_I.copy(), winding)
    else:
        if seed is not None:
            np.random.seed(seed)
        init_im = im_I.copy()

    # Round 1: m only
    imgr.init_next = init_im
    imgr.prior_next = prior
    imgr.transform_next = "mcv"
    imgr.dat_term_next = {"m": 1}
    imgr.reg_term_next = {"hw": 1}
    imgr.systematic_noise = SYSTEMATIC_NOISE
    imgr.make_image_P(grads=True, show_updates=False)

    # Round 2: higher m weight + polarimetric TV
    imgr.init_next = imgr.out_last().blur_circ(0.5*res, 0.5*res)
    imgr.prior_next = imgr.init_next
    imgr.transform_next = "mcv"
    imgr.dat_term_next = {"m": 5}
    imgr.reg_term_next = {"hw": 1, "ptv": 1}
    imgr.make_image_P(grads=True, show_updates=False)

    # Optional round 3: gentle pvis stabilizer
    if use_pvis:
        imgr.init_next = imgr.out_last().blur_circ(0.5*res, 0.5*res)
        imgr.prior_next = imgr.init_next
        imgr.transform_next = "mcv"
        imgr.dat_term_next = {"m": 5, "pvis": 0.1}
        imgr.reg_term_next = {"hw": 1, "ptv": 1}
        imgr.make_image_P(grads=True, show_updates=False)

    return imgr.out_last().copy()


# ── Find ring center ──────────────────────────────────────────────
def find_ring_center(im, cfg):
    # Ring center fixed at (0,0) by the cm regularizer during imaging.
    # Ring radius fixed by black hole mass — not a free parameter.
    return 0.0, 0.0, cfg["ring_r_uas"] * eh.RADPERUAS


# ── Extract P(φ) on annulus ───────────────────────────────────────
def extract_P_of_phi(im, x0, y0, r0):
    I, Q, U = im.imarr(), im.imarr(pol="Q"), im.imarr(pol="U")
    ny, nx = I.shape
    psize = im.psize

    x = (np.arange(nx) - (nx-1)/2.0) * psize
    y = (np.arange(ny) - (ny-1)/2.0) * psize
    X, Y = np.meshgrid(x, y)

    R = np.sqrt((X-x0)**2 + (Y-y0)**2)
    hw = ANNULUS_HW_UAS * eh.RADPERUAS
    annulus = (R >= r0-hw) & (R <= r0+hw)
    phi_pix = np.arctan2(Y-y0, X-x0) % (2*np.pi)

    edges = np.linspace(0, 2*np.pi, NPHI+1)
    phi = 0.5*(edges[:-1] + edges[1:])
    P_phi = np.zeros(NPHI, dtype=complex)
    weights = np.zeros(NPHI)

    k = np.clip(np.digitize(phi_pix[annulus], edges) - 1, 0, NPHI-1)
    w = np.clip(I[annulus], 0, None)
    Pvals = Q[annulus] + 1j*U[annulus]

    for j in range(NPHI):
        m = k == j
        if m.any():
            weights[j] = w[m].sum()
            P_phi[j] = np.average(Pvals[m], weights=w[m])

    absP = np.abs(P_phi)
    good = (absP > P_MASK_FRAC*np.nanmax(absP)) & (weights > 0)
    return phi, P_phi, good


# ── Compute C (phase-increment formula) ──────────────────────────
def compute_C(phi, P, good):
    """C = (1/2π) Σ arg(P[i+1]·P[i]*) — total EVPA winding."""
    P_g = P[good]
    if len(P_g) < 10:
        return np.nan
    return np.sum(np.angle(P_g[1:] * np.conj(P_g[:-1]))) / (2*np.pi)


# ── Compute χ²_m for a polarimetric image ────────────────────────
def compute_chi2_m(im_pol, obs):
    obs_model = im_pol.observe_same(obs, ttype="nfft", add_th_noise=False)
    m_d = (obs.data["qvis"] + 1j*obs.data["uvis"]) / obs.data["vis"]
    m_m = (obs_model.data["qvis"] + 1j*obs_model.data["uvis"]) / obs_model.data["vis"]
    sigma_m = obs.data["sigma"] / np.abs(obs.data["vis"])
    return float(np.nanmean(np.abs(m_d - m_m)**2 / sigma_m**2))


# ── Bayesian posterior from χ² per winding basin ──────────────────
def posterior_from_chi2(chi2_dict):
    """P(C=n) ∝ exp(-Δχ²/2) where Δχ² = χ²(n) - min(χ²)."""
    if not chi2_dict:
        return {n: 1.0/len(C_RANGE) for n in C_RANGE}
    chi2_min = min(chi2_dict.values())
    raw = {n: (np.exp(-0.5*(chi2_dict[n] - chi2_min)) if n in chi2_dict
               else 1e-10)
           for n in C_RANGE}
    total = sum(raw.values())
    return {n: p/total for n, p in raw.items()}


def combine_posteriors(posteriors):
    """Multiply independent posteriors (flat prior)."""
    K = len(C_RANGE)
    prior = {n: 1.0/K for n in C_RANGE}
    log_post = {n: np.log(prior[n]) for n in C_RANGE}
    for post in posteriors:
        for n in C_RANGE:
            log_post[n] += np.log(post.get(n, 1e-10)) - np.log(prior[n])
    mx = max(log_post.values())
    raw = {n: np.exp(lp - mx) for n, lp in log_post.items()}
    total = sum(raw.values())
    return {n: p/total for n, p in raw.items()}


def bayes_factor_FK(posterior):
    """P(C≠0) / P(C=0)."""
    p0 = posterior.get(0, 1e-10)
    return sum(p for n, p in posterior.items() if n != 0) / p0


def band_consistency(post_lo, post_hi):
    """P(C_lo = C_hi) — should be high if FK signal is real."""
    return sum(post_lo.get(n, 0)*post_hi.get(n, 0) for n in C_RANGE)


def fmt_post(posterior, top_n=3):
    s = sorted(posterior.items(), key=lambda x: -x[1])
    return ", ".join(f"C={n}: {p:.1%}" for n, p in s[:top_n] if p > 0.01)


# ── Plotting ──────────────────────────────────────────────────────
def plot_results(im, x0, y0, r0, phi, P, good, C_val, label):
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    fov = im.fovx() / eh.RADPERUAS
    ext = [-fov/2, fov/2, -fov/2, fov/2]

    axes[0].imshow(im.imarr(), origin="lower", cmap="afmhot", extent=ext)
    theta = np.linspace(0, 2*np.pi, 256)
    hw = ANNULUS_HW_UAS * eh.RADPERUAS
    for r in [r0-hw, r0+hw]:
        axes[0].plot((x0+r*np.cos(theta))/eh.RADPERUAS,
                     (y0+r*np.sin(theta))/eh.RADPERUAS, "c-", lw=0.8)
    axes[0].set_title("Stokes I + annulus")
    axes[1].plot(np.degrees(phi), np.abs(P))
    axes[1].axhline(np.nanmax(np.abs(P))*P_MASK_FRAC, color="r", ls="--", lw=0.8)
    axes[1].set(title="|P(φ)|", xlabel="φ (deg)")
    chi_uw = np.unwrap(np.angle(P[good])) / 2
    axes[2].plot(np.degrees(phi[good]), np.degrees(chi_uw), ".-", ms=2)
    axes[2].set(title=f"EVPA(φ), C = {C_val:.3f}", xlabel="φ (deg)")
    im.display(plotp=True, show=False, axis=axes[3])
    axes[3].set_title("EVPA ticks")

    plt.suptitle(label, fontsize=14)
    plt.tight_layout()
    plt.savefig(str(OUT_DIR / f"{label}.png"), dpi=150)
    plt.close()


# ── Run one observation ──────────────────────────────────────────
def run_one(day, band, cfg, frame_sec=None, windings=C_WINDINGS, use_pvis=False):
    label = f"{day}_{band}"
    print(f"\n{'='*60}\n  {label}\n{'='*60}")

    frames = load_and_prepare(day, band, cfg, frame_sec=frame_sec)

    frame_posteriors = []
    C_timeseries = []

    for fi, obs in enumerate(frames):
        times = obs.data["time"].astype(float)
        t_mid = 0.5*(times.min() + times.max())
        if len(frames) > 1:
            print(f"\n  Frame {fi}/{len(frames)-1}  "
                  f"({len(obs.data)} vis, t={times.min():.2f}-{times.max():.2f} hr)")

        try:
            imgr, im_I = image_stokes_I(obs, cfg)
        except Exception as e:
            import traceback
            print(f"  ⚠ Stokes I failed: {e}")
            traceback.print_exc()
            continue

        basin_chi2 = {}
        for n in windings:
            try:
                im_pol = image_polarization(imgr, im_I, winding=n, use_pvis=use_pvis)
                x0, y0, r0 = find_ring_center(im_pol, cfg)
                phi, P, good = extract_P_of_phi(im_pol, x0, y0, r0)
                C_val = compute_C(phi, P, good)
                chi2_m = compute_chi2_m(im_pol, obs)
            except Exception as e:
                print(f"    w={n:+d}: FAILED ({e})")
                continue

            C_int = int(np.round(C_val))
            print(f"    w={n:+d}: C={C_val:+.3f} → basin {C_int:+d}  "
                  f"χ²_m={chi2_m:.3f}  mL={im_pol.lin_polfrac()*100:.1f}%")

            suffix = f"{label}_f{fi}_w{n:+d}" if len(frames) > 1 else f"{label}_w{n:+d}"
            plot_results(im_pol, x0, y0, r0, phi, P, good, C_val, suffix)
            im_pol.save_fits(str(OUT_DIR / f"{suffix}_pol.fits"))

            if C_int not in basin_chi2 or chi2_m < basin_chi2[C_int]:
                basin_chi2[C_int] = chi2_m

        if basin_chi2:
            fp = posterior_from_chi2(basin_chi2)
            frame_posteriors.append(fp)
            C_map_frame = max(fp, key=fp.get)
            C_timeseries.append((t_mid, C_map_frame, fp))
            if len(frames) > 1:
                print(f"  Frame {fi} (t={t_mid:.2f}hr): {fmt_post(fp)}")

    if not frame_posteriors:
        return {"day": day, "band": band, "C_map": 0, "posterior": {}}

    posterior = combine_posteriors(frame_posteriors) if len(frame_posteriors) > 1 \
        else frame_posteriors[0]
    C_map = max(posterior, key=posterior.get)

    if len(C_timeseries) > 1:
        C_maps = [c for _, c, _ in C_timeseries]
        C_avg = np.mean(C_maps)
        C_err = np.std(C_maps) / np.sqrt(len(C_maps))
        print(f"\n  ⟨C⟩ = {C_avg:.3f} ± {C_err:.3f}  (N={len(C_maps)} frames)")

    print(f"  Posterior: {fmt_post(posterior)}")
    print(f"  MAP: C = {C_map}  (P = {posterior[C_map]:.1%})")

    return {"day": day, "band": band, "C_map": C_map,
            "posterior": posterior, "C_timeseries": C_timeseries}


# ── Main ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute Banks-Fischler C observable from EHT data")
    parser.add_argument("--source", default="m87", choices=SOURCES.keys())
    parser.add_argument("--day", default=None)
    parser.add_argument("--band", default=None)
    parser.add_argument("--frame-sec", type=float, default=None,
                        help="Frame duration in seconds (0=full obs). "
                             "Default: per-source setting.")
    parser.add_argument("--use-pvis", action="store_true")
    args = parser.parse_args()

    cfg = SOURCES[args.source]
    days = [args.day] if args.day else list(cfg["day_to_doy"].keys())
    bands = [args.band] if args.band else ALL_BANDS

    results = []
    for day in days:
        for band in bands:
            try:
                results.append(run_one(day, band, cfg,
                                       frame_sec=args.frame_sec,
                                       use_pvis=args.use_pvis))
            except Exception as e:
                print(f"  FAILED {day} {band}: {e}")

    if not results:
        raise SystemExit("No successful observations")

    print(f"\n{'='*60}\n  RESULTS ({args.source.upper()})\n{'='*60}")
    for r in results:
        print(f"  {r['day']:>8}_{r['band']:<2}  MAP=C={r['C_map']:+d}  "
              f"{fmt_post(r['posterior'])}")

    posteriors = [r["posterior"] for r in results if r.get("posterior")]
    if len(posteriors) > 1:
        combined = combine_posteriors(posteriors)
        C_map = max(combined, key=combined.get)
        bf = bayes_factor_FK(combined)
        print(f"\n  Combined: {fmt_post(combined)}")
        print(f"  MAP: C = {C_map}  (P = {combined[C_map]:.1%})")
        print(f"  Bayes factor (C≠0 / C=0): {bf:.4f}")
        print(f"  → {'FK signal detected!' if bf > 3 else 'No FK signal (C=0 preferred)'}")

        for day in set(r["day"] for r in results):
            lo = next((r for r in results if r["day"]==day and r["band"]=="lo"), None)
            hi = next((r for r in results if r["day"]==day and r["band"]=="hi"), None)
            if lo and hi and lo.get("posterior") and hi.get("posterior"):
                bc = band_consistency(lo["posterior"], hi["posterior"])
                print(f"  Band consistency ({day}): {bc:.3f}")
