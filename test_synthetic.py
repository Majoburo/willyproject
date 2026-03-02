"""
Synthetic validation for the C observable pipeline.

Creates ring images with known EVPA winding C, generates synthetic
EHT visibilities, images them, and checks C recovery.
"""

import ehtim as eh
import numpy as np
import matplotlib; matplotlib.use("Agg")
from pathlib import Path

from compute_C import (
    find_ring_center, extract_P_of_phi, compute_C,
    posterior_from_chi2, fmt_post,
    image_polarization, NPIX, FOV_UAS,
)

DATA_DIR = Path("data/m87/hops_data")
UVFITS = list(DATA_DIR.glob("April11/*_101_lo_hops_*.uvfits"))[0]


def make_ring_image(obs, C_input, ring_r_uas=21.0, ring_w_uas=5.0,
                    total_flux=0.6, pol_frac=0.1):
    """Create a model ring with EVPA winding number = C_input."""
    fov = FOV_UAS * eh.RADPERUAS
    im = eh.image.make_square(obs, NPIX, fov)
    psize = im.psize
    x = (np.arange(NPIX) - (NPIX-1)/2.0) * psize
    X, Y = np.meshgrid(x, x)
    R = np.sqrt(X**2 + Y**2)
    phi = np.arctan2(Y, X) % (2*np.pi)

    r0 = ring_r_uas * eh.RADPERUAS
    w = ring_w_uas * eh.RADPERUAS
    I = total_flux * np.exp(-0.5*((R-r0)/w)**2)
    I *= total_flux / I.sum()

    chi = C_input * phi / 2.0
    P = pol_frac * I * np.exp(2j*chi)

    im.imvec = I.ravel()
    im.qvec = np.real(P).ravel()
    im.uvec = np.imag(P).ravel()
    im.vvec = np.zeros_like(im.imvec)
    return im


def image_stokes_I_synthetic(obs_syn):
    """Stokes I imaging + self-cal for synthetic data."""
    fov = FOV_UAS * eh.RADPERUAS
    prior = eh.image.make_square(obs_syn, NPIX, fov)
    prior = prior.add_gauss(0.6, (40*eh.RADPERUAS, 40*eh.RADPERUAS, 0, 0, 0))

    imgr = eh.imager.Imager(obs_syn, prior, prior,
        data_term={"cphase": 1, "logcamp": 1},
        reg_term={"simple": 1, "flux": 100, "cm": 50},
        maxit=300, ttype="direct")
    imgr.make_image_I(grads=True, show_updates=False)
    im = imgr.out_last().copy()

    imgr.init_next = im.blur_circ(20*eh.RADPERUAS)
    imgr.reg_term_next = {"tv": 1, "flux": 100, "cm": 50}
    imgr.make_image_I(grads=True, show_updates=False)
    im = imgr.out_last().copy()

    for blur in [10, 5]:
        obs_sc = eh.self_cal.self_cal(obs_syn, im, method="both",
                                      ttype="direct", processes=1)
        imgr_sc = eh.imager.Imager(obs_sc, im.blur_circ(blur*eh.RADPERUAS), prior,
            data_term={"amp": 1, "cphase": 1},
            reg_term={"tv": 2, "flux": 100, "cm": 50},
            maxit=300, ttype="direct")
        imgr_sc.make_image_I(grads=True, show_updates=False)
        im = imgr_sc.out_last().copy()
        obs_syn = obs_sc

    return imgr_sc


def measure_C(im_pol):
    x0, y0, r0 = find_ring_center(im_pol)
    phi, P, good = extract_P_of_phi(im_pol, x0, y0, r0)
    return compute_C(phi, P, good)


def test_recovery(nseeds=10):
    print("\n" + "="*60)
    print("  C recovery from synthetic data")
    print("="*60)

    obs = eh.obsdata.load_uvfits(str(UVFITS), polrep="stokes")
    obs = obs.flag_sites(["JC"])
    obs.add_scans()
    obs = obs.avg_coherent(120.0, scan_avg=False)

    for C_input in [0, 1, -1, 2]:
        print(f"\n  --- C_input = {C_input} ---")
        im_model = make_ring_image(obs, C_input)
        print(f"  Direct: C = {measure_C(im_model):.4f}")

        obs_syn = im_model.observe_same(obs, ttype="direct",
            add_th_noise=True, ampcal=True, phasecal=True,
            dcal=True, frcal=True, seed=42)

        imgr = image_stokes_I_synthetic(obs_syn)
        im_I = imgr.out_last().copy()

        C_list = []
        for s in range(nseeds):
            im_pol = image_polarization(imgr, im_I, seed=42+s)
            C_val = measure_C(im_pol)
            C_list.append(C_val)
            print(f"    seed {s}: C = {C_val:.4f}")

        # Seed-count posterior (chi2 is identical across windings
        # on synthetic data, so seed counting is the right approach)
        C_ints = [round(c) for c in C_list if np.isfinite(c)]
        counts = {}
        for c in C_ints:
            counts[c] = counts.get(c, 0) + 1
        C_mode = max(counts, key=counts.get)
        ok = "✓" if C_mode == C_input else "✗"
        print(f"  >> MAP={C_mode} ({counts[C_mode]}/{len(C_ints)} seeds) {ok}")


if __name__ == "__main__":
    test_recovery()
