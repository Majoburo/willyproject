"""
Azimuthal Fourier decomposition of P(φ) on the photon ring.

Decomposes the complex polarization P(φ) = Q + iU into azimuthal modes:
    P(φ) = Σ_m  a_m · exp(i·m·φ)

The winding number C corresponds to the dominant mode m.  If m=1 stands
out above the noise floor, it is a coherent signal (consistent with FK);
if the power spectrum is flat or red, it is more likely turbulent Faraday.

Usage:
  # From the command line — reads saved pol FITS:
  python azimuthal_spectrum.py output/April11_hi_C+0_pol.fits --ring-r 21

  # From another script:
  from azimuthal_spectrum import azimuthal_power_spectrum, plot_spectrum
  modes, power, phase = azimuthal_power_spectrum(phi, P, good)
  plot_spectrum(modes, power, phase, "April11_hi")
"""

import argparse
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path


def azimuthal_power_spectrum(phi, P, good, m_max=10):
    """Compute azimuthal Fourier decomposition of P(φ) and exp(2iχ).

    Parameters
    ----------
    phi : 1-d array, azimuthal angles (rad)
    P   : 1-d complex array, P = Q + iU at each φ bin
    good : bool mask — bins to include
    m_max : highest mode to compute

    Returns
    -------
    modes : array of int, [−m_max … +m_max]
    power_P : |a_m|² for P(φ) (normalised)
    power_chi : |a_m|² for exp(2iχ(φ)) — EVPA phase only (normalised)
    a_m_P : complex Fourier coefficients of P
    a_m_chi : complex Fourier coefficients of exp(2iχ)
    """
    phi_g = phi[good]
    P_g = P[good]

    # EVPA phase factor: exp(2iχ) = P/|P|, removes brightness modulation
    absP = np.abs(P_g)
    phase_factor = np.where(absP > 0, P_g / absP, 0.0)

    modes = np.arange(-m_max, m_max + 1)
    a_m_P = np.zeros(len(modes), dtype=complex)
    a_m_chi = np.zeros(len(modes), dtype=complex)

    dphi = np.median(np.diff(phi_g))
    for i, m in enumerate(modes):
        basis = np.exp(-1j * m * phi_g) * dphi / (2 * np.pi)
        a_m_P[i] = np.sum(P_g * basis)
        a_m_chi[i] = np.sum(phase_factor * basis)

    power_P = np.abs(a_m_P) ** 2
    power_chi = np.abs(a_m_chi) ** 2
    tot_P = power_P.sum()
    tot_chi = power_chi.sum()
    if tot_P > 0:
        power_P /= tot_P
    if tot_chi > 0:
        power_chi /= tot_chi

    return modes, power_P, power_chi, a_m_P, a_m_chi


def plot_spectrum(modes, power_P, power_chi, a_m_P, a_m_chi, label,
                  outdir="output"):
    """Plot azimuthal power spectrum for P(φ) and EVPA phase factor."""
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Top row: P(φ)
    ax = axes[0, 0]
    ax.bar(modes, power_P, color="steelblue", edgecolor="k", linewidth=0.5)
    peak_m = modes[np.argmax(power_P)]
    ax.axvline(peak_m, color="r", ls="--", lw=0.8, label=f"peak m={peak_m}")
    ax.set_ylabel("Fractional power")
    ax.set_title("P(φ) = Q+iU  power spectrum")
    ax.legend(fontsize=9)

    ax = axes[0, 1]
    ax.semilogy(modes, power_P, "o-", ms=4, color="steelblue")
    ax.axvline(peak_m, color="r", ls="--", lw=0.8)
    ax.set_title("P(φ) power (log)")

    # Bottom row: exp(2iχ) — EVPA phase only
    ax = axes[1, 0]
    ax.bar(modes, power_chi, color="darkorange", edgecolor="k", linewidth=0.5)
    peak_m_chi = modes[np.argmax(power_chi)]
    ax.axvline(peak_m_chi, color="r", ls="--", lw=0.8,
               label=f"peak m={peak_m_chi}")
    ax.set_xlabel("Azimuthal mode m")
    ax.set_ylabel("Fractional power")
    ax.set_title("exp(2iχ)  power spectrum  (brightness removed)")
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    ax.semilogy(modes, power_chi, "o-", ms=4, color="darkorange")
    ax.axvline(peak_m_chi, color="r", ls="--", lw=0.8)
    ax.set_xlabel("Azimuthal mode m")
    ax.set_title("exp(2iχ) power (log)")

    plt.suptitle(f"Azimuthal decomposition — {label}\n"
                 f"FK signal: C=1 → m=+2 in exp(2iχ)", fontsize=13)
    plt.tight_layout()
    plt.savefig(str(outdir / f"{label}_azspec.png"), dpi=150)
    plt.close()
    print(f"  Saved {outdir / f'{label}_azspec.png'}")

    # Print summary for both
    for name, power, a_m in [("P(φ)", power_P, a_m_P),
                              ("exp(2iχ)", power_chi, a_m_chi)]:
        print(f"  Top modes [{name}]:")
        order = np.argsort(power)[::-1]
        for idx in order[:5]:
            m = modes[idx]
            print(f"    m={m:+d}:  power={power[idx]:.3f}  "
                  f"phase={np.degrees(np.angle(a_m[idx])):+.1f}°")


def extract_P_from_fits(fits_path, ring_r_uas=21.0, annulus_hw_uas=5.0,
                        nphi=360, p_mask_frac=0.05):
    """Extract P(φ) from a polarimetric FITS image (standalone helper)."""
    import ehtim as eh

    im = eh.image.load_fits(str(fits_path))
    I = im.imarr()
    Q = im.imarr(pol="Q")
    U = im.imarr(pol="U")
    ny, nx = I.shape
    psize = im.psize

    x = (np.arange(nx) - (nx - 1) / 2.0) * psize
    y = ((ny - 1) / 2.0 - np.arange(ny)) * psize
    X, Y = np.meshgrid(x, y)

    r0 = ring_r_uas * eh.RADPERUAS
    hw = annulus_hw_uas * eh.RADPERUAS
    R = np.sqrt(X ** 2 + Y ** 2)
    annulus = (R >= r0 - hw) & (R <= r0 + hw)
    phi_pix = np.arctan2(Y, X) % (2 * np.pi)

    edges = np.linspace(0, 2 * np.pi, nphi + 1)
    phi = 0.5 * (edges[:-1] + edges[1:])
    P_phi = np.zeros(nphi, dtype=complex)
    weights = np.zeros(nphi)

    k = np.clip(np.digitize(phi_pix[annulus], edges) - 1, 0, nphi - 1)
    w = np.clip(I[annulus], 0, None)
    Pvals = Q[annulus] + 1j * U[annulus]

    for j in range(nphi):
        m = k == j
        if m.any():
            weights[j] = w[m].sum()
            P_phi[j] = np.average(Pvals[m], weights=w[m])

    absP = np.abs(P_phi)
    good = (absP > p_mask_frac * np.nanmax(absP)) & (weights > 0)
    return phi, P_phi, good


# ── CLI ───────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Azimuthal Fourier decomposition of polarimetric FITS")
    parser.add_argument("fits", nargs="+", help="Polarimetric FITS file(s)")
    parser.add_argument("--ring-r", type=float, default=21.0,
                        help="Ring radius in µas (default: 21 for M87)")
    parser.add_argument("--hw", type=float, default=5.0,
                        help="Annulus half-width in µas")
    parser.add_argument("--m-max", type=int, default=10,
                        help="Highest azimuthal mode (default: 10)")
    parser.add_argument("--outdir", default="output")
    args = parser.parse_args()

    for fpath in args.fits:
        label = Path(fpath).stem.replace("_pol", "")
        print(f"\n{'='*50}\n  {label}\n{'='*50}")
        phi, P, good = extract_P_from_fits(fpath, ring_r_uas=args.ring_r,
                                           annulus_hw_uas=args.hw)
        modes, power_P, power_chi, a_m_P, a_m_chi = \
            azimuthal_power_spectrum(phi, P, good, m_max=args.m_max)
        plot_spectrum(modes, power_P, power_chi, a_m_P, a_m_chi,
                      label, outdir=args.outdir)
