"""
Generate the aggregate-throughput figure S(T_V) for two saturated stations.

Compares:
  - Original Bankov model (EDCA, RTS/CTS, A-MPDU structure, classic Bianchi
    with p = 1 - (1 - tau)^(N-1), ideal channel) against
  - the Adapted DCF model (DCF, no RTS/CTS, fixed frame, error-prone channel),
  - and the hardware measurement (openwifi, two stations).

Output: figures/throughput_vs_period.{pdf,svg,png} (IEEE single-column style).
Run:    python figure_throughput_vs_period.py
"""

import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "models"))
import dcf_model as dcf          # noqa: E402
import bankov_model as bankov    # noqa: E402

# --- IEEE single-column style (Times New Roman, 9 / 8.5 pt) ---
plt.rcParams.update({
    "font.family": "Times New Roman", "mathtext.fontset": "stix",
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})

COLOR_BANKOV = "#d9700a"
COLOR_DCF = "#1f77b4"


def load_measurement(path):
    """Read (period_us, throughput_Mbps) from the measurement CSV."""
    period, thr = [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if not row.get("throughput_Mbps"):
                continue
            period.append(float(row["period_us"]))
            thr.append(float(row["throughput_Mbps"]))
    return np.array(period), np.array(thr)


def main():
    data_path = os.path.join(HERE, "data", "throughput_vs_period.csv")
    W, Y = load_measurement(data_path)
    grid = np.arange(1200.0, 4700.0, 10.0)

    # Original Bankov: EDCA/AIFS, RTS/CTS, ideal channel, classic Bianchi (N-1).
    bankov_params = bankov.Params3(
        N=2, L=1000, W_min=16, W_max=1024, R=7,
        T_e=9.0, T_SIFS=16.0, T_AIFS=43.0, T_lim=680.0,
        T_RTS_meas=52.0, T_CTS_meas=44.0, T_BAck_meas=44.0,
        data_a=0, data_b=484.0, data_table={1: 484.0},
    )
    S_bankov = np.array([bankov.throughput(bankov_params, float(w)) for w in grid])

    # Adapted DCF: no RTS/CTS, fixed frame, error-prone channel.
    dcf_params = dcf.ParamsDCF(N=2, p_err=0.10, T_data=484.0)
    S_dcf = np.array([dcf.throughput(dcf_params, float(w)) for w in grid])

    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    ax.plot(W / 1000, Y, "x", color="black", ms=4, mew=1.0, linestyle="none",
            label="Measurement (2 STAs)")
    ax.plot(grid / 1000, S_bankov, "--", color=COLOR_BANKOV, lw=1.2,
            label="Original Bankov (RTS/CTS, ideal channel)")
    ax.plot(grid / 1000, S_dcf, "-", color=COLOR_DCF, lw=1.2,
            label=r"Adapted DCF model ($p_{\mathrm{err}}=0.10$)")
    ax.set_xlim(W.min() / 1000 - 0.05, W.max() / 1000 + 0.05)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$T_V$ [ms]")
    ax.set_ylabel(r"Aggregate throughput $S$ [Mbit/s]")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()

    out_dir = os.path.join(HERE, "figures")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "throughput_vs_period")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(out + "." + ext, bbox_inches="tight")
    print("Bankov (no R-TWT) = %.2f Mbit/s, DCF (no R-TWT) = %.2f Mbit/s"
          % (bankov.throughput_no_rtwt(bankov_params),
             dcf.throughput_no_rtwt(dcf_params)))
    print("Saved %s.{pdf,svg,png}" % out)


if __name__ == "__main__":
    main()
