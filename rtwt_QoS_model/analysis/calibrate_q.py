#!/usr/bin/env python3
"""
Calibration and validation of the retransmission error model.

Extended model: the first attempt fails with probability p_err, every
retransmission fails with probability q (simplified Gilbert channel), so
    PLR(R) = p_err * q^(R-1).

Calibration (experiment 2, N = 3):
    p_err = PLR(R=1),   q = sqrt(PLR(R=3) / PLR(R=1)).
Validation: predict PLR(R=2) = p_err * q and compare it with the independent
R = 2 run, which is not used for calibration. Also listed: the prediction of the
original model with independent attempts, PLR(R=2) = p_err^2.

Confidence intervals: 95% block bootstrap (blocks of 100 packets), because
losses come in bursts and are not independent.

    python analysis/calibrate_q.py
"""

import argparse
import os

import numpy as np

import common


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=3, help="SP duration used for calibration")
    ap.add_argument("--data-dir", default=os.path.join(common.DATA_DIR, "exp2_sp_duration"))
    args = ap.parse_args()

    traces = common.find_traces(args.data_dir, r"exp2_WI\d+us_R(?P<R>\d+)_N(?P<N>\d+)\.csv")
    plr = {}
    print(f"{'file':<28} {'packets':>8} {'PLR':>7}   95% CI (block bootstrap)")
    for key, path in traces:
        seqs, _, _ = common.read_trace(path)
        p, lo, hi = common.loss_ci_block_bootstrap(seqs)
        n = common.loss_indicator(seqs).size
        plr[(key["R"], key["N"])] = (p, lo, hi)
        print(f"{os.path.basename(path):<28} {n:>8} {p*100:6.2f}%   [{lo*100:5.2f}, {hi*100:5.2f}]%")

    N = args.N
    p1 = plr[(1, N)][0]
    p3 = plr[(3, N)][0]
    q = np.sqrt(p3 / p1)
    print(f"\nCalibration at N={N}:  p_err = PLR(R=1) = {p1:.3f},  "
          f"q = sqrt(PLR(R=3)/PLR(R=1)) = sqrt({p3:.4f}/{p1:.4f}) = {q:.3f}")

    # q from every N with both R=1 and R=3 runs (consistency check)
    Ns = sorted(n for (r, n) in plr if r == 1 and (3, n) in plr)
    qs = [np.sqrt(plr[(3, n)][0] / plr[(1, n)][0]) for n in Ns]
    print("q per N:  " + ", ".join(f"N={n}: {x:.3f}" for n, x in zip(Ns, qs))
          + f"   (mean {np.mean(qs):.3f})")

    if (2, N) in plr:
        p2, lo, hi = plr[(2, N)]
        print(f"\nValidation, R=2, N={N}:")
        print(f"  measured                          {p2*100:6.2f}%  (95% CI {lo*100:.2f}-{hi*100:.2f}%)")
        for label, val in (("extended model  p_err*q", p1 * q),
                           ("  with figure values 0.22*0.69", 0.22 * 0.69),
                           ("original model  p_err^2 (indep.)", p1 ** 2)):
            inside = "inside" if lo <= val <= hi else "OUTSIDE"
            print(f"  {label:<33} {val*100:6.2f}%  -> {inside} the CI")


if __name__ == "__main__":
    main()
