#!/usr/bin/env python3
"""
Experiment 1: performance vs. R-TWT wake interval (WI), N = 3 slots, R = 1 and 3.

Data:   data/exp1_wake_interval/exp1_N3_R<R>_WI<WI>us.csv
Output: figures/exp1_wake_interval_2col.{pdf,svg,png}   (2x2, full page width)
        figures/exp1_wake_interval_1col.{pdf,svg,png}   (4x1, one column)
        results/exp1_wake_interval_metrics.csv          (measured + model values)

The defaults reproduce the figures in the paper:
    python analysis/plot_exp1_wake_interval.py
"""

import argparse
import os

import numpy as np

import common
import plot_style

METRICS = ("mean_delay_ms", "jitter_ms", "loss_prob", "pctl_delay_ms")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--S", type=float, default=340.0, help="slot duration S [us]")
    ap.add_argument("--N", type=int, default=3, help="SP duration [slots]")
    ap.add_argument("--p-err", type=float, default=0.22, help="first-attempt error probability")
    ap.add_argument("--q", type=float, default=0.69,
                    help="retransmission error probability of the extended model")
    ap.add_argument("--inv-lambda", type=float, default=16000.0,
                    help="mean packet inter-arrival time [us]")
    ap.add_argument("--K", type=int, default=20, help="buffer size [packets]")
    ap.add_argument("--pctl", type=float, default=99.0, help="delay percentile in panel (d)")
    ap.add_argument("--keep-stalls", action="store_true",
                    help="do not remove receiver-stall / warm-up packets from delay statistics")
    ap.add_argument("--data-dir", default=os.path.join(common.DATA_DIR, "exp1_wake_interval"))
    args = ap.parse_args()

    model_kw = dict(S=args.S, N=args.N, p_err=args.p_err, pctl=args.pctl,
                    inv_lambda=args.inv_lambda, K=args.K)
    Rs = (1, 3)

    # measurements
    meas = {R: [] for R in Rs}
    for key, path in common.find_traces(args.data_dir, r"exp1_N\d+_R(?P<R>\d+)_WI(?P<WI>\d+)us\.csv"):
        if key["R"] in meas:
            m = common.measured_metrics(path, args.pctl, drop_stalls=not args.keep_stalls)
            meas[key["R"]].append((key["WI"] / 1000.0, m))
    for R in Rs:
        meas[R].sort(key=lambda t: t[0])
        print(f"R={R}: {len(meas[R])} runs, WI = {[x for x, _ in meas[R]]} ms")
    wis = sorted({x for s in meas.values() for x, _ in s})

    # model: smooth curves over the WI range
    sweep = np.geomspace(min(wis), max(wis), 80)
    ana = {R: [common.model_metrics(t * 1000.0, R=R, **model_kw) for t in sweep] for R in Rs}
    ext = {R: [common.model_metrics(t * 1000.0, R=R, p_retry=args.q, **model_kw) for t in sweep]
           for R in Rs if R > 1}

    plot_kw = dict(
        xlabel="R-TWT WI [ms]",
        pctl=args.pctl,
        model={R: (list(sweep), ana[R]) for R in Rs},
        meas={R: ([x for x, _ in s], [m for _, m in s]) for R, s in meas.items()},
        extra=[(f"Extended ($q={args.q:g}$), $R={R}$", plot_style.R_COLORS[R], list(sweep), e)
               for R, e in ext.items()],
        xticks=wis, xscale="log2", xtick_fmt="{:.3f}")
    os.makedirs(common.FIG_DIR, exist_ok=True)
    base = os.path.join(common.FIG_DIR, "exp1_wake_interval")
    outs = plot_style.save(plot_style.two_column_figure(**plot_kw), base + "_2col")
    outs += plot_style.save(plot_style.one_column_figure(**plot_kw), base + "_1col")
    print("Figure -> " + ", ".join(os.path.relpath(o, common.ROOT) for o in outs))

    # table: measured and model values at the measured WIs
    rows = []
    for R in Rs:
        for x, m in meas[R]:
            rows.append(["measured", R, x] + [m[k] for k in METRICS])
        for x in wis:
            m = common.model_metrics(x * 1000.0, R=R, **model_kw)
            rows.append(["model", R, x] + [m[k] for k in METRICS])
            if R > 1:
                m = common.model_metrics(x * 1000.0, R=R, p_retry=args.q, **model_kw)
                rows.append([f"extended_model_q{args.q:g}", R, x] + [m[k] for k in METRICS])
    os.makedirs(common.RES_DIR, exist_ok=True)
    common.write_metrics_csv(
        os.path.join(common.RES_DIR, "exp1_wake_interval_metrics.csv"),
        ["source", "R", "wi_ms", "mean_delay_ms", "jitter_ms", "loss_prob",
         f"p{args.pctl:g}_delay_ms"], rows)


if __name__ == "__main__":
    main()
