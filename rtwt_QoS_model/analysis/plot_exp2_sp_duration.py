#!/usr/bin/env python3
"""
Experiment 2: performance vs. R-TWT SP duration T_SP = N slots, WI = 8.192 ms,
R = 1 and 3.

Data:   data/exp2_sp_duration/exp2_WI8192us_R<R>_N<N>.csv
Output: figures/exp2_sp_duration_2col.{pdf,svg,png}   (2x2, full page width)
        figures/exp2_sp_duration_1col.{pdf,svg,png}   (4x1, one column)
        results/exp2_sp_duration_metrics.csv          (measured + model values)

The R = 2 run (exp2_WI8192us_R2_N3.csv) is not plotted; it is used by
calibrate_q.py to validate the extended error model.

The defaults reproduce the figures in the paper:
    python analysis/plot_exp2_sp_duration.py
"""

import argparse
import os

import common
import plot_style

METRICS = ("mean_delay_ms", "jitter_ms", "loss_prob", "pctl_delay_ms")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--S", type=float, default=340.0, help="slot duration S [us]")
    ap.add_argument("--wi-us", type=float, default=8192.0, help="R-TWT wake interval [us]")
    ap.add_argument("--p-err", type=float, default=0.22, help="first-attempt error probability")
    ap.add_argument("--q", type=float, default=0.69,
                    help="retransmission error probability of the extended model")
    ap.add_argument("--inv-lambda", type=float, default=16000.0,
                    help="mean packet inter-arrival time [us]")
    ap.add_argument("--K", type=int, default=20, help="buffer size [packets]")
    ap.add_argument("--pctl", type=float, default=99.0, help="delay percentile in panel (d)")
    ap.add_argument("--keep-stalls", action="store_true",
                    help="do not remove receiver-stall / warm-up packets from delay statistics")
    ap.add_argument("--data-dir", default=os.path.join(common.DATA_DIR, "exp2_sp_duration"))
    args = ap.parse_args()

    model_kw = dict(S=args.S, p_err=args.p_err, pctl=args.pctl,
                    inv_lambda=args.inv_lambda, K=args.K)
    T = args.wi_us
    Rs = (1, 3)

    # measurements
    meas = {R: [] for R in Rs}
    for key, path in common.find_traces(args.data_dir, r"exp2_WI\d+us_R(?P<R>\d+)_N(?P<N>\d+)\.csv"):
        if key["R"] in meas:
            m = common.measured_metrics(path, args.pctl, drop_stalls=not args.keep_stalls)
            meas[key["R"]].append((key["N"], m))
    for R in Rs:
        meas[R].sort(key=lambda t: t[0])
        print(f"R={R}: {len(meas[R])} runs, N = {[n for n, _ in meas[R]]}")
    Ns = sorted({n for s in meas.values() for n, _ in s})

    # model at the same (integer) N
    ana = {R: [common.model_metrics(T, N=N, R=R, **model_kw) for N in Ns] for R in Rs}
    ext = {R: [common.model_metrics(T, N=N, R=R, p_retry=args.q, **model_kw) for N in Ns]
           for R in Rs if R > 1}

    plot_kw = dict(
        xlabel=r"R-TWT $T_{\mathrm{SP}}$ [slots]",
        pctl=args.pctl,
        model={R: (Ns, ana[R]) for R in Rs},
        meas={R: ([n for n, _ in s], [m for _, m in s]) for R, s in meas.items()},
        extra=[(f"Extended ($q={args.q:g}$), $R={R}$", plot_style.R_COLORS[R], Ns, e)
               for R, e in ext.items()],
        xticks=Ns, model_markers=True)
    os.makedirs(common.FIG_DIR, exist_ok=True)
    base = os.path.join(common.FIG_DIR, "exp2_sp_duration")
    outs = plot_style.save(plot_style.two_column_figure(**plot_kw), base + "_2col")
    outs += plot_style.save(plot_style.one_column_figure(**plot_kw), base + "_1col")
    print("Figure -> " + ", ".join(os.path.relpath(o, common.ROOT) for o in outs))

    # table
    rows = []
    for R in Rs:
        rows += [["measured", R, N] + [m[k] for k in METRICS] for N, m in meas[R]]
        rows += [["model", R, N] + [m[k] for k in METRICS] for N, m in zip(Ns, ana[R])]
        if R in ext:
            rows += [[f"extended_model_q{args.q:g}", R, N] + [m[k] for k in METRICS]
                     for N, m in zip(Ns, ext[R])]
    os.makedirs(common.RES_DIR, exist_ok=True)
    common.write_metrics_csv(
        os.path.join(common.RES_DIR, "exp2_sp_duration_metrics.csv"),
        ["source", "R", "sp_slots", "mean_delay_ms", "jitter_ms", "loss_prob",
         f"p{args.pctl:g}_delay_ms"], rows)


if __name__ == "__main__":
    main()
