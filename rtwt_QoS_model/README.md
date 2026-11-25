# R-TWT on openwifi: measurements vs. analytical model

Measurement data and analysis code for the experimental validation of the
analytical model of dedicated Restricted Target Wake Time (R-TWT) for real-time
applications:

> A. Belogaev, X. Shen, C. Pan, X. Jiang, C. Blondia, J. Famaey,
> "Dedicated Restricted Target Wake Time for Real-Time Applications in Wi-Fi 7,"
> IEEE WCNC 2024.

The model is extended with correlated retransmission errors (see
[Extended error model](#extended-error-model)) and compared with measurements on
an openwifi testbed.

## Repository layout

```
rtwt-openwifi-validation/
├── README.md
├── requirements.txt
├── model/
│   └── rtwt_model.py                 analytical model (Markov chain) + error-model extension
├── analysis/
│   ├── common.py                     trace reading, stall filter, metrics, model evaluation
│   ├── plot_style.py                 IEEE figure style (2-column and 1-column layouts)
│   ├── plot_exp1_wake_interval.py    Experiment 1 figure and table
│   ├── plot_exp2_sp_duration.py      Experiment 2 figure and table
│   └── calibrate_q.py                calibration of q and validation with R = 2
├── data/
│   ├── exp1_wake_interval/           Experiment 1 packet traces (10 runs)
│   └── exp2_sp_duration/             Experiment 2 packet traces (13 runs)
├── figures/                          generated figures (PDF, SVG, PNG)
└── results/                          generated tables of measured and model metrics
```

## Quick start

Python 3.10+ is required.

```bash
pip install -r requirements.txt
python analysis/plot_exp1_wake_interval.py
python analysis/plot_exp2_sp_duration.py
python analysis/calibrate_q.py
```

Run the commands from the repository root. With the default arguments they
reproduce the figures and tables in `figures/` and `results/`. Run any script
with `--help` to list its parameters.

## Experiments

Common setup: one RTA station sends Poisson traffic (200-byte UDP payload,
mean inter-arrival time 1/λ = 16 ms) to the receiver over openwifi, HT MCS 2,
20 MHz. The station transmits only inside the R-TWT service period (SP). Each run
contains 10 000 sent packets.

| Experiment | Varied parameter | Fixed parameters | Runs |
|---|---|---|---|
| 1 | wake interval WI = 1.024, 2.048, 4.096, 8.192, 16.384 ms | N = 3 slots | R = 1 and R = 3 |
| 2 | SP duration T_SP = N = 1 … 6 slots | WI = 8.192 ms | R = 1 and R = 3, plus R = 2 at N = 3 |

R is the maximum number of transmission attempts (R − 1 retransmissions).
WI is a multiple of 512 µs (TSF-based scheduling).

The runs with N = 3, WI = 8.192 ms belong to both experiments: the files
`exp1_N3_R{1,3}_WI08192us.csv` and `exp2_WI8192us_R{1,3}_N3.csv` are identical.

### File naming

```
data/exp1_wake_interval/exp1_N<N>_R<R>_WI<WI in µs>us.csv
data/exp2_sp_duration/exp2_WI<WI in µs>us_R<R>_N<N>.csv
```

### Trace format

One CSV file per run, one row per **received** packet:

| Column | Description |
|---|---|
| `seq` | sender sequence number (0 … 9999); missing numbers are lost packets |
| `tx_ts_ns` | send timestamp at the sender [ns] |
| `rx_ts_ns` | receive timestamp at the receiver [ns] |
| `payload_bytes` | UDP payload size [bytes] |

The sender and receiver clocks are synchronised, so the packet delay is
`rx_ts_ns − tx_ts_ns`. The absolute timestamps are not wall-clock time (the
boards' clocks are reset at boot) and only differences are meaningful.

## Metrics

| Metric | Definition |
|---|---|
| Mean delay | mean of `rx_ts_ns − tx_ts_ns` over received packets |
| Jitter | standard deviation of the delay |
| Packet loss probability | 1 − (received unique packets / sent packets), from `seq` |
| 99th percentile delay | 99th percentile of the delay |

**Receiver-stall filter.** A few runs contain receiver stalls: a gap > 100 ms
in reception followed by a burst of ≥ 3 packets arriving < 0.2 ms apart (read
from a buffer at once), the first with a delay > 50 ms. Some runs also start
with a warm-up (leading packets with a delay > 50 ms). These packets are
excluded from the delay statistics (not from the loss). The scripts print every
removed event. Use `--keep-stalls` to disable the filter.

## Model parameters

| Parameter | Value | Meaning |
|---|---|---|
| S | 340 µs | slot = airtime of one transmission attempt |
| 1/λ | 16 ms | mean packet inter-arrival time |
| K | 20 | buffer size [packets] |
| p_err | 0.22 | error probability of the first attempt (measured loss at R = 1) |
| q | 0.69 | error probability of a retransmission (extended model) |

## Extended error model

The original model assumes independent attempts, so the packet loss ratio is
PLR(R) = p_err^R. The measurements show much less gain from retransmissions:
with p_err ≈ 0.22 the original model predicts about 1% loss for R = 3, but about
10% is measured.

The extended model (simplified Gilbert channel) uses p_err for the first
attempt and q for every retransmission:

```
PLR(R) = p_err · q^(R−1)
```

An error leaves the channel in a "bad" state, in which the next attempt also
fails with probability q. With `p_retry=None` the model in `model/rtwt_model.py`
reduces exactly to the original one.

**Calibration** (`calibrate_q.py`, Experiment 2, N = 3):
p_err = PLR(R=1) = 0.227 and q = √(PLR(R=3) / PLR(R=1)) = 0.695.
Values of q for the other N are 0.65–0.70 (mean 0.68).
The figures use the rounded values p_err = 0.22 and q = 0.69.

**Validation** with the R = 2 run, which is not used for calibration:

| | PLR (R = 2, N = 3) |
|---|---|
| measured | 14.81% (95% CI 13.95–15.67%) |
| extended model, 0.22 · 0.69 | 15.18% |
| extended model, calibrated 0.227 · 0.695 | 15.75% |
| original model, 0.227² | 5.14% |

Confidence intervals come from a block bootstrap (blocks of 100 packets, 2000
resamples). Losses come in bursts, so the binomial formula, which assumes
independent packets, would understate the interval.

**Known limitation.** The model gives every attempt its own slot. When N is
small or the WI is long, openwifi can retransmit within the same SP, so the
model overestimates the delay for R = 3 at N ≤ 3.

## Outputs

| File | Content |
|---|---|
| `figures/exp1_wake_interval_2col.*` | Experiment 1, 2×2 layout, 7.16 in wide (full page width) |
| `figures/exp1_wake_interval_1col.*` | Experiment 1, 4×1 layout, 3.5 in wide (one column) |
| `figures/exp2_sp_duration_2col.*` | Experiment 2, 2×2 layout |
| `figures/exp2_sp_duration_1col.*` | Experiment 2, 4×1 layout |
| `results/exp1_wake_interval_metrics.csv` | measured, model and extended-model metrics |
| `results/exp2_sp_duration_metrics.csv` | measured, model and extended-model metrics |

Include the figures at their native width (`\textwidth` or `\columnwidth`) to
keep the font sizes.
