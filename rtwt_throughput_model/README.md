# R-TWT Throughput Model

Analytical throughput models for Restricted Target Wake Time (R-TWT) under
saturated traffic, compared against openwifi hardware measurements for two
stations. This repository contains everything needed to reproduce the
aggregate-throughput figure `S(T_V)`.

## Contents

```
rtwt_throughput_model/
├── models/
│   ├── bankov_model.py          Original Bankov R-TWT model (EDCA, RTS/CTS,
│   │                            A-MPDU structure, classic Bianchi, ideal channel)
│   └── dcf_model.py             Adapted model for openwifi (DCF, no RTS/CTS,
│                                fixed frame, error-prone channel)
├── data/
│   ├── throughput_vs_period.csv Clean 2-column measurement used by the figure
│   └── throughput_2sta.csv      Full measurement export (all metrics)
├── figure_throughput_vs_period.py   Builds the comparison figure
├── requirements.txt
└── README.md
```

## Quick start

```bash
pip install -r requirements.txt
python figure_throughput_vs_period.py
```

This writes `figures/throughput_vs_period.{pdf,svg,png}` (IEEE single-column,
Times New Roman). Each model file also runs a self-test when executed directly:

```bash
python models/dcf_model.py
python models/bankov_model.py
```

## The two models

| Aspect | `bankov_model.py` (original) | `dcf_model.py` (adapted) |
| --- | --- | --- |
| Channel access | EDCA (AIFS) | DCF (DIFS) |
| RTS/CTS | yes (short collision) | no (collision = full frame) |
| Frame | A-MPDU structure | fixed-size frame |
| Collision prob. | `1 - (1 - tau)^(N-1)` | `1 - (1 - tau)^(N-1)*(1 - p_err)` |
| Channel | ideal | error-prone (`p_err`) |
| Retry / CW | R = 7, CW 16/1024 | R = 3, CW 15/1023 |

Both compute the saturation throughput over an R-TWT period with a dynamic
program that only counts transmissions finishing before the period boundary.

Note: the Bankov paper writes the conditional collision probability as
`p = 1 - (1 - tau)^N`. This code uses the classic Bianchi form
`p = 1 - (1 - tau)^(N-1)` (N-1 other stations); for N = 2 the difference in
throughput is below 1%.

## Measurement data

Two saturated stations, uplink UDP (ports 9001 / 9002), 1000 B payload,
MCS 2 (802.11n HT, 20 MHz, long GI).

`data/throughput_vs_period.csv` is the clean file read by the figure:

| Column | Meaning |
| --- | --- |
| `period_us` | Effective R-TWT window `T_V` [us] (scheduled period plus one-frame overrun) |
| `throughput_Mbps` | Aggregate windowed throughput of frames finishing inside the window |

`data/throughput_2sta.csv` is the full export it was extracted from; its
`throughput_Mbps` column is the `filt_win_Mbps` column below:

| Column | Meaning |
| --- | --- |
| `W_us` | Same as `period_us` above |
| `meas_sum_Mbps` | Aggregate throughput over the whole period, all frames |
| `meas_win_Mbps` | Aggregate throughput normalised to the window, all frames |
| `filt_sum_Mbps` | Over the whole period, only frames finishing inside the window |
| `filt_win_Mbps` | Normalised to the window, only frames finishing inside it (used in the figure) |

## Physical parameters (MCS 2, 802.11n HT, 20 MHz, long GI)

| Quantity | Value |
| --- | --- |
| Nominal rate | 19.5 Mbit/s |
| Data-frame airtime | 484 us |
| One-frame channel time (`T_s`, no RTS/CTS) | 578 us |
| Slot time | 9 us |
| SIFS | 16 us |
| DIFS | 34 us |
