# Experimental Evaluation of Analytical Models for R-TWT in Wi-Fi 7

Measurement data and code for the paper:

> K. Szczech, P. Avila-Campos, J. Haxhibeqiri, J. Hoebeke, S. Szott,
> "Experimental Evaluation of Analytical Models for Restricted Target Wake Time in Wi-Fi 7,"
> submitted to *Computer Networks*.

The repository compares two analytical models of Restricted Target Wake Time
(R-TWT) with measurements from an openwifi SDR testbed.

| Folder | Model | Paper |
|---|---|---|
| [`rtwt_QoS_model/`](rtwt_QoS_model) | QoS of RTA traffic (Belogaev et al., WCNC 2024), extended with a retransmission failure probability *q* | Figs. 2–3, Table 3 |
| [`rtwt_throughput_model/`](rtwt_throughput_model) | Non-RTA throughput (Bankov et al., 2024), original and adapted to openwifi DCF | Fig. 4 |

Each folder contains the model code, the raw measurements, scripts that
reproduce the figures, and its own README with details.

## Quick start

Python 3.10+.

```bash
cd rtwt_QoS_model
pip install -r requirements.txt
python analysis/plot_exp1_wake_interval.py   # Fig. 2
python analysis/plot_exp2_sp_duration.py     # Fig. 3
python analysis/calibrate_q.py               # Table 3 (calibration of q)

cd ../rtwt_throughput_model
pip install -r requirements.txt
python figure_throughput_vs_period.py        # Fig. 4
```

## Testbed

openwifi (Xilinx ZC702 + AD FMCOMMS2/3/4), IEEE 802.11n, MCS 2, 20 MHz,
PTP-synchronized nodes, imec Industry 4.0 testbed (Ghent, Belgium).

## Citation

If you use this data or code, please cite the paper above
(citation details will be updated after publication).

## License

Code: MIT. Data: CC BY 4.0.
