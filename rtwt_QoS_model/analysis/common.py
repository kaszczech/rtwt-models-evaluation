"""
Shared helpers: reading packet traces, filtering receiver stalls, computing the
measured metrics and evaluating the analytical model.

Trace format (one CSV per run, one row per received packet):
    seq,tx_ts_ns,rx_ts_ns,payload_bytes
seq is the sender sequence number; missing numbers are lost packets.
tx_ts_ns / rx_ts_ns are the send / receive timestamps in ns (synchronised clocks).
"""

import csv
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "model"))
import rtwt_model as model  # noqa: E402

DATA_DIR = os.path.join(ROOT, "data")
FIG_DIR = os.path.join(ROOT, "figures")
RES_DIR = os.path.join(ROOT, "results")

# ---------------------------------------------------------------------------
# Receiver stall / warm-up filter
# ---------------------------------------------------------------------------
# A receiver stall is a gap in reception > STALL_GAP_MS followed by a burst of
# >= STALL_MIN_BURST packets arriving < STALL_BURST_GAP_MS apart (packets read
# from a buffer at once), the first of them with a delay > STALL_MIN_DELAY_MS.
# Single long gaps (normal for Poisson traffic with losses) are NOT stalls.
# Warm-up: leading packets with a delay > WARMUP_DELAY_MS (a queue built up
# before the link started transmitting).
STALL_GAP_MS = 100.0
STALL_BURST_GAP_MS = 0.2
STALL_MIN_BURST = 3
STALL_MIN_DELAY_MS = 50.0
WARMUP_DELAY_MS = 50.0


def warmup_count(delays_ms):
    """Number of warm-up packets at the start of the trace."""
    k = 0
    while k < delays_ms.size and delays_ms[k] > WARMUP_DELAY_MS:
        k += 1
    return k


def stall_mask(rx_ns, delays_ms):
    """Boolean mask (True = packet from a stall or warm-up) and a list of
    events (gap_ms, n_packets); warm-up is reported as gap 0."""
    mask = np.zeros(delays_ms.size, dtype=bool)
    events = []
    k = warmup_count(delays_ms)
    if k:
        mask[:k] = True
        events.append((0.0, int(k)))
    if rx_ns.size < 2:
        return mask, events
    gap = np.diff(rx_ns) / 1e6
    for i in np.where(gap > STALL_GAP_MS)[0]:
        j = i + 1
        while j + 1 < delays_ms.size and gap[j] < STALL_BURST_GAP_MS:
            j += 1
        n = j - i
        if n >= STALL_MIN_BURST and delays_ms[i + 1] > STALL_MIN_DELAY_MS:
            mask[i + 1:j + 1] = True
            events.append((float(gap[i]), int(n)))
    return mask, events


# ---------------------------------------------------------------------------
# Measured metrics
# ---------------------------------------------------------------------------
def read_trace(path):
    """Return (seq, rx_ns, delay_ms) arrays in file (reception) order."""
    seqs, rxs, delays = [], [], []
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        cols = [h.strip().lower() for h in next(reader)]
        i_seq, i_tx, i_rx = cols.index("seq"), cols.index("tx_ts_ns"), cols.index("rx_ts_ns")
        for r in reader:
            try:
                s, tx, rx = int(r[i_seq]), int(r[i_tx]), int(r[i_rx])
            except (ValueError, IndexError):
                continue
            seqs.append(s); rxs.append(rx); delays.append(rx - tx)
    return (np.asarray(seqs, dtype=np.int64), np.asarray(rxs, dtype=np.int64),
            np.asarray(delays, dtype=np.float64) / 1e6)


def loss_indicator(seqs):
    """0/1 array over all sent sequence numbers (1 = lost)."""
    s = np.unique(seqs)
    lost = np.ones(int(s.max() - s.min() + 1), dtype=bool)
    lost[s - s.min()] = False
    return lost


def measured_metrics(path, pctl, drop_stalls=True, verbose=True):
    """Mean delay, jitter (std. dev. of delay), packet loss and delay percentile.
    Packet loss is computed from all packets; delay statistics exclude stall and
    warm-up packets when drop_stalls is True."""
    seqs, rx_ns, d_ms = read_trace(path)
    loss = float(loss_indicator(seqs).mean())
    if drop_stalls:
        mask, events = stall_mask(rx_ns, d_ms)
        if events and verbose:
            print(f"  [stalls] {os.path.basename(path)}: removed {int(mask.sum())} packets ("
                  + ", ".join(f"gap {g:.0f} ms / {n} pkt" for g, n in events) + ")")
        d_ms = d_ms[~mask]
    return {
        "mean_delay_ms": float(np.mean(d_ms)),
        "jitter_ms": float(np.std(d_ms, ddof=1)),
        "loss_prob": loss,
        "pctl_delay_ms": float(np.percentile(d_ms, pctl)),
    }


def loss_ci_block_bootstrap(seqs, block=100, n_boot=2000, level=95.0, seed=0):
    """Confidence interval of the packet loss ratio by block bootstrap.
    Losses come in bursts (they are not independent), so the trace is cut into
    blocks of `block` consecutive packets that are resampled with replacement."""
    x = loss_indicator(seqs).astype(float)
    nb = x.size // block
    blocks = x[:nb * block].reshape(nb, block).mean(axis=1)
    rng = np.random.default_rng(seed)
    boot = [blocks[rng.integers(0, nb, nb)].mean() for _ in range(n_boot)]
    a = (100.0 - level) / 2
    lo, hi = np.percentile(boot, [a, 100.0 - a])
    return float(x.mean()), float(lo), float(hi)


# ---------------------------------------------------------------------------
# Model metrics
# ---------------------------------------------------------------------------
def model_metrics(T_us, S, N, R, p_err, pctl, p_retry=None, inv_lambda=16000.0, K=20):
    p = model.Params(S=S, inv_lambda=inv_lambda, p_err=p_err, R=R, N=N, K=K,
                     p_retry=p_retry)
    m = model.metrics(p, T_us)
    return {
        "mean_delay_ms": m["mean_delay_us"] / 1000.0,
        "jitter_ms": m["jitter_us"] / 1000.0,
        "loss_prob": m["packet_loss"],
        "pctl_delay_ms": model.delay_percentile_us(p, T_us, pctl) / 1000.0,
    }


# ---------------------------------------------------------------------------
# Data files
# ---------------------------------------------------------------------------
def find_traces(directory, pattern):
    """List of (dict of named regex groups as int, path) for files matching pattern."""
    rx = re.compile(pattern)
    out = []
    for name in sorted(os.listdir(directory)):
        m = rx.fullmatch(name)
        if m:
            out.append(({k: int(v) for k, v in m.groupdict().items()},
                        os.path.join(directory, name)))
    return out


def write_metrics_csv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for row in rows:
            w.writerow([f"{v:.4f}" if isinstance(v, float) else v for v in row])
    print(f"Table  -> {os.path.relpath(path, ROOT)}")
