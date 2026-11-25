"""
Analytical model of a dedicated R-TWT service for one real-time (RTA) flow.

Based on:
  A. Belogaev, X. Shen, C. Pan, X. Jiang, C. Blondia, J. Famaey,
  "Dedicated Restricted Target Wake Time for Real-Time Applications in Wi-Fi 7,"
  IEEE WCNC 2024.

A single RTA flow uses dedicated, contention-free R-TWT service periods (SPs).
The continuous system is approximated by a slotted batch-arrival queue (a packet
needing r transmissions == a batch of size r). A 2-D Markov chain over
(queue length k, slot index n) gives the stationary distribution, from which the
packet delay distribution, mean delay, jitter, packet loss and delay percentiles
follow. Equation numbers refer to the paper. Delay is computed in slots
(multiply by S for microseconds).

Extension (this work): correlated retransmission errors. The first attempt fails
with probability p_err, every retransmission (attempts 2..R) fails with
probability p_retry = q. This is a simplified Gilbert channel: an error moves the
channel to a "bad" state in which the next attempt fails with probability q.
With p_retry=None the attempts are independent, exactly as in the paper.

Requires numpy (scipy is used if available for a faster sparse solve).
"""

from dataclasses import dataclass
from math import exp, sqrt
import numpy as np

try:
    import scipy.sparse as _sp
    import scipy.sparse.linalg as _spla
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


@dataclass
class Params:
    S: float = 114.4              # slot = airtime of one transmission attempt (us)
    inv_lambda: float = 16000.0   # mean packet inter-arrival time 1/lambda (us)
    p_err: float = 0.10           # error probability of the first attempt
    R: int = 3                    # max number of transmission attempts
    K: int = 20                   # buffer size (packets)
    N: int = 3                    # SP duration in slots (SP = N*S)
    p_retry: float = None         # error probability of attempts 2..R (q);
                                  # None -> p_err (independent attempts, as in the paper)

    def attempt_err(self, r):
        """Error probability of attempt r (1-based)."""
        if r == 1 or self.p_retry is None:
            return self.p_err
        return self.p_retry

    def fail_prob(self):
        """Probability that all R attempts fail (packet loss ratio)."""
        out = 1.0
        for r in range(1, self.R + 1):
            out *= self.attempt_err(r)
        return out

    def M_slots(self, T_us):
        """Vacation length in slots for a period (wake interval) T: T = (N+M)*S."""
        total = round(T_us / self.S)
        return max(0, total - self.N)


# ---------------------------------------------------------------------------
# Batch-arrival probabilities (Section III-A, eqs. 1-2)
# ---------------------------------------------------------------------------
def batch_probs(p: Params):
    lamS = p.S / p.inv_lambda                 # lambda * S
    b0 = exp(-lamS)                           # P(no packet arrives in a slot)
    b = 1.0 - b0                              # P(a packet arrives)
    # P(exactly r attempts, the last one successful) = prod_{i<r} e_i * (1 - e_r),
    # e_i = attempt_err(i). With p_retry=None this is (1-p_err) p_err^(r-1).
    # reach[r] = P(attempt r takes place at all)
    reach = [1.0] * (p.R + 2)
    for r in range(1, p.R + 1):
        reach[r + 1] = reach[r] * p.attempt_err(r)
    # b_hat[r] = P(a batch of size r arrives), regardless of success (eq. 2)
    b_hat = [0.0] * (p.R + 1)
    b_hat[0] = b0
    for r in range(1, p.R):
        b_hat[r] = b * reach[r] * (1 - p.attempt_err(r))
    b_hat[p.R] = b * reach[p.R]
    # successful-batch probabilities b_r and loss probability b_fail (eq. 1)
    b_succ = [0.0] * (p.R + 1)
    for r in range(1, p.R + 1):
        b_succ[r] = b * reach[r] * (1 - p.attempt_err(r))
    b_fail = b * reach[p.R + 1]
    return b0, b, b_hat, b_succ, b_fail


# ---------------------------------------------------------------------------
# Markov chain transitions and stationary distribution (Section III-B, eq. 3)
# ---------------------------------------------------------------------------
def stationary(p: Params, T_us):
    M = p.M_slots(T_us)
    N, K, R = p.N, p.K, p.R
    n_slots = N + M                                  # slots per period
    b0, b, b_hat, b_succ, b_fail = batch_probs(p)

    def idx(k, n):
        return k * n_slots + n

    n_states = (K + 1) * n_slots
    rows, cols, vals = [], [], []

    def add(k, n, k2, n2, prob):
        if prob <= 0:
            return
        rows.append(idx(k, n)); cols.append(idx(k2, n2)); vals.append(prob)

    for n in range(n_slots):
        n_next = (n + 1) % n_slots
        in_SP = (n < N)
        for k in range(K + 1):
            if not in_SP:
                # --- vacation ---
                # (a) no batch, or batch too big to fit -> k unchanged
                p_stay = b_hat[0] + sum(b_hat[r] for r in range(K - k + 1, R + 1))
                add(k, n, k, n_next, p_stay)
                # (b) batch of size r fits (1<=r<=K-k) -> k+r
                for r in range(1, min(R, K - k) + 1):
                    add(k, n, k + r, n_next, b_hat[r])
            else:
                # --- SP: one attempt is transmitted (k -> k-1) unless k==0 ---
                if k == 0:
                    # stay at 0 (nothing to transmit) unless a batch arrives and fits
                    p_stay = b_hat[0] + sum(b_hat[r] for r in range(K + 1, R + 1))
                    add(0, n, 0, n_next, p_stay)
                    for r in range(1, min(R, K) + 1):        # batch at empty queue: -> r-1
                        add(0, n, r - 1, n_next, b_hat[r])
                else:
                    # no batch / batch too big -> k-1
                    p_down = b_hat[0] + sum(b_hat[r] for r in range(K - k + 1, R + 1))
                    add(k, n, k - 1, n_next, p_down)
                    # batch of size r fits -> k+r-1 (one attempt transmitted this slot)
                    for r in range(1, min(R, K - k) + 1):
                        add(k, n, k + r - 1, n_next, b_hat[r])

    # solve pi P = pi, sum(pi) = 1
    if _HAVE_SCIPY:
        P = _sp.csr_matrix((vals, (rows, cols)), shape=(n_states, n_states))
        A = (P.T - _sp.eye(n_states)).tolil()
        A[-1, :] = 1.0
        rhs = np.zeros(n_states); rhs[-1] = 1.0
        pi = _spla.spsolve(A.tocsr(), rhs)
    else:
        P = np.zeros((n_states, n_states))
        for r_, c_, v_ in zip(rows, cols, vals):
            P[r_, c_] += v_
        A = (P.T - np.eye(n_states))
        A[-1, :] = 1.0
        rhs = np.zeros(n_states); rhs[-1] = 1.0
        pi = np.linalg.solve(A, rhs)

    pi = np.maximum(pi, 0.0)
    pi /= pi.sum()
    return pi.reshape(K + 1, n_slots), M


# ---------------------------------------------------------------------------
# Packet delay distribution (Section III-C, eqs. 4-9)
# ---------------------------------------------------------------------------
def _D_V(q, N, M):
    """Vacation slots incurred by a queue of size q (eq. 5)."""
    if q <= 0:
        return 0
    if q % N != 0:
        return (q // N) * M
    return (q // N - 1) * M


def delay_distribution(p: Params, T_us):
    """Return (delays_slots, probs, S_us, meta) for successfully delivered packets."""
    pi, M = stationary(p, T_us)
    N, K, R = p.N, p.K, p.R
    n_slots = N + M
    b0, b, b_hat, b_succ, b_fail = batch_probs(p)
    denom = b * (1.0 - b_fail)

    dist = {}   # delay in slots -> probability
    for n in range(n_slots):
        in_SP = (n < N)
        for k in range(K + 1):
            pk = pi[k, n]
            if pk <= 0:
                continue
            for r in range(1, R + 1):
                p_tilde = pk * b_succ[r] / denom          # eq. (4)
                if p_tilde <= 0:
                    continue
                if in_SP:                                  # eq. (7)
                    q_tilde = k + r - min(N - n, k + r)
                    H = 1 if q_tilde > 0 else 0
                    D = (k + r) + H + _D_V(q_tilde, N, M)
                else:                                      # eq. (6)
                    D = (n_slots - n) + (k + r) + _D_V(k + r, N, M)
                dist[D] = dist.get(D, 0.0) + p_tilde

    tot = sum(dist.values())
    if tot > 0:
        for d in dist:
            dist[d] /= tot
    delays = np.array(sorted(dist))
    probs = np.array([dist[d] for d in delays])
    return delays, probs, p.S, {"M": M, "b_fail": b_fail, "b": b}


# ---------------------------------------------------------------------------
# Performance indicators
# ---------------------------------------------------------------------------
def delay_percentile_us(p: Params, T_us, pctl):
    """Delay percentile pctl (in %, e.g. 99 or 99.9) in microseconds:
    the smallest delay d with CDF(d) >= pctl/100."""
    delays, probs, S, _ = delay_distribution(p, T_us)
    d_us = delays * S
    idx = int(np.searchsorted(np.cumsum(probs), pctl / 100.0))
    return float(d_us[min(idx, len(d_us) - 1)])


def metrics(p: Params, T_us):
    """Mean delay, jitter and packet loss (times in microseconds)."""
    delays, probs, S, meta = delay_distribution(p, T_us)
    d_us = delays * S
    mean = float(np.sum(d_us * probs))
    var = float(np.sum((d_us - mean) ** 2 * probs))
    return {
        "T_us": (p.N + meta["M"]) * p.S,     # effective period
        "mean_delay_us": mean,
        "jitter_us": sqrt(var),
        "packet_loss": p.fail_prob(),        # all R attempts fail (buffer overflow negligible)
    }


if __name__ == "__main__":
    # Sanity check against the paper (Table II: 1/lambda=16 ms, K=20, p_err=10%, S=114.4 us):
    # the 99.9th percentile should be ~10 ms at T=6 ms for N=3, R=3.
    print(f"{'R':>2} {'T[ms]':>6} {'mean[ms]':>9} {'jitter[ms]':>11} {'PLR':>10} {'P99.9[ms]':>10}")
    for R in (1, 3):
        for T_ms in (2, 4, 6, 8, 10, 12, 16):
            pp = Params(R=R, N=3)
            m = metrics(pp, T_ms * 1000)
            print(f"{R:>2} {m['T_us']/1000:>6.1f} {m['mean_delay_us']/1000:>9.3f} "
                  f"{m['jitter_us']/1000:>11.3f} {m['packet_loss']:>10.1e} "
                  f"{delay_percentile_us(pp, T_ms * 1000, 99.9)/1000:>10.3f}")
