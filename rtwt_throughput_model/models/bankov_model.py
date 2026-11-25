"""
Analytical model of paper 3:
  D. V. Bankov, A. I. Lyakhov, E. A. Stepanova, E. M. Khorov,
  "Performance Evaluation of Wi-Fi 7 Networks with Restricted Target Wake Time,"
  Problems of Information Transmission, vol. 60, no. 3, 2024.

MEASURED-INPUT VARIANT
----------------------
This version does NOT compute OFDM frame timings from K_sym / T_sym / preamble.
Instead you supply the airtimes you measure on your hardware:
  * the data-frame airtime as a function of the number of aggregated packets k,
  * the fixed control-frame airtimes (RTS, CTS, BAck).
The only thing that MUST stay k-dependent is the data-frame airtime -- that
dependence is the physical source of the non-monotonic throughput curve.
Provide it either as two numbers (linear model  T_data(k) = a + b*k) or as a
lookup table {k: airtime_us} you measured directly.

Computes the saturation throughput S of N non-RTA (user) stations under R-TWT,
where a periodic R-TWT moment every T_p may not be crossed by a transmission
(the last transmission in a period shrinks its A-MPDU to fit, or is postponed).

Equation numbers refer to the paper. Units: time us, throughput Mbit/s.
"""

from dataclasses import dataclass, field
import sys
from math import ceil, floor, lgamma, log, exp
sys.setrecursionlimit(100000)


@dataclass
class Params3:
    N: int = 10
    W_min: int = 16
    W_max: int = 1024
    R: int = 7
    L: int = 1500                 # packet payload (bytes) -- turns packets into bits
    T_e: float = 9.0              # empty backoff slot (us)
    T_SIFS: float = 16.0
    T_AIFS: float = 43.0
    T_lim: float = 1000.0         # TXOP limit (us)

    # ---- MEASURED data-frame airtime T_data(k) ----
    # Option A (default): linear model  T_data(k) = data_a + data_b * k
    #   data_a = fixed overhead (preamble + PHY/MAC headers), us
    #   data_b = airtime added per aggregated packet, us
    # Measure two points (e.g. k=1 and k=10) and solve:
    #   data_b = (airtime_k10 - airtime_k1) / 9
    #   data_a =  airtime_k1 - data_b
    data_a: float = 44.0          # <-- MEASURE AND REPLACE
    data_b: float = 95.2          # <-- MEASURE AND REPLACE
    # Option B: exact measured table {k: airtime_us}. Overrides the linear model
    # for the k values present; other k fall back to linear.
    data_table: dict = field(default_factory=dict)   # e.g. {1: 139.2, 2: 234.4}

    # ---- MEASURED control-frame airtimes (fixed sizes -> constants) ----
    T_RTS_meas: float = 32.0      # <-- MEASURE AND REPLACE
    T_CTS_meas: float = 28.0      # <-- MEASURE AND REPLACE
    T_BAck_meas: float = 36.0     # <-- MEASURE AND REPLACE

    @property
    def T_RTS(self):  return self.T_RTS_meas
    @property
    def T_CTS(self):  return self.T_CTS_meas
    @property
    def T_BAck(self): return self.T_BAck_meas

    @property
    def T_EIFS(self): return self.T_SIFS + self.T_BAck + self.T_AIFS
    @property
    def T_c(self):    return self.T_RTS + self.T_EIFS          # eq. (7)


# ---- data-frame airtime from MEASURED input ----
def T_data(p, k):
    """Airtime (us) of a data frame aggregating k packets, from measurement.
    Uses the exact table if k is present, else the linear model a + b*k.
    The k-dependence is essential -- it produces the non-monotonic throughput."""
    if k <= 0:
        return 0.0
    if k in p.data_table:
        return p.data_table[k]
    return p.data_a + p.data_b * k


def T_s_exchange(p, k):        # successful exchange WITHOUT the trailing AIFS
    return p.T_RTS + p.T_CTS + T_data(p, k) + p.T_BAck + 3 * p.T_SIFS

def T_s(p, T):                 # successful virtual slot (eq. 6, includes AIFS)
    return T_s_exchange(p, k_of_T(p, T)) + p.T_AIFS

def T_s_min(p):                # min exchange (k=1), no AIFS
    return T_s_exchange(p, 1)

def k_of_T(p, T, k_cap=10000):  # eq. (10): max k with exchange <= min(T, T_lim)
    budget = min(T, p.T_lim)
    if T_s_exchange(p, 1) > budget:
        return 0
    k = 1
    while k < k_cap and T_s_exchange(p, k + 1) <= budget:
        k += 1
    return k


# ---- Bianchi/Vishnevsky saturation (eqs. 2-5) ----
def bianchi(p):
    W = [min((2 ** r) * p.W_min, p.W_max) for r in range(p.R + 1)]
    def tau_of(pc):
        num = sum(pc ** r for r in range(p.R + 1))
        den = sum(((W[r] + 1) / 2.0) * pc ** r for r in range(p.R + 1))
        return num / den
    pc = 0.0
    for _ in range(100000):
        tau = tau_of(pc)
        pc_new = 1 - (1 - tau) ** (p.N - 1)
        if abs(pc_new - pc) < 1e-14:
            pc = pc_new; break
        pc = pc_new
    tau = tau_of(pc)
    pi_e = (1 - tau) ** p.N
    pi_s = p.N * tau * (1 - tau) ** (p.N - 1)
    pi_c = 1 - pi_e - pi_s
    return pi_e, pi_s, pi_c




def _make_P_tail(p, pi_e, pi_c):
    """P_tail(T*) via DP: prob. that a slot sequence filling T* has no success.
    f(T) = pi_e*f(T-Te) + pi_c*f(T-Tc) for T>=T_s^min, else 1  (== eq. 13)."""
    Tmin = T_s_min(p)
    memo = {}
    def f(T):
        if T < Tmin - 1e-9:
            return 1.0
        key = round(T, 3)
        v = memo.get(key)
        if v is not None:
            return v
        v = pi_e * f(T - p.T_e) + pi_c * f(T - p.T_c)
        memo[key] = v
        return v
    return f


def _mean_payload_and_shrink(p, T_p_eff, pi_e, pi_s, pi_c):
    Ts_lim = T_s(p, p.T_lim)
    Tmin = T_s_min(p)
    k_lim = k_of_T(p, p.T_lim)
    f = _make_P_tail(p, pi_e, pi_c)
    le = log(pi_e) if pi_e > 0 else -1e9
    ls = log(pi_s) if pi_s > 0 else -1e9
    lc = log(pi_c) if pi_c > 0 else -1e9
    L_bar = 0.0
    T_bar = 0.0
    s_max = int((T_p_eff - Tmin) // Ts_lim)
    C_CAP = 25          # pi_c**c negligible beyond this
    for s in range(0, s_max + 1):
        base_s = T_p_eff - Ts_lim * s
        if base_s - Tmin < 0:
            break
        for c in range(0, C_CAP + 1):
            base_sc = base_s - p.T_c * c
            if base_sc - Tmin < 0:      # even e=0 too tight -> stop increasing c
                break
            e = 0
            while True:
                T_tail = base_sc - p.T_e * e
                if T_tail < Tmin:       # increasing e only lowers T_tail
                    break
                t = s + c + e
                logw = (lgamma(t + 1) - lgamma(s + 1) - lgamma(c + 1) - lgamma(e + 1)
                        + e * le + s * ls + c * lc)
                w = pi_s * exp(logw)
                Ts_tail = T_s(p, T_tail)
                ptail = 1.0 if T_tail <= Ts_tail else f(T_tail - Ts_tail)
                w *= ptail
                if w > 1e-18:
                    k_tail = k_of_T(p, T_tail)
                    L_bar += w * (k_lim * s + k_tail) * p.L
                    T_bar += w * ((Ts_tail - T_tail) if T_tail < Ts_tail else p.T_e / 2.0)
                e += 1
    return L_bar, T_bar


def throughput(p, T_p, iters=1, PER=0.0):
    """Non-RTA throughput S (Mbit/s) at R-TWT period T_p (us).
    `iters`>0 applies the effective-period correction T_p* = T_p - T_bar (eq. 15);
    the correction is ~Te/2, so iters=1 is plenty."""
    pi_e, pi_s, pi_c = bianchi(p)
    T_p_eff = T_p
    for _ in range(iters):
        L_bar, T_bar = _mean_payload_and_shrink(p, T_p_eff, pi_e, pi_s, pi_c)
        T_p_eff = T_p - T_bar
    L_bar, _ = _mean_payload_and_shrink(p, T_p_eff, pi_e, pi_s, pi_c)
    return 8 * L_bar / T_p * (1 - PER)          # eq. (11)


def throughput_no_rtwt(p, PER=0.0):
    """Reference saturation throughput without R-TWT (the horizontal line)."""
    pi_e, pi_s, pi_c = bianchi(p)
    k_lim = k_of_T(p, p.T_lim)
    num = 8 * pi_s * k_lim * p.L
    den = pi_e * p.T_e + pi_s * T_s(p, p.T_lim) + pi_c * p.T_c
    return num / den * (1 - PER)


if __name__ == "__main__":
    # Defaults reproduce the paper's MCS5 case via the LINEAR model
    # (data_a=44, data_b=95.2 fitted to the paper's T_data(k)).
    p = Params3()
    print("Using MEASURED-input model.")
    print("  T_data(1)=%.1f  T_data(5)=%.1f  T_data(10)=%.1f us"
          % (T_data(p, 1), T_data(p, 5), T_data(p, 10)))
    print("  T_RTS=%.1f  T_CTS=%.1f  T_BAck=%.1f  T_c=%.1f"
          % (p.T_RTS, p.T_CTS, p.T_BAck, p.T_c))
    print("  T_s^min=%.1f us" % T_s_min(p))
    print("  reference S without R-TWT = %.1f Mbit/s\n" % throughput_no_rtwt(p))
    print("S vs R-TWT period T_p:")
    print("%9s %9s" % ("T_p [ms]", "S [Mbps]"))
    for tp_ms in (0.3, 0.5, 0.7, 0.9, 0.95, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0):
        print("%9s %9.1f" % (tp_ms, throughput(p, tp_ms * 1000)))

    print("\nSame but with a measured T_data lookup table (illustrative):")
    p2 = Params3(data_table={1: 139.2, 2: 234.4, 3: 329.6})
    for tp_ms in (0.5, 1.0, 2.0):
        print("  T_p=%sms -> S=%.1f Mbps" % (tp_ms, throughput(p2, tp_ms * 1000)))
