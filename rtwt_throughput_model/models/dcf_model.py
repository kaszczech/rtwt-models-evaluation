"""
Adapted DCF throughput model for R-TWT on openwifi.

Independent implementation (does NOT import bankov_model). It adapts the Bankov
R-TWT throughput model to the openwifi hardware reality:

  1. No RTS/CTS -> a collision occupies the WHOLE data frame: T_c ~ T_s
     (full transmission + EIFS), not a short RTS. Bankov's eq. (7) does not hold.
  2. No aggregation -> fixed-size frame (payload bytes). The A-MPDU shrinking
     case disappears: T_s^min = T_s^max, so the model reduces to "fits / does
     not fit" within the period.
  3. DCF instead of EDCA: DIFS instead of AIFS; CWmin/CWmax and R from openwifi.
  4. Conditional failure probability (classic Bianchi, N-1 other stations) with
     an error-prone channel:  p = 1 - (1 - tau)^(N-1) * (1 - p_err).
  5. A slot is successful only if exactly one station transmits AND no error.

Throughput over an R-TWT period W is computed with a dynamic program: the
expected number of successful transmissions that FINISH before the period
boundary (the "do not cross T" rule). The effective window W = T + overrun is
supplied by the caller. Units: time in us, throughput in Mbit/s.
"""

from dataclasses import dataclass
import sys


@dataclass
class ParamsDCF:
    N: int = 2                 # number of saturated stations
    aCWmin: int = 15           # CWmin (window size = aCWmin + 1 = 16)
    aCWmax: int = 1023         # CWmax (window size = 1024)
    R: int = 3                 # transmission attempts per frame (openwifi)
    payload: int = 1000        # payload bytes (fixed-size frame)
    p_err: float = 0.0         # frame error probability (non-ideal channel)
    # durations [us]
    T_data: float = 484.0      # data-frame airtime (MCS2, 1086 B PSDU)
    T_ACK: float = 44.0        # ACK frame
    T_SIFS: float = 16.0
    T_slot: float = 9.0        # empty backoff slot

    @property
    def T_DIFS(self):          # DCF: DIFS = SIFS + 2*slot
        return self.T_SIFS + 2 * self.T_slot          # 34 us

    @property
    def T_EIFS(self):
        return self.T_SIFS + self.T_ACK + self.T_DIFS

    @property
    def T_s(self):             # successful slot: DIFS + DATA + SIFS + ACK
        return self.T_DIFS + self.T_data + self.T_SIFS + self.T_ACK

    @property
    def T_c(self):             # collision/error: full frame + EIFS (no valid ACK)
        return self.T_data + self.T_EIFS


def bianchi(prm):
    """Return (tau, p, P_e, P_s, P_f) for saturated DCF with an error channel."""
    W0 = prm.aCWmin + 1
    Wmax = prm.aCWmax + 1
    W = [min((2 ** r) * W0, Wmax) for r in range(prm.R + 1)]

    def tau_of(p):
        num = sum(p ** r for r in range(prm.R + 1))
        den = sum(((W[r] + 1) / 2.0) * (p ** r) for r in range(prm.R + 1))
        return num / den

    p = 0.0
    for _ in range(100000):
        tau = tau_of(p)
        # classic Bianchi: N-1 other stations, combined with the error channel
        p_new = 1.0 - (1.0 - tau) ** (prm.N - 1) * (1.0 - prm.p_err)
        if abs(p_new - p) < 1e-14:
            p = p_new
            break
        p = p_new
    tau = tau_of(p)
    P_e = (1.0 - tau) ** prm.N
    P_s = prm.N * tau * (1.0 - tau) ** (prm.N - 1) * (1.0 - prm.p_err)
    P_f = 1.0 - P_e - P_s
    return tau, p, P_e, P_s, P_f


def throughput_no_rtwt(prm):
    """Saturation throughput without R-TWT [Mbit/s]."""
    _, _, P_e, P_s, P_f = bianchi(prm)
    num = 8 * prm.payload * P_s
    den = P_e * prm.T_slot + P_s * prm.T_s + P_f * prm.T_c
    return num / den


def expected_successes(prm, W):
    """E[number of successful transmissions in a window W], with the boundary rule.

    In each virtual slot: empty (T_slot), success (T_s) or failure/collision
    (T_c). A transmission is only started if it can FINISH before the boundary
    (remaining time >= its duration)."""
    _, _, P_e, P_s, P_f = bianchi(prm)
    Te, Ts, Tc = prm.T_slot, prm.T_s, prm.T_c

    sys.setrecursionlimit(1000000)
    memo = {}

    def f(t):
        if t < Te:                     # nothing more fits
            return 0.0
        key = round(t, 2)
        v = memo.get(key)
        if v is not None:
            return v
        res = P_e * f(t - Te)          # empty slot always fits (t >= Te)
        if t >= Ts:                    # a success is only possible if it fits
            res += P_s * (1.0 + f(t - Ts))
        if t >= Tc:                    # a collision/error must also fit
            res += P_f * f(t - Tc)
        memo[key] = res
        return res

    return f(float(W))


def throughput(prm, W):
    """Throughput [Mbit/s] over an R-TWT window of length W [us]."""
    ns = expected_successes(prm, W)
    return 8 * prm.payload * ns / W


if __name__ == "__main__":
    print("=== dcf_model self-test (N=2, MCS2, no RTS/CTS, DCF) ===")
    for perr in (0.0, 0.10):
        p = ParamsDCF(N=2, p_err=perr)
        tau, pc, Pe, Ps, Pf = bianchi(p)
        print("\np_err=%s: tau=%.4f p=%.4f | P_e=%.3f P_s=%.3f P_f=%.3f"
              % (perr, tau, pc, Pe, Ps, Pf))
        print("  T_s=%.0f T_c=%.0f DIFS=%.0f | no R-TWT S=%.2f Mbps"
              % (p.T_s, p.T_c, p.T_DIFS, throughput_no_rtwt(p)))
        print("  %7s %9s" % ("W[us]", "S[Mbps]"))
        for W in (1200, 1600, 2000, 2400, 3000, 4000, 4600):
            print("  %7d %9.2f" % (W, throughput(p, W)))
