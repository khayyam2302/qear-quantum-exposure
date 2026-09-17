"""QEAR analytic resource model used for the revised manuscript.

All quantities are closed-form functions of the word size b and the message
equations of Table XI. Run `python3 qear_model.py` to write numbers.json.
"""
import math, json
from itertools import combinations
import numpy as np
from scipy.stats import kendalltau

# ---------------------------------------------------------------- cost rules
def swaps(b):
    """Controlled swaps of one in-place barrel rotation (stage j shifts by 2^j < b)."""
    return sum(b - math.gcd(b, 1 << j) for j in range(math.ceil(math.log2(b))) if (1 << j) < b)

def op_cost(op, b):
    """(X, CNOT, Toffoli, carry ancillas, gate depth, Toffoli depth)."""
    s = swaps(b)
    return {
        "XOR":      (0,     b,       0,     0, 1,         0),
        "COPY":     (0,     b,       0,     0, 1,         0),
        "AND":      (0,     0,       b,     0, 1,         1),
        "OR":       (5 * b, 0,       b,     0, 3,         1),
        "ADD":      (0,     4 * b,   2 * b, 1, 5 * b + 1, 2 * b),
        "ROT":      (0,     4 * s + b, 2 * s, 0, 6 * s + 1, 2 * s),   # rotate, XOR out, rotate back
        "ROT_NOUNDO": (0,   2 * s + b, s,   0, 3 * s + 1, s),         # no rotate-back (incorrect)
        "ROT_STATIC": (0,   b,       0,     0, 1,         0),         # relabelling
    }[op]

PROTOCOLS = {   # inputs, key registers, secret nonces, messages
    "LMAP": (7, 4, 2, {"A": ["XOR"] * 3, "B": ["OR", "ADD"], "C": ["COPY", "ADD", "ADD"], "D": ["COPY", "ADD", "XOR"]}),
    "EMAP": (7, 4, 2, {"A": ["XOR"] * 3, "B": ["OR", "XOR"], "C": ["AND", "XOR"], "D": ["OR", "XOR"]}),
    "M2AP": (6, 3, 2, {"A": ["XOR"] * 3, "B": ["OR", "ADD"], "C": ["COPY", "ADD", "ADD"]}),
    "SASI": (5, 2, 2, {"A": ["XOR"] * 3, "B": ["OR", "ADD"], "C": ["ROT", "ROT", "ADD"]}),
    "URAP": (5, 2, 2, {"M1": ["ROT", "XOR", "XOR"], "M2": ["ROT", "XOR", "XOR"],
                       "M3": ["ROT", "XOR", "XOR"], "M4": ["ROT", "XOR", "XOR"]}),
}
ORDER = ["LMAP", "EMAP", "M2AP", "SASI", "URAP", "Ascon-Hash256"]
ASCON_CALLS = 8            # one Ascon-Hash256 on a 256-bit input with 256-bit output

def ulw_circuit(name, b=32, rot="ROT"):
    n_in, n_key, n_nonce, msgs = PROTOCOLS[name]
    X = CX = TOF = anc = 0; D = []; TD = []; nl = 0
    for ops in msgs.values():
        d = td = 0
        for op in ops:
            op = rot if op == "ROT" else op
            x, c, t, a, gd, tdep = op_cost(op, b)
            X += x; CX += c; TOF += t; anc += a; d += gd; td += tdep
            nl += op in ("AND", "OR", "ADD", "ROT", "ROT_NOUNDO")
        D.append(d); TD.append(td)
    Q = n_in * b; A = len(msgs) * b + anc
    return dict(Q=Q, A=A, W=Q + A, X=X, CNOT=CX, Toffoli=TOF, G=X + CX + TOF, D=max(D),
                QC=max(D) * (Q + A), T=7 * TOF, TD=3 * max(TD), TD_serial=3 * sum(TD),
                k=n_key * b, m=len(msgs) * b, m_eff=(len(msgs) - n_nonce) * b, ops=sum(map(len, msgs.values())),
                nonlinear=nl)

def ascon_round():
    """Bitsliced S-box with compute-uncompute of the five product ancillas (reused each round)."""
    slices = 64
    return dict(X=slices * 21 + 8, CNOT=slices * 11 + 640, Toffoli=slices * 10, D=45, ToffD=10)

def ascon_circuit(calls=ASCON_CALLS):
    r = ascon_round(); rounds = 12 * calls
    Q, A = 320 + 256 + 256, 320
    X, CX, TOF = r["X"] * rounds, r["CNOT"] * rounds, r["Toffoli"] * rounds
    D, TD = r["D"] * rounds, 3 * r["ToffD"] * rounds
    return dict(Q=Q, A=A, W=Q + A, X=X, CNOT=CX, Toffoli=TOF, G=X + CX + TOF, D=D, QC=D * (Q + A),
                T=7 * TOF, TD=TD, TD_serial=TD, k=128, m=256, m_eff=256, ops=3, nonlinear=calls)

def attack(c, td_key="TD"):
    k, m = c["k"], c["m_eff"]; r = math.ceil(k / m)
    n = math.floor(math.pi / 4 * 2 ** (k / 2))
    t_or = 2 * r * c["T"] + 14 * (r * m - 1)
    per = t_or + 14 * (k - 1)
    td_it = 2 * c[td_key] + 3 * (r * m - 1) + 3 * (k - 1)
    return dict(r=r, log2_Niter=math.log2(n), T_oracle=t_or, per_iter=per,
                log2_Ttotal=math.log2(n * per), log2_TDtotal=math.log2(n * td_it),
                depth_iter=td_it, W_oracle=k + r * c["W"])

def nist(log2T, log2TD, md):
    cost = log2T + max(0.0, log2TD - md)
    return cost, cost - (170 - md)

def machines(log2TD, md):
    return max(0.0, 2 * (log2TD - md))

def surface(Wor, log2TD, p=1e-3, pth=1e-2, tc=1e-6):
    for d in range(3, 999, 2):
        lg = math.log10(0.1) + (d + 1) / 2 * math.log10(p / pth) + math.log10(Wor) + log2TD * math.log10(2) + math.log10(d)
        if lg < 0:
            return d, 2 * d * d * Wor, 2 ** log2TD * d * tc / 3.156e7

def circuits(b=32, rot="ROT"):
    c = {n: ulw_circuit(n, b, rot) for n in PROTOCOLS}
    c["Ascon-Hash256"] = ascon_circuit()
    return c

QES_KEYS = ["Q", "A", "W", "G", "D", "QC", "T", "TD"]
def qes(cs, names, w=None):
    X = np.array([[cs[n][k] for k in QES_KEYS] for n in names], float)
    rng = X.max(0) - X.min(0); rng[rng == 0] = 1
    N = 1 - (X - X.min(0)) / rng
    w = np.full(len(QES_KEYS), 1 / len(QES_KEYS)) if w is None else w
    return N @ w, N

def ranks_desc(v):
    return list((-np.asarray(v)).argsort().argsort() + 1)

if __name__ == "__main__":
    import os
    os.makedirs("results", exist_ok=True)
    out = {}
    C = circuits(); Cn = circuits(rot="ROT_NOUNDO"); Cs = circuits(rot="ROT_STATIC")
    out["circuits"] = C
    out["rotate_back"] = {n: dict(static=Cs[n], noundo=Cn[n], undo=C[n]) for n in ("SASI", "URAP")}
    # SASI block table
    b = 32; s = swaps(b)
    out["sasi_blocks"] = [("IDS xor K1 xor n1", 0, 96, 0, 3), ("IDS or K2", 160, 0, 32, 3),
                          ("+ n2 mod 2^32", 0, 128, 64, 161),
                          ("Rot(K1,n2) with rotate-back", 0, 4 * s + b, 2 * s, 6 * s + 1),
                          ("Rot(K2,n1) with rotate-back", 0, 4 * s + b, 2 * s, 6 * s + 1),
                          ("+ K1 mod 2^32", 0, 128, 64, 161)]
    out["ascon_round"] = ascon_round()
    A = {n: attack(C[n]) for n in ORDER}; As = {n: attack(C[n], "TD_serial") for n in ORDER}
    out["attack"] = A; out["attack_serial"] = As
    # QES and rankings
    q, N = qes(C, ORDER); rq = ranks_desc(q)
    cost = [A[n]["log2_Ttotal"] for n in ORDER]; rc = list(np.argsort(np.argsort(cost)) + 1)
    out["qes"] = dict(zip(ORDER, q.round(3).tolist())); out["qes_rank"] = dict(zip(ORDER, map(int, rq)))
    out["cost_rank"] = dict(zip(ORDER, map(int, rc)))
    out["tau"] = float(kendalltau(rq, rc)[0])
    ind = [QES_KEYS.index(k) for k in ("Q", "A", "G", "D", "T", "TD")]
    qi = N[:, ind].mean(1); out["qes_ind"] = dict(zip(ORDER, qi.round(3).tolist()))
    out["qes_ind_rank"] = dict(zip(ORDER, map(int, ranks_desc(qi))))
    # leave-one-out
    loo = {}
    for drop in ORDER:
        names = [n for n in ORDER if n != drop]; qq, _ = qes(C, names)
        loo[drop] = dict(zip(names, qq.round(3).tolist()))
    out["loo"] = loo
    # sensitivity: w = u / sum(u), u ~ U(0,1)^8
    rng = np.random.default_rng(20260916); n_w = 50000
    u = rng.random((n_w, len(QES_KEYS))); Wt = u / u.sum(1, keepdims=True)
    S = Wt @ N.T; R = (-S).argsort(1).argsort(1) + 1
    eq = np.array(rq)
    taus = np.array([kendalltau(eq, r)[0] for r in R])
    out["sens"] = dict(mean_rank=dict(zip(ORDER, R.mean(0).round(2).tolist())),
                       min_rank=dict(zip(ORDER, R.min(0).tolist())), max_rank=dict(zip(ORDER, R.max(0).tolist())),
                       holds=dict(zip(ORDER, ((R == eq).mean(0) * 100).round(1).tolist())),
                       tau_mean=float(taus.mean()), tau_min=float(taus.min()),
                       identical=float((R == eq).all(1).mean() * 100),
                       dist={n: [float((R[:, i] == k).mean() * 100) for k in range(1, 7)] for i, n in enumerate(ORDER)})
    np.save("results/sens_ranks.npy", R)
    # NIST, machines, surface code
    out["nist"] = {n: {md: nist(A[n]["log2_Ttotal"], A[n]["log2_TDtotal"], md)[1] for md in (40, 64, 96)} for n in ORDER}
    out["machines"] = {n: {md: machines(A[n]["log2_TDtotal"], md) for md in (40, 64, 96)} for n in ORDER}
    out["surface"] = {n: surface(A[n]["W_oracle"], A[n]["log2_TDtotal"]) for n in ORDER}
    minb = {}
    for n in PROTOCOLS:
        minb[n] = {}
        for md in (40, 64, 96):
            for bb in range(8, 257, 8):
                a = attack(ulw_circuit(n, bb))
                if nist(a["log2_Ttotal"], a["log2_TDtotal"], md)[1] >= 0:
                    minb[n][md] = bb; break
    out["min_b"] = minb
    sweep = {n: {bb: {md: nist(attack(ulw_circuit(n, bb))["log2_Ttotal"], attack(ulw_circuit(n, bb))["log2_TDtotal"], md)[1]
                      for md in (40, 64, 96)} for bb in range(16, 129, 8)} for n in PROTOCOLS}
    out["sweep"] = sweep
    # Pareto frontier under both depth assumptions
    def front(att):
        P = {n: (att[n]["W_oracle"], att[n]["depth_iter"]) for n in ORDER}
        return [n for n, (w, d) in P.items() if not any(w2 <= w and d2 <= d and (w2, d2) != (w, d) for w2, d2 in P.values())]
    out["frontier_msgpar"] = front(A); out["frontier_serial"] = front(As)
    # range of oracle T and iteration count
    out["oracle_T_bits"] = math.log2(max(A[n]["T_oracle"] for n in ORDER) / min(A[n]["T_oracle"] for n in ORDER))
    out["oracle_T_bits_ulw"] = math.log2(max(A[n]["T_oracle"] for n in PROTOCOLS) / min(A[n]["T_oracle"] for n in PROTOCOLS))
    out["iter_bits"] = max(A[n]["log2_Niter"] for n in ORDER) - min(A[n]["log2_Niter"] for n in ORDER)
    json.dump(out, open("results/numbers.json", "w"), indent=1, default=lambda o: o if not isinstance(o, np.generic) else o.item())
    # console summary
    print("QEV:"); [print(f"  {n:14s}", {k: C[n][k] for k in ['Q','A','W','G','X','CNOT','Toffoli','D','T','TD','QC']}) for n in ORDER]
    print("Attack:"); [print(f"  {n:14s} r={A[n]['r']} Toracle={A[n]['T_oracle']} per={A[n]['per_iter']} T=2^{A[n]['log2_Ttotal']:.1f} TD=2^{A[n]['log2_TDtotal']:.1f} (ser 2^{As[n]['log2_TDtotal']:.1f}) Wor={A[n]['W_oracle']} dit={A[n]['depth_iter']}") for n in ORDER]
    print("QES", out["qes"], "rank", out["qes_rank"]); print("QESind", out["qes_ind"], out["qes_ind_rank"])
    print("cost rank", out["cost_rank"], "tau", round(out["tau"], 3))
    print("LOO"); [print("  drop", d, v) for d, v in loo.items()]
    print("sens", {k: v for k, v in out["sens"].items() if k != "dist"})
    print("NIST", {n: {m: round(v, 1) for m, v in d.items()} for n, d in out["nist"].items()})
    print("machines", {n: {m: round(v, 1) for m, v in d.items()} for n, d in out["machines"].items()})
    print("surface", {n: (d, f"{q:.2e}", f"{y:.2e}") for n, (d, q, y) in out["surface"].items()})
    print("min_b", minb); print("frontier", out["frontier_msgpar"], out["frontier_serial"])
    print("oracle T bits all/ulw", round(out["oracle_T_bits"], 2), round(out["oracle_T_bits_ulw"], 2), "iter bits", round(out["iter_bits"], 1))
    print("rotate-back", {n: (Cs[n]['Toffoli'], Cn[n]['Toffoli'], C[n]['Toffoli'], Cs[n]['D'], Cn[n]['D'], C[n]['D']) for n in ('SASI','URAP')})
