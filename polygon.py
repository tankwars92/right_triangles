from math import comb
from kernel import H, P, STATES

NEEDED = [(i, j) for j in range(0, 4) for i in range(0, 4) if i + j <= 3]

def fsums(L, R, a, b, m):
    if R < L: return {s: 0 for s in NEEDED}
    cnt = R - L + 1
    if a >= 0: base, sgn, a1, b1 = L, 1, a, a * L + b
    else:      base, sgn, a1, b1 = R, -1, -a, a * R + b
    c = b1 // m; r = b1 - c * m
    Hd = H(cnt, m, a1, r)
    Pc = [P(k, cnt) for k in range(4)]
    def T(k, kk): return Pc[k] if kk == 0 else Hd[(k, kk)]
    out = {}
    for (i, j) in NEEDED:
        tot = 0
        for k in range(i + 1):               
            cp = comb(i, k) * base ** (i - k) * (sgn ** k)
            for kk in range(j + 1):          
                tot += cp * comb(j, kk) * c ** (j - kk) * T(k, kk)
        out[(i, j)] = tot
    return out

def last_true(lo, hi, pred):
    if lo > hi or not pred(lo): return lo - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if pred(mid): lo = mid
        else: hi = mid - 1
    return lo

def G_pair(al, be, N, cap=None):
    qC = N // al
    if cap is not None: qC = min(qC, cap)
    pcap = N // be
    if cap is not None: pcap = min(pcap, cap)
    fdec = lambda p: (N - al * p) // be
    qcap = lambda p: min(qC, fdec(p))
    q0 = lambda p: (be * p + al - 1) // al
    if pcap < 1 or qC < 1: return (0, 0, 0)
    pmax = last_true(1, pcap, lambda p: qcap(p) >= 1)
    if pmax < 1: return (0, 0, 0)
    pT = last_true(1, pmax, lambda p: fdec(p) >= qC)
    pA = last_true(1, pmax, lambda p: q0(p) <= qcap(p))
    pA1 = last_true(1, pmax, lambda p: q0(p) - 1 <= qcap(p))
    cuts = sorted({0, pT, pA, pA1, pmax} & set(range(0, pmax + 1)))
    CONST = (0, qC, 1); DEC = (-al, N, be); Q0M = (be, -1, al)
    acc = [0, 0, 0]                           # шкала 6
    for idx in range(1, len(cuts)):
        s1, s2 = cuts[idx - 1] + 1, cuts[idx]
        if s1 > s2: continue
        qc_form = CONST if s1 <= pT else DEC
        ub_form = Q0M if s1 <= pA1 else qc_form
        D = fsums(s1, s2, *ub_form)
        g = lambda i, j: D[(i, j)]
        acc[0] += 6 * g(0, 1)
        acc[1] += 6 * (al + be) * g(1, 1) + 3 * be * g(0, 2) + 3 * be * g(0, 1)
        acc[2] += 6 * al * be * g(2, 1) + 3 * be * be * g(1, 2) + 3 * be * be * g(1, 1)
        if s1 <= pA:
            for form, sign in ((qc_form, 1), (Q0M, -1)):
                D = fsums(s1, s2, *form)
                g = lambda i, j: D[(i, j)]
                acc[0] += sign * 6 * g(0, 1)
                acc[1] += sign * (6 * al * g(1, 1) + 3 * (al + be) * g(0, 2) + 3 * (al + be) * g(0, 1))
                acc[2] += sign * (3 * al * al * g(1, 2) + 3 * al * al * g(1, 1)
                                  + al * be * (2 * g(0, 3) + 3 * g(0, 2) + g(0, 1)))
    assert all(x % 6 == 0 for x in acc), acc
    return tuple(x // 6 for x in acc)

def G_brute(al, be, N, cap=None):
    r = [0, 0, 0]
    pm = N // be if cap is None else min(N // be, cap)
    for p in range(1, pm + 1):
        qm = N // al if cap is None else min(N // al, cap)
        for q in range(1, qm + 1):
            X = al * p + be * q; Y = max(al * q, be * p)
            if X <= N and Y <= N:
                r[0] += 1; r[1] += X + Y; r[2] += X * Y
    return tuple(r)
