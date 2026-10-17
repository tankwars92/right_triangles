from fractions import Fraction as Fr
from math import comb, lcm

STATES = [(0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (0, 3)]

def P(r, n):
    if n <= 0: return 0
    if r == 0: return n
    if r == 1: return n * (n - 1) // 2
    if r == 2: return (n - 1) * n * (2 * n - 1) // 6
    return (n * (n - 1) // 2) ** 2

def _pmul(A, B):
    R = {}
    for (a1, b1), c1 in A.items():
        for (a2, b2), c2 in B.items():
            k = (a1 + a2, b1 + b2); R[k] = R.get(k, 0) + c1 * c2
    return R

def _w(j):
    return {(k, 0): Fr(comb(j, k)) for k in range(j)}

def _Pg(i):
    if i == 0: return {(0, 1): Fr(1), (0, 0): Fr(1)}
    if i == 1: return {(0, 2): Fr(1, 2), (0, 1): Fr(1, 2)}
    return {(0, 3): Fr(1, 3), (0, 2): Fr(1, 2), (0, 1): Fr(1, 6)}

_RECIP = {}
for (i, j) in STATES:
    poly = _pmul(_w(j), _Pg(i))
    den = 1
    for c in poly.values(): den = lcm(den, c.denominator)
    _RECIP[(i, j)] = (den, [(int(c * den), k) for k, c in poly.items() if c != 0])

def H(n, m, a, b):
    if n <= 0: return {s: 0 for s in STATES}
    if a == 0:
        c = b // m
        return {(i, j): c ** j * P(i, n) for (i, j) in STATES}
    if a >= m or b >= m:
        A, B = a // m, b // m
        bb = H(n, m, a % m, b % m)
        Pn = [P(r, n) for r in range(4)]
        out = {}
        for (i, j) in STATES:
            tot = 0
            for k in range(j + 1):
                coef = comb(j, k)
                for s in range(j - k + 1):
                    c2 = coef * comb(j - k, s) * A ** s * B ** (j - k - s)
                    if k == 0: tot += c2 * Pn[i + s]
                    else: tot += c2 * bb[(i + s, k)]
            out[(i, j)] = tot
        return out
    Y = (a * (n - 1) + b) // m
    if Y == 0: return {s: 0 for s in STATES}
    G = H(Y, a, m, m - b - 1)
    Q = [P(r, Y) for r in range(4)]
    out = {}
    for (i, j) in STATES:
        den, terms = _RECIP[(i, j)]
        s = 0
        for c, (ta, tb) in terms:
            s += c * (G[(ta, tb)] if tb >= 1 else Q[ta])
        out[(i, j)] = P(i, n) * Y ** j - s // den if s % den == 0 else None
    return out
