from math import isqrt
from polygon import G_pair

def Z_fun(N, counter=None):
    if N < 1: return (0, 0, 0)
    B = isqrt(N)
    tot = [0, 0, 0]
    for al in range(1, B + 1):
        for be in range(al, B + 1):
            w = 1 if al == be else 2
            g1 = G_pair(al, be, N)
            g2 = G_pair(al, be, N, B)
            if counter is not None: counter[0] += 1
            for k in range(3):
                tot[k] += w * (2 * g1[k] - g2[k])
    return tuple(tot)

def Z_brute(N):
    r = [0, 0, 0]
    for al in range(1, N + 1):
        for be in range(1, N + 1):
            for p in range(1, N + 1):
                for q in range(1, N + 1):
                    X = al * p + be * q; Y = max(al * q, be * p)
                    if X <= N and Y <= N:
                        r[0] += 1; r[1] += X + Y; r[2] += X * Y
    return tuple(r)

def mobius_prefix(M):
    mu = [1] * (M + 1); is_c = bytearray(M + 1); mu[0] = 0
    for i in range(2, M + 1):
        if not is_c[i]:
            for j in range(i, M + 1, i):
                is_c[j] = 1 if j > i else is_c[j]; mu[j] = -mu[j]
            for j in range(i * i, M + 1, i * i): mu[j] = 0
    P0 = [0] * (M + 1); P1 = [0] * (M + 1); P2 = [0] * (M + 1)
    for d in range(1, M + 1):
        P0[d] = P0[d - 1] + mu[d]; P1[d] = P1[d - 1] + d * mu[d]; P2[d] = P2[d - 1] + d * d * mu[d]
    return P0, P1, P2

def R_nearlinear(n, counter=None):
    if n < 2: return 0
    S = n * (n - 1) // 2
    P0, P1, P2 = mobius_prefix(n)
    W = 0
    d = 1
    while d <= n - 1:
        N = (n - 1) // d
        dh = (n - 1) // N
        m0 = P0[dh] - P0[d - 1]; m1 = P1[dh] - P1[d - 1]; m2 = P2[dh] - P2[d - 1]
        if m0 or m1 or m2:
            Z0, Z1, Z2 = Z_fun(N, counter)
            W += n * n * m0 * Z0 - n * m1 * Z1 + m2 * Z2
        d = dh + 1
    return 4 * S * S + 4 * W
