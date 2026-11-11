"""On |z| = 1, where is |G(z) z^{-m}|^{1/m} largest?  G is the biased MAJ generating factor of class (i, j) = (a m, b m):
    g(phi) = |z-1|^(a+b) |(1-P)+Pz|^(1-a) |P+(1-P)z|^(1-b),  z = e^{i phi},  P = 3/4.
For each (a, b) in a window around (1/4, 1/4) prints the maximising phi, the max, the value at phi = pi (the saddle
w = -1) and the second-largest local max, to check that the conjugate saddle pair is the only peak.
Usage: python circle_max.py [halfwidth]"""
import sys
from math import cos, sin, pi, log, exp, acos

P = 0.75
N = 200000


def g(phi, a, b):
    x, y = cos(phi), sin(phi)
    t1 = (x - 1) ** 2 + y ** 2
    t2 = ((1 - P) + P * x) ** 2 + (P * y) ** 2
    t3 = (P + (1 - P) * x) ** 2 + ((1 - P) * y) ** 2
    return exp(0.5 * ((a + b) * log(t1) + (1 - a) * log(t2) + (1 - b) * log(t3))) if t1 > 0 else 0.0


hw = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
print(f"P = {P}; peaks of g on 0 < phi < pi (g is even in phi)")
print("   a      b    phi_max/pi   g_max    g(pi)   other local maxima")
for a in (0.25 - hw, 0.25, 0.25 + hw):
    for b in (0.25 - hw, 0.25, 0.25 + hw):
        vals = [g(pi * (n + 0.5) / N, a, b) for n in range(N)]
        best = max(range(N), key=lambda n: vals[n])
        loc = [n for n in range(1, N - 1) if vals[n] > vals[n - 1] and vals[n] >= vals[n + 1] and n != best]
        other = ", ".join(f"{vals[n]:.4f}@{(n + 0.5) / N:.3f}pi" for n in loc) or "none"
        print(f"{a:.3f}  {b:.3f}  {(best + 0.5) / N:.6f}  {vals[best]:.6f}  {g(pi, a, b):.6f}  {other}")
print(f"arccos(1/3)/pi = {acos(1 / 3) / pi:.6f}")
