"""Two-peak Laplace formula for the biased MAJ tail q_m(i, j) against its exact value, uniformly over a window of classes.
H(z) = (z-1)^(i+j) X(z)^(m-i) Y(z)^(m-j+1), X = (1-P)+Pz, Y = P+(1-P)z, and the tail sum_{c>=m+1} [z^c] H = q_m is
the z^m coefficient of the polynomial Q = H/(z-1), so q_m = (1/2pi) int_{-pi}^{pi} Q(e^{i phi}) e^{-i m phi} d phi.
With Lam(phi) = log Q(e^{i phi}) - i m phi, the stationary point phi* near arccos(1/3) gives
    q_m ~ 2 Re[ exp(Lam(phi*)) / sqrt(2 pi (-Lam''(phi*))) ]
(the conjugate peak gives the conjugate). Prints, per k, the largest error relative to the envelope
|exp(Lam)| / sqrt(2 pi |Lam''|), over |i - m/4|, |j - m/4| <= sqrt(m), for this complex-saddle form and for the
real-point form T of the paper (Lemma 6; i, j there are u, v), expanded at the real maximiser phi0 of |Q|.
Usage: python laplace_check.py k [k ...]"""
import cmath
import sys
from fractions import Fraction
from math import isqrt, pi, sqrt
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from forced_exact import pmul, ppow

P = 0.75
Pf = Fraction(3, 4)


def lam_derivs(phi, i, j, m):
    z = cmath.exp(1j * phi)
    X, Y = (1 - P) + P * z, P + (1 - P) * z
    lam = (i + j - 1) * cmath.log(z - 1) + (m - i) * cmath.log(X) + (m - j + 1) * cmath.log(Y) - 1j * m * phi
    # d/dphi = i z d/dz
    s1 = (i + j - 1) / (z - 1) + (m - i) * P / X + (m - j + 1) * (1 - P) / Y          # d/dz of log Q
    s2 = -(i + j - 1) / (z - 1) ** 2 - (m - i) * P ** 2 / X ** 2 - (m - j + 1) * (1 - P) ** 2 / Y ** 2
    d1 = 1j * z * s1 - 1j * m
    d2 = -z * s1 - z * z * s2                                                          # d/dphi (i z s1)
    return lam, d1, d2


def laplace(i, j, m):
    phi = complex(1.2309594, 0.0)
    for _ in range(60):
        _, d1, d2 = lam_derivs(phi, i, j, m)
        phi -= d1 / d2
    lam, d1, d2 = lam_derivs(phi, i, j, m)
    term = cmath.exp(lam) / cmath.sqrt(2 * pi * (-d2))
    return 2 * term.real, abs(term), phi


def laplace_real(i, j, m):
    """The form used in the proof: expand at the real maximiser phi0 of |Q(e^{i phi})|, x0 = cos phi0 =
    (3 n2 - 5 n1)/(6m), keeping the linear phase beta = Im Lam'(phi0) inside the Gaussian."""
    n1, n2 = i + j - 1, 2 * m + 1 - i - j
    phi0 = cmath.acos((3 * n2 - 5 * n1) / (6 * m)).real
    lam, d1, d2 = lam_derivs(complex(phi0, 0.0), i, j, m)
    A = -d2
    term = cmath.exp(lam) / cmath.sqrt(2 * pi * A) * cmath.exp(d1 * d1 / (2 * A))   # d1 = i beta, so d1^2 = -beta^2
    return 2 * term.real, abs(term), abs(d1.real)


def exact(i, j, m):
    u, w = Pf.numerator, Pf.denominator
    poly = pmul(pmul(ppow([w - u, u], m - i), ppow([u, w - u], m - j + 1)), ppow([-1, 1], i + j))
    return float(Fraction(sum(poly[m + 1:]), w ** (2 * m - i - j + 1)))


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [41, 81, 141]:
        m = (k - 1) // 2
        c, r = m // 4, isqrt(m)
        worst, worst_r, re_d1, cnt, phis = 0.0, 0.0, 0.0, 0, []
        for i in range(max(1, c - r), c + r + 1):
            for j in range(max(1, c - r), c + r + 1):
                approx, env, phi = laplace(i, j, m)
                ex = exact(i, j, m)
                worst = max(worst, abs(approx - ex) / env)
                ar, envr, rd = laplace_real(i, j, m)
                worst_r = max(worst_r, abs(ar - ex) / envr)
                re_d1 = max(re_d1, rd)
                cnt += 1
                phis.append(phi)
        print(f"k={k} real-phi0 form: max |Laplace - exact| / |T| = {worst_r:.5f}  (max |Re Lam'(phi0)| = {re_d1:.1e})")
        print(f"k={k} m={m} classes={cnt} max |Laplace - exact| / envelope = {worst:.5f}  (times m: {worst * m:.3f})"
              f"  phi* real in [{min(p.real for p in phis):.4f}, {max(p.real for p in phis):.4f}],"
              f" |imag| <= {max(abs(p.imag) for p in phis):.4f}")
        sys.stdout.flush()
