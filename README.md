# A forced coordinate with nearly maximal erasure entropy

Yitzchak Shmalo

This repository contains a short paper on the Fourier entropy-influence (FEI) conjecture of Friedgut and Kalai,
together with the programs behind every numerical value in it.

Paper: [main.pdf](main.pdf) (source [main.tex](main.tex)).
Zenodo: [doi:10.5281/zenodo.24139723](https://doi.org/10.5281/zenodo.24139723).

## Summary

For a Boolean function f and a coordinate x, the erasure entropy E_x(f) is the conditional entropy of the event
"x is in U" given U minus x, where U is drawn from the spectral distribution. It never exceeds the binary entropy of
the influence, so E_x/Inf_x is at most log2(e/Inf_x).

F_k is the monotone function on 2k+1 variables obtained from MAJ_k(y) OR (x AND [y >= y0]), with the point
y0 + e_{m+1} removed from the majority (k = 2m+1), by replacing y_1, ..., y_m with the OR of two fresh variables and
y_{m+1}, ..., y_k with the AND of two. The coordinate x is the only variable lying both in a shortest minterm and in a
shortest maxterm. The paper proves

    E_x(F_k)/Inf_x(F_k) = (1 + o(1)) log2(1/Inf_x(F_k)) = (1 + o(1)) (k-1) log2(4/3),

so the ratio at x is asymptotic to its largest possible value. The relevant biased Fourier coefficients of the
majority are integrals over the unit circle; on the circle the integrand has an explicit modulus with one pair of
peaks, and its phase advances by arccos(1/3), an irrational multiple of pi, from one class of coefficients to the next.

## Code

All programs are Python 3; `forced_rule.py` and `check_gk.py` also need NumPy.

* `code/forced_exact.py` computes E_x/Inf_x exactly (rational arithmetic) from the double sum of Lemma 3. Output for
  odd k <= 141: `code/forced_exact_out.txt`.
* `code/forced_rule.py` (imports `code/check_gk.py`) builds F_k for k = 3, 5, 7, 9 from its full truth table, checks
  that x is the only variable in both a shortest minterm and a shortest maxterm, and compares the ratio from the
  Walsh-Hadamard transform with the biased formula. Output: `code/forced_rule_out.txt`.
* `code/laplace_check.py` compares the Laplace approximation of Lemma 6 with the exact coefficients over a window of
  classes. Output: `code/laplace_check_out.txt`.
* `code/circle_max.py` scans the modulus of the integrand on the unit circle for local maxima. Output:
  `code/circle_max_out.txt`.

## Related papers

1. [The one-coordinate step in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fourier-entropy-influence-one-coordinate)
2. [Good coordinates in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-good-coordinates)
3. [A variance form of Han's bound on the Fourier entropy](https://github.com/yspennstate/fei-variance-han-bound)
4. [Four equivalent forms of the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-equivalent-forms)
5. [A filtration bound for the erasure entropies of Boolean functions](https://github.com/yspennstate/fei-filtration-witness)
6. [Erasure entropies of sparse linear DNFs](https://github.com/yspennstate/fei-sparse-linear-dnf)
7. [The Fourier entropy of linear monotone DNFs on few variables](https://github.com/yspennstate/fei-linear-dnf-census)
8. [Hidden coordinates and selection rules in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-hidden-coordinates)
9. [The erasure entropy of a forced coordinate under every bias](https://github.com/yspennstate/fei-forced-coordinate-bias)

## License

MIT, see [LICENSE](LICENSE).
