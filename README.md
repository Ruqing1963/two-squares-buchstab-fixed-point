# Sums of two squares and Heegner norms with few prime factors, XIV

**The Buchstab fixed point of a vector sieve as an obstacle problem, causal computation, and variable dispersion levels**

Ruqing Chen (GUT Geoservice Inc., Montreal). Draft; continues Part XIII (doi:10.5281/zenodo.24125302).

- `paper/two_squares_part14.tex`: the paper.
  - Section 2: the fixed point of the recursion of Part XIII (Lemma 4.3 there) in the gauge u = a^{-1/2} Phi, v = a^{-1/2} G is an
    obstacle problem; causality in the level gives a single O(n^2) sweep (Proposition 2.1); half-level scaling (Proposition 2.2);
    the free boundary degenerates (optimal cut = smallest admissible cut on a < X - 1.5 e).
  - Section 3: level 1/2, one-sided certification on n up to 12800: r_0 = 7 / 9 / 15 (Theorem 1.1), D = 23 table 9/8, 11/10, 17/16.
  - Section 4: Hypothesis L (Lichtman-type level for {q1 q2 r}, residue N) and its reduction to Lichtman's Propositions 3.1/3.2;
    the sieve constant G_1 (recomputed 5.581 against Li's 6.069); r_0 = 6 / 8 / 14 (Theorem 1.2), P_6.
  - Section 5: the problem of Nath-Xie (p - 1 sum of two squares, p + 2 = P_r): their computation files integrate the linear sieve
    without the factor 1/s; with the correct functions their inequality gives r = 84, the fixed point r = 65 (one-sided) / 56.
- `code/`: `hjb_sieve_engine.py` (engine; modes verify, phase, app2), `causal_fixed_point.py` (modes check, run, boundary, lowmem),
  `lichtman_level.py`, `nath_xie_app.py`, `nath_xie_check.py`, `switch_probe.py`, with `semilinear_sieve.py` and
  `iterated_buchstab.py` from Part XIII.
- `results/`, `figures/`: outputs.

Run with `python <script>.py [mode]` (numpy, matplotlib).
