#!/usr/bin/env python3
"""Margin of Theorem 1.3 (omega(p+2) <= 15): the corrected Nath-Xie inequality (switched term 0.579909 C(theta1) W(theta2))
with the fixed point on finer grids.  One-sided rule = certified; trapezoid for comparison.
Usage: python nath_xie_margin.py [n ...]   (default 3200).   Output: results/nath_xie_margin.txt"""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
NS = [int(a) for a in sys.argv[1:]] or [3200]
sys.argv = ["x", "lib"]
sys.path.insert(0, HERE)
import hjb_sieve_engine as E          # noqa: E402
out = open(os.path.join(HERE, "..", "results", "nath_xie_margin.txt"), "a", encoding="utf8")
T0 = time.time()


def say(s=""):
    print(s, flush=True)
    out.write(s + "\n")
    out.flush()


g = np.euler_gamma
A, SW = 0.216213, 0.579909
fa, Fa = E.sieve_tables(1.0)


def Cth(th1):
    bb = np.linspace(0, np.sqrt(1 - 2 * th1), 4001)[1:]
    b = bb ** 2
    return np.trapz(2 * np.log((1 - b - th1) / th1) / (1 - b), bb)


def W(th2):
    s2 = np.arange(1.0, 0.49 / th2, 0.005)
    return float(((np.exp(-g) * E.lk(Fa, s2) / 2) / (th2 * (0.5 - s2 * th2))).min())


def best_r(Phi, G, n, th1, th2, w):
    h = 0.5 / n
    j = int(round(th1 / h))
    Xg = np.arange(1, n + 1) * h
    main = A * Phi[n, j] / (np.sqrt(th1) * th2)
    T = SW * Cth(th1) * w
    if main <= T:
        return None
    best = None
    for th in np.linspace(0.15, 0.25, 201):
        al = np.geomspace(th2, th, 800)
        I = np.trapz((1 - al / th) * np.interp(0.5 - al, Xg, G[1:, j]) / al, al)
        lam = (main - T) / (A * I / (np.sqrt(th1) * th2))
        V = 1 / lam + 1 / th
        if best is None or V < best[0]:
            best = (V, th, lam, main, T)
    return best


for n in NS:
    for th2 in (0.007, 0.008, 0.009):
        eng = E.Engine(th2, n, kappa_aux=1.0, nS=200)
        w = W(th2)
        for rule in ("one-sided", "trapezoid"):
            Phi, G = eng.solve(rule)
            best = None
            h = 0.5 / n
            for jj in range(int(round(0.390 / h)), int(round(0.408 / h)) + 1):
                t1 = jj * h
                b = best_r(Phi, G, n, t1, th2, w)
                if b and (best is None or b[0] < best[0]):
                    best = b + (t1,)
            say("n=%d theta2 %.4f %-9s: min 1/lam+1/th = %.4f => omega(p+2) <= %d (theta1 %.5f theta %.4f lam %.5f main %.4f sw %.4f) [%4.0f s]"
                % (n, th2, rule, best[0], int(np.ceil(best[0])) - 1, best[5], best[1], best[2], best[3], best[4], time.time() - T0))
            del Phi, G
