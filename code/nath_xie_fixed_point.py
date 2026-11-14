#!/usr/bin/env python3
"""Section 5 (same computation as code/part14_fixed_point.py of the repository nath-xie-sieve-audit): the Nath-Xie inequality with the Buchstab fixed point on
the inert primes of p - 1 (Part XIV, Section 5 / Part XIII, Lemma 4.3) in place of the rectangle bounds f, F, with
the switched term of Proposition 4.2 in the three versions R0, R1, R2 of verify_nath_xie.py.

Uses hjb_sieve_engine.py (and semilinear_sieve.py) in this folder.  Supersedes nath_xie_app.py (switched term as printed
in Nath-Xie, theta2 >= 0.002 only) and the H-table of nath_xie_check.py.
Output: results/nath_xie_fixed_point.txt.   Run time: about 15 minutes."""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
P14 = HERE
sys.path.insert(0, P14)
_argv = sys.argv[:]
sys.argv = ["x", "lib"]
import hjb_sieve_engine as E          # noqa: E402
sys.argv = _argv
out = open(os.path.join(HERE, "..", "results", "nath_xie_fixed_point.txt"), "w", encoding="utf8")
T0 = time.time()


def say(s=""):
    print(s, flush=True)
    out.write(s + "\n")
    out.flush()


def clock():
    return "[%6.0f s]" % (time.time() - T0)


g = np.euler_gamma
# constants from data/verify_nath_xie.txt (primes to 1e7)
A = 0.216213                      # C kappa e^{-3 gamma/2} / 4
SW = {"R0": 1.111744, "R1": 0.579909, "R2": 0.579909}
fa, Fa = E.sieve_tables(1.0)      # linear F2 (correct)


def Cth(th1):
    bb = np.linspace(0, np.sqrt(1 - 2 * th1), 4001)[1:]
    b = bb ** 2
    return np.trapz(2 * np.log((1 - b - th1) / th1) / (1 - b), bb)


def W(th1, th2, var):
    if var in ("R0", "R1"):
        return 1.0 / (th2 * min(th1, 0.5 - th2))
    s2 = np.arange(1.0, 0.49 / th2, 0.005)
    return float(((np.exp(-g) * E.lk(Fa, s2) / 2) / (th2 * (0.5 - s2 * th2))).min())


def best_r(Phi, G, n, th1, th2, var):
    h = 0.5 / n
    j = int(round(th1 / h))
    Xg = np.arange(1, n + 1) * h
    main = A * Phi[n, j] / (np.sqrt(th1) * th2)
    T = SW[var] * Cth(th1) * W(th1, th2, var)
    if main <= T:
        return None
    best = None
    for th in np.linspace(th2 + 0.005, min(0.45, 0.5 - th2 - 0.005), 240):
        al = np.geomspace(th2, th, 400)
        Gt = np.interp(0.5 - al, Xg, G[1:, j])
        I = np.trapz((1 - al / th) * Gt / al, al)
        lam = (main - T) / (A * I / (np.sqrt(th1) * th2))
        V = 1 / lam + 1 / th
        if best is None or V < best[0]:
            best = (V, th, lam, main, T, Phi[n, j])
    return best


GRID = {"R0": [0.0005, 0.0006, 0.0008, 0.001, 0.0015],
        "R1": [0.0008, 0.001, 0.0012, 0.0015, 0.002, 0.003],
        "R2": [0.004, 0.005, 0.006, 0.007, 0.008, 0.010, 0.012]}
say("Nath-Xie inequality with the Part XIV fixed point; min of 1/lambda + 1/theta subject to H >= 0  %s" % clock())
n = 400
BEST = {}
for var in ("R2", "R1", "R0"):
    for rule in ("rect", "trapezoid", "one-sided"):
        glob = None
        for th2 in GRID[var]:
            eng = E.Engine(th2, n, kappa_aux=1.0, nS=200)
            if rule == "rect":
                Phi = np.zeros((n + 1, n + 1))
                G = np.zeros((n + 1, n + 1))
                for i in range(1, n + 1):
                    Phi[i, 1:n], G[i, 1:n] = eng.rect_row(i)
            else:
                Phi, G = eng.solve(rule)
            for jj in range(int(0.30 / (0.5 / n)), n):
                t1 = jj * 0.5 / n
                if t1 + 2 * th2 > 0.5:
                    break
                b = best_r(Phi, G, n, t1, th2, var)
                if b and (glob is None or b[0] < glob[0]):
                    glob = b + (t1, th2)
        if glob:
            BEST[(var, rule)] = glob
            say("  %s %-9s n=%d: min 1/lam+1/th = %.3f => omega(p+2) <= %d  (th1 %.5f th2 %.4f th %.4f lam %.4f main %.3f sw %.3f Phi %.5f) %s"
                % (var, rule, n, glob[0], int(np.ceil(glob[0])) - 1, glob[6], glob[7], glob[1], glob[2], glob[3], glob[4], glob[5], clock()))
        else:
            say("  %s %-9s n=%d: no admissible point %s" % (var, rule, n, clock()))

say("\nrefinement at n = 1600 (one-sided is the certified rule): each rule at the theta2 of its own n = 400 optimum"
    "\nand at that of the other rule; theta1 rescanned on the fine grid within 0.005 of the n = 400 optimum")
nn = 1600
for var in ("R2", "R1", "R0"):
    pts = sorted({(BEST[(var, r)][6], BEST[(var, r)][7]) for r in ("trapezoid", "one-sided") if (var, r) in BEST})
    for (t1c, th2) in pts:
        eng = E.Engine(th2, nn, kappa_aux=1.0, nS=200)
        for rule in ("trapezoid", "one-sided"):
            Phi, G = eng.solve(rule)
            best = None
            for jj in range(int(round((t1c - 0.005) / (0.5 / nn))), int(round((t1c + 0.005) / (0.5 / nn))) + 1):
                t1 = jj * 0.5 / nn
                if t1 + 2 * th2 > 0.5:
                    break
                b = best_r(Phi, G, nn, t1, th2, var)
                if b and (best is None or b[0] < best[0]):
                    best = b + (t1,)
            if best:
                say("  %s n=%d %-9s theta2 %.4f: min 1/lam+1/th = %.3f => omega(p+2) <= %d (th1 %.5f th %.4f lam %.4f main %.3f sw %.3f) %s"
                    % (var, nn, rule, th2, best[0], int(np.ceil(best[0])) - 1, best[6], best[1], best[2], best[3], best[4], clock()))
say("done %s" % clock())