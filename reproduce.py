"""reproduce.py (version 2.0.0)

Reproduces every numerical result, Tables 1-2 and Figure 1 of

    M. Hawarey, "Exactly Three Single-Amplification Cycles in GCD-Augmented
    Fibonacci Recurrences", AIR Journal of Mathematics and Computational
    Sciences (2026), revised version (AIR-2026-000962-V2).

Output blocks are labelled with the section, table, figure or proposition of the
manuscript in which the values appear, in manuscript order. The script also writes
figure1.png (300 dpi) and numbers.json (all quoted values, keyed Q0..Q24).

Exact integer arithmetic throughout. Requires Python 3.8+ and matplotlib.
Usage:  python reproduce.py        (runtime about one minute)
In the author's project files the same script is kept as figures_and_numbers.py.
"""
import json
from math import gcd, sqrt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = {}


def put(key, value):
    OUT[key] = value


# ------------------------------------------------------------------ basics
def fib_lucas(n_max):
    F = {-1: 1, 0: 0, 1: 1}
    for n in range(2, n_max + 1):
        F[n] = F[n - 1] + F[n - 2]
    L = {q: F[q - 1] + F[q + 1] for q in range(n_max)}
    return F, L


def seq(a, b, n_max):
    """C(1..n_max) for C(n) = C(n-1) + C(n-2) + gcd(C(n-1), C(n-2))."""
    C = [None, a, b]
    while len(C) <= n_max:
        C.append(C[-1] + C[-2] + gcd(C[-1], C[-2]))
    return C


def R(al, be):
    """Cofactor map R(alpha, beta) = (beta/g, (alpha+beta+1)/g), step multiplier g = gcd(beta, alpha+1)."""
    g = gcd(be, al + 1)
    return be // g, (al + be + 1) // g, g


def orbit(state):
    """States of the cycle through `state` and their step multipliers."""
    orb, gs, t = [state], [], state
    while True:
        x, y, g = R(*t)
        gs.append(g)
        t = (x, y)
        if t == state:
            return orb, gs
        orb.append(t)


def prod(xs):
    p = 1
    for x in xs:
        p *= x
    return p


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def fmt_factor(f):
    return "*".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(f.items()))


def first_lock_in(C, q, G, n_max):
    """Least n with C(m+q) = G*C(m) for all n <= m < n_max."""
    return next(n for n in range(1, n_max) if all(C[m + q] == G * C[m] for m in range(n, n_max)))


F, L = fib_lucas(200)
PHI = (1 + sqrt(5)) / 2
THREE = {(1, 2): "fixed point (1,2)", (3, 5): "2-cycle (3,5)", (29, 47): "4-cycle (29,47)"}


def show(title):
    print(f"\n== {title}")


# ------------------------------------------------------------------ Section 4
show("Section 4, after Theorem 7: candidate states, q = 1..8")
cands = []
for q in range(1, 9):
    if q % 2 == 0:
        g, st = L[q], ((L[q] - 1) * F[q + 1] - 1, (L[q] - 1) * F[q + 2] - 1)
    else:
        g, st = L[q] + 1, (F[q + 1], F[q + 2])
    cands.append([q, g, list(st)])
    print(f"q = {q}: cycle multiplier g = {g}, candidate state = {st}")
put("Q0", cands)

# ------------------------------------------------------------------ Section 5
show("Section 5, after Lemma 8: non-coprime candidates")
lem8 = {}
for q in (6, 8):
    st = ((L[q] - 1) * F[q + 1] - 1, (L[q] - 1) * F[q + 2] - 1)
    lem8[f"q{q}"] = [list(st), gcd(*st)]
    print(f"q = {q}: state = {st}, gcd = {gcd(*st)}")
put("Q12", lem8)

show("Table 1: the three single-amplification cycles")
table1 = []
for q, g, st in [(1, 2, (1, 2)), (2, 3, (3, 5)), (4, 7, (29, 47))]:
    orb, gs = orbit(st)
    u, v = st[0] + 1, st[1] + 1
    D = g * g - L[q] * g + (-1) ** q
    N = v * v - u * v - u * u
    assert len(orb) == q and gs[:-1] == [1] * (q - 1) and gs[-1] == g and N == -(g - 1) ** 2
    table1.append({"q": q, "g": g, "cycle": [list(o) for o in orb], "uv": [u, v], "D": D, "N": N})
    print(f"q = {q}, cycle multiplier g = {g}: cycle {orb}, (u, v) = {(u, v)}, D = {D}, N(u, v) = {N}")
put("Q13", table1)

show("Section 5, boundary cases (A083658; closed forms of the author's earlier work), p <= 49")
def state_at(C, n):
    d = gcd(C[n - 1], C[n])
    return C[n - 1] // d, C[n] // d
C = seq(1, 1, 60)
bc11 = state_at(C, 4) == (3, 5) and all(C[n + 2] == 3 * C[n] for n in range(3, 55))
even_ok, odd_ok = True, True
for p in range(1, 50):
    C = seq(1, p, 60)
    if p % 2 == 0:
        even_ok &= state_at(C, 4) == (1, 2) and all(C[n + 1] == 2 * C[n] for n in range(3, 55))
    elif p % 3:
        odd_ok &= state_at(C, 6) == (3, 5) and all(C[n + 2] == 3 * C[n] for n in range(5, 55))
# published values: [2, Table 3] (p = 1, 2, 5, 7, 10; n <= 12) and [3, Table 2] (odd p <= 49 prime to 3;
# its columns labelled C(5), C(6) hold C(6), C(7))
T3 = {1: [1, 1, 3, 5, 9, 15, 27, 45, 81, 135, 243, 405], 2: [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048],
      5: [1, 5, 7, 13, 21, 35, 63, 105, 189, 315, 567, 945], 7: [1, 7, 9, 17, 27, 45, 81, 135, 243, 405, 729, 1215],
      10: [1, 10, 12, 24, 48, 96, 192, 384, 768, 1536, 3072, 6144]}
t3_ok = all(seq(1, p, 12)[1:] == row for p, row in T3.items())
T2b = {1: (15, 27), 5: (35, 63), 7: (45, 81), 11: (65, 117), 13: (75, 135), 17: (95, 171), 19: (105, 189),
       23: (125, 225), 25: (135, 243), 29: (155, 279), 37: (195, 351), 41: (215, 387), 43: (225, 405),
       47: (245, 441), 49: (255, 459)}
t2b_ok = all(tuple(seq(1, p, 7)[6:8]) == v for p, v in T2b.items())
print(f"seed (1,1): state (3,5) at index 4 and C(n+2) = 3C(n) for n >= 3: {bc11}")
print(f"even p <= 49: state (1,2) at index 4 and C(n+1) = 2C(n) for n >= 3: {even_ok}")
print(f"odd p <= 49 prime to 3: state (3,5) at index 6 and C(n+2) = 3C(n) for n >= 5: {odd_ok}")
print(f"agreement with [2, Table 3]: {t3_ok}; with [3, Table 2] (columns = C(6), C(7)): {t2b_ok}")
put("Q17", {"(1,1)": bc11, "even": even_ok, "odd_coprime3": odd_ok, "MRS1_T3": t3_ok, "MRS2_T2": t2b_ok})

show("Proposition 12: seeds p = 9, 15, 27")
prop12 = {}
for p in (9, 15, 27):
    C = seq(1, p, 420)
    n_p = first_lock_in(C, 4, 7, 400)
    d = gcd(C[n_p], C[n_p + 1])
    base = [C[n_p + r] for r in range(4)]
    assert base == [d * c for c in (29, 47, 77, 125)]
    prop12[p] = {"n_p": n_p, "d_p": d, "base": base}
    print(f"p = {p}: C(n+4) = 7C(n) for all n >= {n_p} (and not for n = {n_p - 1}); d_p = {d} = 3^{factor(d).get(3, 0)}; "
          f"C_p(n_p..n_p+3) = {base}")
put("Q8", {p: prop12[p]["n_p"] for p in prop12})
put("Q9", {p: prop12[p]["d_p"] for p in prop12})
put("Q10", {p: prop12[p]["base"] for p in prop12})
C = seq(1, 15, 7)
d = gcd(C[6], C[7])
print(f"worked example p = 15: C(1..7) = {C[1:]}, gcd(C(6), C(7)) = {d}, state at index 7 = {(C[6] // d, C[7] // d)}")
put("Q19", {"C15_1_7": C[1:8], "gcd_C6_C7": d, "sigma7": [C[6] // d, C[7] // d]})

show("Figure 1 caption: exact four-step ratio 7 before lock-in (6 <= n < n_p)")
early = {}
for p in (9, 15, 27):
    C = seq(1, p, 60)
    early[p] = [n for n in range(6, prop12[p]["n_p"]) if C[n + 4] == 7 * C[n]]
    print(f"p = {p}: C(n+4) = 7C(n) holds at n = {early[p] if early[p] else 'none'} before n_p = {prop12[p]['n_p']}")
put("Q23", early)


# ------------------------------------------------------------------ census engine
def limit(a, b, steps, keep=()):
    """Follow the states from seed (a, b) for at most `steps` applications of R.
    Returns ('cycle', least state, period, step multipliers, k) where k is the step at which a state
    first repeats, or ('none', {t: state after t steps for t in keep})."""
    d = gcd(a, b)
    s, seen, snap = (a // d, b // d), set(), {}
    for k in range(steps + 1):
        if k in keep:
            snap[k] = s
        if s in seen:
            orb, gs = orbit(s)
            return ("cycle", min(orb), len(orb), gs, k)
        seen.add(s)
        if k < steps:
            s = R(*s)[:2]
    return ("none", snap)


def label(res):
    if res[0] == "none":
        return "none"
    m = res[1]
    return THREE.get(m, f"period-{res[2]} cycle through {m}")


# one pass over the box 200 at 2,400 steps; the 1,200-step census and the box 120 are read off it
RUNS = {}
for a in range(1, 201):
    for b in range(1, 201):
        if gcd(a, b) == 1:
            RUNS[(a, b)] = limit(a, b, 2400, keep=(1200, 2400))


def census(n_box, steps):
    cnt, seeds = {}, {}
    for (a, b), res in RUNS.items():
        if a > n_box or b > n_box:
            continue
        k = "none" if (res[0] == "none" or res[4] > steps) else label(res)
        cnt[k] = cnt.get(k, 0) + 1
        seeds.setdefault(k, []).append((a, b))
    return cnt, seeds


ORDER = ["fixed point (1,2)", "2-cycle (3,5)", "4-cycle (29,47)", "period-22 cycle through (97, 157)",
         "period-65 cycle through (4003, 6477)", "none"]


def print_census(title, cnt):
    tot = sum(cnt.values())
    show(title)
    keys = ORDER + sorted(k for k in cnt if k not in ORDER)
    for k in keys:
        v = cnt.get(k, 0)
        print(f"{k:38s} {v:7,d}  {100 * v / tot:5.1f}%")
    single = sum(cnt.get(k, 0) for k in ORDER[:3])
    print(f"{'total':38s} {tot:7,d};  three single-amplification cycles: {single:,} ({100 * single / tot:.1f}%)")
    return tot, single


pct = lambda c, t, k: round(100 * c.get(k, 0) / t, 1)
c120, seeds120 = census(120, 1200)
t120, s120 = print_census("Table 2, first column: coprime seeds 1 <= a, b <= 120, at most 1,200 steps", c120)
put("Q1", t120)
put("Q2", [c120["fixed point (1,2)"], pct(c120, t120, "fixed point (1,2)")])
put("Q3", [c120["2-cycle (3,5)"], pct(c120, t120, "2-cycle (3,5)")])
put("Q4", [c120["4-cycle (29,47)"], pct(c120, t120, "4-cycle (29,47)")])
other120 = {k: v for k, v in c120.items() if k not in ORDER[:3] and k != "none"}
put("Q5", [sum(other120.values()), round(100 * sum(other120.values()) / t120, 1), other120])
put("Q6", [c120["none"], pct(c120, t120, "none")])
put("Q7", [s120, round(100 * s120 / t120, 1)])

c200, seeds200 = census(200, 1200)
t200, s200 = print_census("Table 2, second column: coprime seeds 1 <= a, b <= 200, at most 1,200 steps", c200)
put("Q20", {"total": t200, "counts": c200, "pct": {k: pct(c200, t200, k) for k in c200},
            "single": [s200, round(100 * s200 / t200, 1)]})

show("Section 5: robustness (step limit doubled to 2,400)")
c120b, _ = census(120, 2400)
c200b, _ = census(200, 2400)
print(f"box 120: counts unchanged at 2,400 steps: {c120b == c120}")
print(f"box 200: counts unchanged at 2,400 steps: {c200b == c200}")
put("Q18", {"N120_steps2400_equal_to_steps1200": c120b == c120, "N200_steps2400_equal_to_steps1200": c200b == c200,
            "N200_total": t200, "N200_29_47_pct": pct(c200, t200, "4-cycle (29,47)"),
            "N200_3_5_pct": pct(c200, t200, "2-cycle (3,5)")})

show("Section 5: growth of the seeds with no detected cycle (box 120)")
none120 = [RUNS[ab] for ab in seeds120["none"]]
def digits_stats(t):
    dg = sorted(len(str(res[1][t][1])) for res in none120)
    return [dg[0], dg[len(dg) // 2], dg[-1]]
g1200, g2400 = digits_stats(1200), digits_stats(2400)
print(f"{len(none120)} seeds; digits of beta after 1,200 steps (min, median, max): {g1200}")
print(f"{len(none120)} seeds; digits of beta after 2,400 steps (min, median, max): {g2400}")
d45 = [len(str(RUNS[(1, 45)][1][t][1])) for t in (1200, 2400)]
print(f"seed (1,45): no repeated state; digits of beta after 1,200 and 2,400 steps: {d45}")
put("Q22", {"n_none": len(none120), "digits_1200": g1200, "digits_2400": g2400, "seed_1_45": d45})

show("Section 5: other odd multiples of 3 up to 49")
others = {}
for p in (3, 21, 33, 39, 45):
    res = limit(1, p, 1200)
    lab = label(res)
    if lab in ("2-cycle (3,5)", "4-cycle (29,47)"):
        q, G = (2, 3) if lab == "2-cycle (3,5)" else (4, 7)
        n0 = first_lock_in(seq(1, p, 420), q, G, 400)
        others[p] = [lab, n0]
        print(f"p = {p}: reaches the {lab}; C(n+{q}) = {G}C(n) for n >= {n0}")
    else:
        others[p] = [lab, None]
        print(f"p = {p}: {'no repeated state within 1,200 steps' if lab == 'none' else lab}")
put("Q15", others)

show("Section 5: values for the record (cf. the correction notices of [2] and [3])")
C21, C27, C5 = seq(1, 21, 10), seq(1, 27, 10), seq(1, 5, 7)
print(f"C_21(9), C_21(10) = {C21[9]}, {C21[10]};  C_27(7..10) = {C27[7:11]}")
print(f"p = 5: C(6), C(7) = {C5[6]}, {C5[7]} = 5(p+2), 9(p+2)")
put("Q16", {"C21_n9_10": [C21[9], C21[10]], "C27_n7_10": C27[7:11], "MRS2_T2_shift_example_p5": [C5[6], C5[7]]})

show("Section 5: the multi-amplification cycles found in the census")
multi = {}
for name, start in (("period-22", (1, 87)), ("period-65", (3, 187))):
    res = RUNS[start]
    m, per, gs = res[1], res[2], res[3]
    G = prod(gs)
    amps = [g for g in gs if g > 1]
    key = label(res)
    multi[name] = {"least_state": list(m), "period": per, "amplifications": len(amps), "multipliers": amps,
                   "G": G, "factorization": fmt_factor(factor(G)), "n_box120": c120.get(key, 0),
                   "n_box200": c200.get(key, 0)}
    print(f"{name} cycle reached from seed {start}: least state {m}, period {per}, {len(amps)} amplifications "
          f"with step multipliers {amps}; cycle multiplier G = {G:,} = {fmt_factor(factor(G))}")
seeds65 = seeds200.get("period-65 cycle through (4003, 6477)", [])
print(f"seeds in the box 200 that reach the period-65 cycle: {seeds65}")
multi["period-65"]["seeds_box200"] = seeds65
put("Q11", {"min_state": multi["period-22"]["least_state"], "period": 22,
            "amplifications": multi["period-22"]["multipliers"], "multiplier": multi["period-22"]["G"],
            "factorization": multi["period-22"]["factorization"],
            "check": multi["period-22"]["G"] == 3 ** 6 * 5 * 11})
put("Q21", multi["period-65"])


# ------------------------------------------------------------------ Proposition 13
show("Proposition 13: cycle equation and the bound G > phi^q for every cycle found")
def matpow(q):
    A, Q = [[1, 0], [0, 1]], [[0, 1], [1, 1]]
    for _ in range(q):
        A = [[A[i][0] * Q[0][j] + A[i][1] * Q[1][j] for j in range(2)] for i in range(2)]
    return A


def prop13(state):
    orb, gs = orbit(state)
    q = len(orb)
    last = max(i for i in range(q) if gs[i] > 1)          # rotate to start after an amplification
    orb, gs = orb[last + 1:] + orb[:last + 1], gs[last + 1:] + gs[:last + 1]
    pos = [i for i in range(q) if gs[i] > 1]
    qs = [pos[0] + 1] + [pos[i] - pos[i - 1] for i in range(1, len(pos))]
    g = [gs[i] for i in pos]
    G = prod(g)
    u = (orb[0][0] + 1, orb[0][1] + 1)
    Qq = matpow(q)
    lhs = (G * u[0] - Qq[0][0] * u[0] - Qq[0][1] * u[1], G * u[1] - Qq[1][0] * u[0] - Qq[1][1] * u[1])
    rhs, pre = [0, 0], 1
    for i in range(len(g)):
        s = sum(qs[i + 1:])
        rhs[0] += (g[i] - 1) * pre * F[s + 1]
        rhs[1] += (g[i] - 1) * pre * F[s + 2]
        pre *= g[i]
    bound = G >= (L[q] if q % 2 == 0 else L[q] + 1)
    D = G * G - L[q] * G + (-1) ** q
    return {"q": q, "k": len(g), "segments": qs, "multipliers": g, "G": G, "D": D, "eq_holds": lhs == tuple(rhs),
            "G_ge_Lucas_bound": bound, "G_over_phi_q": round(G / PHI ** q, 4)}


p13 = {}
for st in [(1, 2), (3, 5), (29, 47), tuple(multi["period-22"]["least_state"]), tuple(multi["period-65"]["least_state"])]:
    r = prop13(st)
    p13[str(st)] = r
    print(f"cycle through {st}: q = {r['q']}, k = {r['k']}, cycle equation holds: {r['eq_holds']}, "
          f"G >= {'L_q' if r['q'] % 2 == 0 else 'L_q + 1'}: {r['G_ge_Lucas_bound']}, G/phi^q = {r['G_over_phi_q']}, "
          f"D = {r['D']:,}")
put("Q24", p13)

# ------------------------------------------------------------------ Figure 1
show("Figure 1")
NMIN, NMAX = 6, 40
style = {9: ("#2a78d6", "o"), 15: ("#eb6834", "s"), 27: ("#1baf7a", "^")}
plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix",
                     "font.size": 10, "axes.linewidth": 0.6})
fig, axes = plt.subplots(3, 1, figsize=(6.5, 4.6), dpi=300, sharex=True, sharey=True)
for ax, (p, (col, mk)) in zip(axes, style.items()):
    C = seq(1, p, NMAX + 5)
    ns = list(range(NMIN, NMAX + 1))
    rs = [C[n + 4] / C[n] for n in ns]
    n_p = prop12[p]["n_p"]
    ax.axhline(7, color="#9a9a9a", lw=0.6, zorder=1)
    ax.axvline(n_p, color="#9a9a9a", lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.plot(ns, rs, color=col, lw=1.2, zorder=2)
    pre = [(n, r) for n, r in zip(ns, rs) if n < n_p]
    post = [(n, r) for n, r in zip(ns, rs) if n >= n_p]
    if pre:  # hollow markers before lock-in
        ax.plot(*zip(*pre), ls="none", marker=mk, ms=4.2, mfc="white", mec=col, mew=1.0, zorder=3)
    ax.plot(*zip(*post), ls="none", marker=mk, ms=4.2, color=col, zorder=3)
    ax.text(0.99, 0.95, f"$p = {p}$:  $C_p(n+4) = 7\\,C_p(n)$ for all $n \\geq {n_p}$",
            transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color="#1a1a1a")
    ax.set_yticks([6.9, 7.0, 7.1])
    ax.grid(axis="y", color="#ececec", lw=0.5, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_ylim(6.84, 7.16)
axes[-1].set_xlim(NMIN - 0.8, NMAX + 0.8)
axes[-1].set_xlabel("$n$")
axes[1].set_ylabel("$C_p(n+4)\\,/\\,C_p(n)$")
fig.tight_layout(h_pad=0.6)
fig.savefig("figure1.png", dpi=300)
print("written: figure1.png")

with open("numbers.json", "w") as fh:
    json.dump({"numbers": OUT}, fh, indent=1, default=str)
print("written: numbers.json")
