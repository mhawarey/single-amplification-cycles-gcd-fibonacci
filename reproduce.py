"""reproduce.py

Reproduces the numerical results, Tables 1-2 and Figure 1 of

    M. Hawarey, "Exactly three single-amplification cycles in GCD-augmented
    Fibonacci recurrences", AIR Journal of Mathematics and Computational
    Sciences (2026).

Exact integer arithmetic throughout. Requires Python 3.8+ and matplotlib.
Usage:  python reproduce.py      (prints all values; writes figure1.png)
"""
from math import gcd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
    """Cofactor map R(alpha, beta) = (beta/g, (alpha+beta+1)/g), g = gcd(beta, alpha+1)."""
    g = gcd(be, al + 1)
    return be // g, (al + be + 1) // g, g


def orbit(state):
    """Cycle through `state` and its step multipliers."""
    orb, gs, t = [state], [], state
    while True:
        x, y, g = R(*t)
        gs.append(g)
        t = (x, y)
        if t == state:
            return orb, gs
        orb.append(t)


def limit_cycle(a, b, steps):
    """(least state, period, multipliers) of the cycle reached from seed (a, b), or None."""
    d = gcd(a, b)
    s, seen = (a // d, b // d), set()
    for _ in range(steps + 1):
        if s in seen:
            orb, gs = orbit(s)
            return min(orb), len(orb), gs
        seen.add(s)
        s = R(*s)[:2]
    return None


def first_lock_in(C, q, G, n_max):
    """Least n with C(m+q) = G*C(m) for all n <= m < n_max."""
    return next(n for n in range(1, n_max) if all(C[m + q] == G * C[m] for m in range(n, n_max)))


THREE = {(1, 2): "(1,2)", (3, 5): "(3,5)", (29, 47): "(29,47)"}
F, L = fib_lucas(120)


def show(title):
    print(f"\n== {title}")


# ------------------------------------------------------------------ Section 4
show("Section 4: candidate states of Theorem 7, q = 1..8")
for q in range(1, 9):
    if q % 2 == 0:
        g, st = L[q], ((L[q] - 1) * F[q + 1] - 1, (L[q] - 1) * F[q + 2] - 1)
    else:
        g, st = L[q] + 1, (F[q + 1], F[q + 2])
    print(f"q = {q}: g = {g}, state = {st}")

show("Section 5: Lemma 8 examples (non-coprime candidates)")
for q in (6, 8):
    st = ((L[q] - 1) * F[q + 1] - 1, (L[q] - 1) * F[q + 2] - 1)
    print(f"q = {q}: state = {st}, gcd = {gcd(*st)}")

# ------------------------------------------------------------------ Table 1
show("Table 1: the three single-amplification cycles")
for q, g, st in [(1, 2, (1, 2)), (2, 3, (3, 5)), (4, 7, (29, 47))]:
    orb, gs = orbit(st)
    u, v = st[0] + 1, st[1] + 1
    D = g * g - L[q] * g + (-1) ** q
    N = v * v - u * v - u * u
    assert len(orb) == q and gs[:-1] == [1] * (q - 1) and gs[-1] == g and N == -(g - 1) ** 2
    print(f"q = {q}, g = {g}: cycle {orb}, (u, v) = {(u, v)}, D = {D}, N = {N}")

# ------------------------------------------------------------------ boundary cases
show("Section 5: boundary cases (A083658; earlier closed forms), p <= 49")
def state_at(C, n):
    d = gcd(C[n - 1], C[n])
    return C[n - 1] // d, C[n] // d
C = seq(1, 1, 60)
ok = state_at(C, 4) == (3, 5) and all(C[n + 2] == 3 * C[n] for n in range(3, 55))
for p in range(1, 50):
    C = seq(1, p, 60)
    if p % 2 == 0:
        ok &= state_at(C, 4) == (1, 2) and all(C[n + 1] == 2 * C[n] for n in range(3, 55))
    elif p % 3:
        ok &= state_at(C, 6) == (3, 5) and all(C[n + 2] == 3 * C[n] for n in range(5, 55))
print("seed (1,1) -> (3,5); even p -> (1,2); odd p prime to 3 -> (3,5):", "confirmed" if ok else "FAILED")

# ------------------------------------------------------------------ Proposition 12
show("Proposition 12: seeds p = 9, 15, 27")
n_lock = {}
for p in (9, 15, 27):
    C = seq(1, p, 420)
    n_p = first_lock_in(C, 4, 7, 400)
    d = gcd(C[n_p], C[n_p + 1])
    base = [C[n_p + r] for r in range(4)]
    assert base == [d * c for c in (29, 47, 77, 125)]
    n_lock[p] = n_p
    print(f"p = {p}: C(n+4) = 7C(n) for all n >= {n_p}; d_p = {d}; C_p(n_p..n_p+3) = {base}")
C = seq(1, 15, 7)
d = gcd(C[6], C[7])
print(f"worked example p = 15: C(1..7) = {C[1:]}, gcd(C(6), C(7)) = {d}, state = {(C[6] // d, C[7] // d)}")

show("Section 5: other odd multiples of 3 up to 49")
for p in (3, 21, 33, 39, 45):
    res = limit_cycle(1, p, 1200)
    lab = THREE.get(res[0], "other") if res else "no repeated state within 1,200 steps"
    if lab in ("(3,5)", "(29,47)"):
        q, G = (2, 3) if lab == "(3,5)" else (4, 7)
        print(f"p = {p}: enters {lab}; C(n+{q}) = {G}C(n) for n >= {first_lock_in(seq(1, p, 420), q, G, 400)}")
    else:
        print(f"p = {p}: {lab}")

show("Section 5: values of C_21 and C_27 (for the record)")
C21, C27 = seq(1, 21, 10), seq(1, 27, 10)
print(f"C_21(9), C_21(10) = {C21[9]}, {C21[10]};  C_27(7..10) = {C27[7:11]}")
C5 = seq(1, 5, 7)
print(f"p = 5: C(6), C(7) = {C5[6]}, {C5[7]} = 5(p+2), 9(p+2)")

# ------------------------------------------------------------------ Table 2
def census(n_box, steps):
    c = {"(1,2)": 0, "(3,5)": 0, "(29,47)": 0, "other cycle": 0, "no repeated state": 0}
    for a in range(1, n_box + 1):
        for b in range(1, n_box + 1):
            if gcd(a, b) == 1:
                res = limit_cycle(a, b, steps)
                c[THREE.get(res[0], "other cycle") if res else "no repeated state"] += 1
    return c

show("Table 2: census of coprime seeds 1 <= a, b <= 120 (1,200 steps)")
c120 = census(120, 1200)
tot = sum(c120.values())
for k, v in c120.items():
    print(f"{k:18s} {v:6d}  {100 * v / tot:5.1f}%")
single = c120["(1,2)"] + c120["(3,5)"] + c120["(29,47)"]
print(f"{'total':18s} {tot:6d};  three cycles together: {single} ({100 * single / tot:.1f}%)")
print("counts unchanged with 2,400 steps:", census(120, 2400) == c120)
c200 = census(200, 1200)
t200 = sum(c200.values())
print(f"box 1 <= a, b <= 200: {t200} seeds; (29,47) {100 * c200['(29,47)'] / t200:.1f}%, "
      f"(3,5) {100 * c200['(3,5)'] / t200:.1f}%")

# ------------------------------------------------------------------ Remark 6
show("Remark 6: the cycle reached from seed (1, 87)")
m, per, gs = limit_cycle(1, 87, 2000)
mult = 1
for g in gs:
    mult *= g
print(f"least state {m}, period {per}, amplifications {[g for g in gs if g > 1]}, "
      f"multiplier {mult} = 3^6*5*11: {mult == 3**6 * 5 * 11}")

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
    n_p = n_lock[p]
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
