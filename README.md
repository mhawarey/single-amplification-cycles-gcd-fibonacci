# Exactly Three Single-Amplification Cycles in GCD-Augmented Fibonacci Recurrences

Reproduction code for

> Hawarey, M. (2026). Exactly Three Single-Amplification Cycles in GCD-Augmented Fibonacci Recurrences. *AIR Journal of Mathematics and Computational Sciences* (revised version AIR-2026-000962-V2, under review).

`reproduce.py` recomputes every numerical result of the paper, Tables 1 and 2, and Figure 1, using exact integer arithmetic. It prints the values in manuscript order, with each block labelled by the section, table, figure or proposition where the values appear, and writes `figure1.png` (300 dpi) and `numbers.json` (all quoted values).

```
pip install matplotlib
python reproduce.py
```

Runtime: about one minute.

## What the script reproduces

| Output block | Manuscript |
|---|---|
| Candidate states, q = 1..8 | Section 4, after Theorem 7 |
| Non-coprime candidates | Section 5, after Lemma 8 |
| The three single-amplification cycles | Table 1 |
| Boundary cases and published tables | Section 5, after Corollary 11 |
| Seeds p = 9, 15, 27; worked example | Proposition 12 |
| Exact four-step ratios before lock-in | Figure 1 caption |
| Census, box 120 and box 200, every cycle reached | Table 2 |
| Robustness (2,400 steps) | Section 5 |
| Growth of seeds with no detected cycle | Section 5 |
| Other odd multiples of 3 up to 49 | Section 5 |
| Values for the record | Section 5 |
| Cycles of periods 22 and 65 | Section 5 |
| Cycle equation and the bound G > φ^q | Proposition 13 |
| Figure 1 | Figure 1 |

## Versions

- **2.0.0** — revised manuscript (AIR-2026-000962-V2): output labels aligned with the manuscript, full census for both boxes with every cycle reached (including the period-65 cycle), growth of the seeds with no detected cycle, Figure 1 caption values, checks for Proposition 13.
- **Single_Amp_Cycle_GCD_Fibonacci** (first archived release, DOI [10.5281/zenodo.23265695](https://doi.org/10.5281/zenodo.23265695)) — code for the version first submitted.

## License

MIT (see `LICENSE`).
