# The Silent Group: Track 2 submission

`solver.py` is the file to submit. It defines `solve(oracle) -> {"order", "secret"}` and uses only the
standard library, `qiskit` and the provided `silentgroup` package.

## Contents

| File | What it is |
| --- | --- |
| `solver.py` | Final solver (submit this) |
| `results.txt` | Official checker output, 10 runs per public meeting |
| `OPEN_TRACK.md` | Open track: query-cost model compared with measured costs |
| `slides.pptx` | 8-minute presentation |
| `solver_v1.py` | First version, kept for the "what failed" part of the talk |

## Run

```bash
pip install -e <path to qff-hackathon>
python -m silentgroup.check solver.py --repeat 10
```

## Results (`results.txt`)

| Meeting | Correct | Mean queries | Mean score |
| --- | --- | --- | --- |
| public_1 (valid, ord ≈ 2^56) | 10/10 | ~28k (best 19.6k) | 1.32 |
| public_2 (forgery, ord ≈ 2^35, very smooth) | 10/10 | ~33k (best 18.9k) | 1.29 |

The previous version (`solver_v1.py`) was correct 7 times out of 10 at about 42–75k queries per run.

## Algorithm

Each shot samples one character χ of G. `oracle_gate(h, k)` adds phase k·χ(h) (in turns). Every gate in a
shot shares the same χ.

1. **Order.** A semiclassical phase estimation reads χ(g) = j/n to 126 bits. It uses one qubit that is
   measured, reset and corrected with `if_test` phase gates. Continued fractions then give d = n/gcd(j, n).
   An identity test checks g^d = e. If that fails, g^d has small order t = gcd(j, n), and two cheap reads
   of g^d recover t.
2. **Logs.** With n known, one shot reads χ(g), χ(A), χ(B) and χ(C) to bitlen(n)+7 bits and rounds each to a
   multiple of 1/n: k_g = j and k_X = x·j. Then x = k_X · k_g⁻¹ mod n. If gcd(k_g, n) = s > 1, the solver
   tries up to 4 candidates or combines two samples with CRT.
3. **Checks.** An identity test runs `H · oracle · H · measure` for 16 shots and accepts if at most one shot
   reads 1. If the element is e, the qubit reads 0 except for noise (~0.2%). Otherwise each shot reads 1 with
   probability exactly 1/2, so a false accept happens with probability 17/65536.
   - A spare factor in n always divides gcd(k_g, n). The solver tests g^(n/p) for those primes p, which
     catches noise-corrupted orders such as lcm(n, 2^13).
   - An element is called forged only after three independent samples fail, so noise is not mistaken
     for a forgery.
