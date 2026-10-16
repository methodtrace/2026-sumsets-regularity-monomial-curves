# Sumsets regularity versus Castelnuovo-Mumford regularity: computations

This archive reproduces every computation in the paper

> P. Chaudhary and K. Khandelwal, *Sumsets regularity versus Castelnuovo-Mumford regularity for
> monomial curves*.

Archive DOI: [10.5281/zenodo.23194227](https://doi.org/10.5281/zenodo.23194227)

The theorems are proved in the paper. A few finite assertions (that no four-element set with
d < 131 has reg < sigma and exactly six have d = 131, Remark 4.7, the least member of the band family with
reg < sigma, and that these counterexamples do not depend on the reading of the structure decomposition) are established by these computations; the other statements listed in Table 1 of the
paper are checked on the finite ranges given there. The scripts use only the Python standard library.

## Requirements

Python 3.8 or later. No other software.

## Running

```
python3 run_all.py
```

This runs every script in `checks/`, writes each script's output to `results/<name>.txt`, and
compares the SHA-256 digest of every result file with `MANIFEST.sha256`. The result files depend
only on the code, not on the machine, so a correct run ends with

```
all results match MANIFEST.sha256 byte for byte
```

Progress and timing are printed to the terminal only. The full run takes about 4 to 8 minutes on a
current desktop computer; most of it is `theorem_a` (d <= 150) and `theorem_b` (d <= 131).
Each script can also be run on its own, for example `python3 checks/theorem_b.py 60` for a
smaller range. Do not run with `python3 -O`: the checks are `assert` statements, and `run_all.py`
refuses an optimized interpreter.

## What each script checks

| Script | Statements in the paper |
|---|---|
| `selftest.py` | the method of Section 6: tau, c1, c2 from the homogenized Apery set agree with direct sumsets (four elements, d <= 40; five elements, d <= 24) |
| `theorem_a.py` | Lemmas 3.2 (including the full formula in its proof) and 3.3, the non-ACM criterion nu < mu, Theorems 3.5 (Theorem A) and 3.6, Remark 3.7 (four elements, d <= 150) |
| `case3_identities.py` | Case 3 of the proof of Theorem 3.5 (exact polynomial identities) |
| `theorem_b.py` | Theorem 4.1, Proposition 4.2, (4.5), (4.6), Theorem 4.3 (Theorem B), Corollary 4.4, Remark 4.7 (including the six sets attaining its least value), and the six sets with reg < sigma at d = 131 (ACM four-element sets, d <= 131) |
| `family_p3.py` | Theorem 4.5 (Theorem C): all pairwise coprime I, J, K <= 11, and the odd subfamily |
| `sigma_literal.py` | the paragraph after Definition 2.2: every set with reg < sigma from `family_p3.py` and `band_family.py` also has reg < sigma under the literal reading of (1.1), with an empty central interval allowed |
| `defect.py` | Proposition 2.6 (four elements, d <= 60; five elements, d <= 50) |
| `band_family.py` | Theorem 5.3 (Theorem D) for coprime m, n <= 18, Corollary 5.4, and the least member of the family with reg < sigma (every coprime pair with d <= 101) |
| `gluing.py` | Theorem 5.1 (2) to (4) and Proposition 5.2, for the gluings gB' with d' added, of every four-element B' with d' <= 40 and g in {2, 3} prime to d' |

The library in `curves/` computes the homogenized Apery set level by level with integers used as
bitsets (`core.py`), the invariants tau, reg, c1, c2, N and the defect E from it, a direct sumset
computation used by the self-test, and the fan of a four-element set (`fan.py`).

## License

Code: Apache License 2.0 (`LICENSE`). Result files in `results/`: CC BY 4.0 (`LICENSE-DATA`). See `NOTICE`.

## Citation

See `CITATION.cff`.
