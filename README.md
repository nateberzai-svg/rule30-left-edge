# The left edge of Rule 30 without Rule 30

Code, data and logs accompanying

> N. Berzai, *The left edge of Rule 30 without Rule 30: the stripe recursion, Wolfram's number 2,107,985,255, and the forks beyond it*, 30 September 2026 — `paper/wolfram_number_note.pdf`.

Every number stated in the note is recomputed by `paper/verify_claims.py` from the files here. Every numbered
result in the note carries a status label (Proved / Exact / Measured / Known / Open); nothing is claimed beyond
its label.

## Layout

| folder | contents |
|---|---|
| `paper/` | the note (`.tex`, `.pdf`), its figures and the script that draws them, the data the verifier reads (`data/labels_800k.npz`: the eventual patterns of the seed's first 800,000 diagonals; `data/centre_column_1130k.npz`: the centre column to row 1,130,000), the verifier and its last output |
| `scripts/` | the programs: `wall_recursion.c` / `wall_recursion128.c` (the stripe recursion with 64- and 128-bit patterns), `edge_machine.c` (the coin-driven edge machine), `label_distance.c` (the fork law), `stall_lookback_trace.py` / `stall_lookback.py` / `stall_census.py` (Theorem 6.5 and Remark 6.6), `frozen_column.py` / `frozen_column_autocorr.py` (Section 7), `verify_frontier_segments.py` (re-checks every step of a frontier log with an independently written program) |
| `bin/` | Windows builds of the three C programs (see `bin/README.md`) |
| `results/wall_recursion_logs/` | the recursion's logs: the run from the doubling to 1.508·10¹¹, the second machine's run from there to 10¹² with its checkpoints, the recomputed steps, the timing |
| `results/edge_machine_runs/` | the twenty coin-driven runs across the fork (`br_s*`), the four period-64 runs (`p64_s*`), the stall tails, the input-bias sweep (`eps_*`), the depth sweep (`depth_W*`) |
| `results/fork_law/` | the distance counts at period 32 and period 64 |
| `results/fork_mechanism/` | the labels around the fork, the eight machines across the fork and across the doubling, the census of front histories |
| `results/frozen_column/` | the tests of the background column b_B (Section 7); the 10⁹-bit column itself (125 MB) is not included and is regenerated in 35 s, see below |
| `results/front_800k_derived.npz` | the front M_d and settling times τ_d of the seed to d = 800,000, from a direct simulation |
| `ledger/` | the two working-record entries the verifier cites for the edge machine's validations |

## Reproducing

Python 3 with `numpy` (and `matplotlib` for the figures); a C compiler (gcc) for the full run.

```
python paper/verify_claims.py --quick     # every check that needs no compiler: 43 PASS, 0 FAIL, 2 SKIP (about 6 minutes)
python paper/verify_claims.py             # also recompiles the thirteen frontier steps with both programs: 45 PASS
```

The C programs build with `gcc -O2 -o wall_recursion scripts/wall_recursion.c` (likewise the other two).
To regenerate the background column of Section 7 and rerun its tests:

```
./wall_recursion128 A_HEX B_HEX 90000 1000000000 "" 0 - results/frozen_column/frozen_column_90k_1e9.bin
python scripts/frozen_column.py results/frozen_column/frozen_column_90k_1e9.bin 90000 paper/data/centre_column_1130k.npz
python scripts/frozen_column_autocorr.py results/frozen_column/frozen_column_90k_1e9.bin 90000
```

with `A_HEX B_HEX` the labels at 89,999 and 90,000 as 128-bit hex words (`bin/README.md` gives the
64-bit pair `88d1fb8488d1fb84 7cb53b6b7cb53b6b`, to be repeated twice for 128 bits).
The frontier run from the checkpoint in `results/wall_recursion_logs/frontier_checkpoint.txt` to 10¹² takes
about eight hours on one core; the command is in `bin/README.md`.

## Licence

Mozilla Public License 2.0 (`LICENSE`).
