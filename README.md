# The left edge of Rule 30 without Rule 30

Code, data and logs accompanying

> N. Berzai, *The left edge of Rule 30 without Rule 30: the stripe recursion, Wolfram's number 2,107,985,255, and the forks beyond it*, 30 September 2026 — `paper/wolfram_number_note.pdf`.

The note cites a longer paper, also here:

> N. Berzai, *The centre column of Rule 30: excluded periods, decimation of the Duhamel class, and the front of the diagonal defects*, 2026 — `paper/RULE30_PAPER.pdf` (95 pages; source `paper/RULE30_PAPER.tex`).

The note is self-contained; nothing it proves depends on the longer paper. Every number stated in the note is recomputed by `paper/verify_claims.py` from the files here. The longer paper's own scripts are not all in this repository; it names them by file. Every numbered
result in the note carries a status label (Proved / Exact / Measured / Known / Open); nothing is claimed beyond
its label.

## Prior public work

The eighth frozen diagonal, d = 1,420,878,968, and the first period-64 onsets on its two branches (2,107,985,255 and,
on the other side, 9,958,700,505) were, to our knowledge, first made public on 8 September 2026 in two independent
code repositories:

- [Dibujaron/rule30](https://github.com/Dibujaron/rule30) (Drew McPolstra), commit `94833d2`, 8 September 2026.
  It gives the eighth all-white diagonal at 1,420,878,968, with two complementary continuations. On one, period 64
  starts at exactly 2,107,985,255, which it identifies with NKS p. 871. On the other, the next all-white diagonals
  are 3,340,408,059 and 4,989,445,007, and there is no doubling below 5.16·10⁹.
- [Patto1155/rule30-foundry](https://github.com/Patto1155/rule30-foundry), commit `04f30f6` (pull request 32),
  merged 8 September 2026 from a run of about 30 August 2026. It gives the first branch point at 1,420,878,969 with
  continuations `0x93425cac` and `0x6cbda353`, doublings at 2,107,985,255 and 9,958,700,505, and the whole branch
  tree below 1.2·10¹⁰. On 19 August 2026 it had given the exact periods of the left diagonals for d < 10⁶.

We obtained the same integers independently on 24 September 2026 and learned of these repositories on
1 October 2026. We found the numbers in no paper, preprint or OEIS entry. We claim no priority for them.
The logs in this repository follow the branch `93425cac` only to its next frozen diagonal, 3,340,408,059; our run to
9,958,700,505 is not included here.

To our knowledge no public source determines which branch the pattern grown from one cell takes at 1,420,878,969.
Neither repository does: the first infers the branch from the match with NKS and lists it as not checked; the
second records the actual diagonal as unresolved. A gist by
[szymon-lania](https://gist.github.com/szymon-lania/2535d62468a7912137b5002d900ebdd6) (21 August 2026) gives a
starting state that takes the other branch at diagonal 53,209 (53,208 in the indexing used here); we checked
this by direct simulation on 7 October 2026: before the frozen diagonal 53,207 its left diagonals match the seed's
up to a shift of rows, after it they differ and no shift matches them. So a finite starting word can take the other
branch, and the paper's left-half universality conjecture is false as stated (it is marked refuted there). In a comment on Eric Rowland's video *The Hidden Structure of Rule 30* the author also reports that more
than 2¹⁷ random starting states all took the seed's branch there.

## Layout

| folder | contents |
|---|---|
| `paper/` | the note (`.tex`, `.pdf`), the longer paper (`RULE30_PAPER.tex`, `.pdf`), the note's figures and the script that draws them, the data the verifier reads (`data/labels_800k.npz`: the eventual patterns of the seed's first 800,000 diagonals; `data/centre_column_1130k.npz`: the centre column to row 1,130,000), the verifier and its last output |
| `scripts/` | the programs: `wall_recursion.c` / `wall_recursion128.c` (the stripe recursion with 64- and 128-bit patterns), `edge_machine.c` (the coin-driven edge machine), `label_distance.c` (the fork law), `stall_lookback_trace.py` / `stall_lookback.py` / `stall_census.py` (Theorem 6.5 and Remark 6.6), `frozen_column.py` / `frozen_column_autocorr.py` (Section 7), `verify_frontier_segments.py` (re-checks every step of a frontier log with an independently written program), `check_other_branch_state.py` (evolves a given starting row beside the single seed and compares their left diagonals across the fork after 53,207; the check behind the refuted universality conjecture) |
| `bin/` | Windows builds of the three C programs (see `bin/README.md`) |
| `results/wall_recursion_logs/` | the recursion's logs: the run from the doubling to 1.508·10¹¹, the second machine's run from there to 10¹² with its checkpoints, the recomputed steps, the timing |
| `results/edge_machine_runs/` | the twenty coin-driven runs across the fork (`br_s*`), the four period-64 runs (`p64_s*`), the stall tails, the input-bias sweep (`eps_*`), the depth sweep (`depth_W*`) |
| `results/fork_law/` | the distance counts at period 32 and period 64 |
| `results/fork_mechanism/` | the labels around the fork, the eight machines across the fork and across the doubling, the census of front histories |
| `results/frozen_column/` | the tests of the background column b_B (Section 7); the 10⁹-bit column itself (125 MB) is not included and is regenerated in 35 s, see below |
| `results/front_800k_derived.npz` | the front M_d and settling times τ_d of the seed to d = 800,000, from a direct simulation |
| `ledger/` | the two working-record entries the verifier cites for the edge machine's validations (written 24–25 September 2026; each opens with a note, added 7 October 2026, on the prior public work above) |

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
