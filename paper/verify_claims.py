#!/usr/bin/env python3
"""Recompute every number stated in paper/wolfram_number_note.tex from the data in this repository.

Each check prints PASS, FAIL or SKIP, the claim as the paper states it, and the value computed here.
Run from the repository root:   python paper/verify_claims.py [--quick]
--quick skips the depth-8 census (about a minute) and the recomputation of thirteen steps of the frontier run
(half a minute), and reads the autocorrelation of the 10^9-bit frozen column from
results/frozen_column/autocorrelation.log instead of recomputing it. The recomputation (about two minutes)
needs the column file, 125 MB, which is not in the repository; it is regenerated in about 35 s by
  wall_recursion128 88d1fb8488d1fb84 7cb53b6b7cb53b6b 90000 1000000000 "" 0 - results/frozen_column/frozen_column_90k_1e9.bin
(or put its path in the environment variable FROZEN_COLUMN).
Data: paper/data/labels_800k.npz (the labels r_d, d <= 800,000, from a direct simulation; bit t mod per of
colword[d] is r_d(t)), paper/data/centre_column_1130k.npz, results/front_800k_derived.npz (the fronts M_d),
and the logs under results/.
"""
import glob, math, os, re, shutil, subprocess, sys, tempfile
import numpy as np

_open = open
def open(f, mode='r', *a, **k):        # every text file here is UTF-8, whatever the platform's default
    if 'b' not in mode and not a and 'encoding' not in k:
        k['encoding'] = 'utf-8'
    return _open(f, mode, *a, **k)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QUICK = '--quick' in sys.argv
results = []

def check(ok, claim, value):
    tag = 'SKIP' if ok is None else ('PASS' if ok else 'FAIL')
    results.append(tag)
    print(f"[{tag}] {claim}\n        computed: {value}")

def path(*p): return os.path.join(ROOT, *p)
def text(*p): return open(path(*p)).read()

# ----------------------------------------------------------------------------------------------- data
L = np.load(path('paper', 'data', 'labels_800k.npz'))
per = [int(x) for x in L['per']]; cw = [int(x) for x in L['colword']]; tau = [int(x) for x in L['tau']]
M = [int(x) for x in np.load(path('results', 'front_800k_derived.npz'))['M']]
def r(d, t): return (cw[d] >> (t % per[d])) & 1
MASK = (1 << 128) - 1
def word(d):                                       # label d as a 128-bit periodic word (bit t = r_d(t))
    x = 0
    for t in range(128): x |= r(d, t) << t
    return x
def rotl(x, k): return ((x << k) | (x >> (128 - k))) & MASK
def step(a, b):
    """next label from (r_{d-1}, r_d); two candidates at a frozen diagonal"""
    if b == 0:
        X = 0; s = 0
        for t in range(128):
            if s: X |= 1 << t
            s ^= (a >> t) & 1
        return (X, ~X & MASK)
    z = ~b & MASK; n = 0
    while z: z &= rotl(z, 1); n += 1
    c = 0
    for _ in range(n + 1): c = rotl(a ^ b ^ (~b & MASK & c), 1)
    assert rotl(a ^ b ^ (~b & MASK & c), 1) == c
    return (c,)
def weight(x, p): return bin(x & ((1 << p) - 1)).count('1')
def period(x):
    for p in (1, 2, 4, 8, 16, 32, 64, 128):
        if rotl(x, p) == x: return p

# ------------------------------------------------------------------ Section 3: validation of the recursion
a, b = word(19), word(20); ok = 0; bad = 0; frozen_seen = []
for d in range(20, 800000):
    out = step(a, b); nxt = word(d + 1)
    if len(out) == 2:
        frozen_seen.append(d)
        if nxt not in out: bad += 1
        c = nxt
    else:
        c = out[0]
        if c == nxt: ok += 1
        else: bad += 1
    a, b = b, c
check(bad == 0 and ok + len(frozen_seen) == 799980,
      "the recursion from (r_19, r_20), with the seed's branch at each frozen diagonal, reproduces the labels at all 799,980 diagonals to 800,000",
      f"{ok} deterministic steps agree, {len(frozen_seen)} frozen diagonals {frozen_seen} each offer the seed's label, {bad} disagreements")

# ------------------------------------------------------------------ Table 1, Corollary 3.2, Theorem 4.1
frozen = [d for d in range(1, 800001) if cw[d] == 0]
rows = [(d, per[d - 1], weight(cw[d - 1], per[d - 1])) for d in frozen]
check([x[0] for x in rows] == [2, 7, 28, 399, 53207, 58286, 87866]
      and [(x[1], x[2]) for x in rows] == [(1, 1), (2, 1), (4, 3), (8, 3), (16, 6), (16, 10), (16, 5)],
      "Table 1 to 800,000: frozen diagonals 2, 7, 28, 399, 53,207, 58,286, 87,866 with (period, weight) of r_{d-1} = (1,1) (2,1) (4,3) (8,3) (16,6) (16,10) (16,5)", rows)
runmax = []; mx = 0
for d in range(0, 800001):
    if per[d] > mx: mx = per[d]; runmax.append((d, mx))
check([x[0] for x in runmax[1:]] == [3, 8, 29, 400, 87867],
      "the running maximum of p_d first reaches 2, 4, 8, 16, 32 at d = 3, 8, 29, 400, 87,867", runmax)
twins = all(cw[d - 1] == cw[d - 2] and per[d - 1] == per[d - 2] for d in frozen if d >= 2)
non_frozen_twins = [d for d in range(2, 800001) if cw[d] != 0 and cw[d - 1] == cw[d - 2] and per[d - 1] == per[d - 2]]
check(twins and not non_frozen_twins, "Corollary 3.2(i) on the data: frozen exactly where r_{d-2} = r_{d-1}",
      f"frozen diagonals all twins: {twins}; twins that are not frozen: {non_frozen_twins[:5]}")

logs = {f: text('results', 'wall_recursion_logs', f) for f in ('wall_recursion_3e9.txt', 'wall_recursion_b0.txt', 'wall_recursion_b1.txt')}
allog = '\n'.join(logs.values())
fz = re.findall(r'd=(\d+): frozen diagonal \(r_d = 0\), p_\(d-1\)=(\d+), weight (\d+), candidates ([0-9a-f]+) / ([0-9a-f]+)', allog)
fz = sorted(set((int(d), int(p), int(w), x[:8], y[:8]) for d, p, w, x, y in fz))
check((1420878968, 32, 20, '93425cac', '6cbda353') in fz and (2107985254, 32, 17, 'b6e14797', '491eb868') in fz
      and (3340408059, 32, 16, '743a866e', '8bc57991') in fz,
      "Theorem 4.1: d* = 1,420,878,968 (period 32, weight 20, candidates 93425cac / 6cbda353); 2,107,985,254 (weight 17); 3,340,408,059 (weight 16)", fz)
b1 = logs['wall_recursion_b1.txt']; b0 = logs['wall_recursion_b0.txt']
first_frozen_after = re.findall(r'd=(\d+): frozen', logs['wall_recursion_3e9.txt'])
check('branch 1 taken' in b1 and '2107985254' in b1 and 'branch 0 taken' in b0 and '3340408059' in b0
      and first_frozen_after and int(first_frozen_after[0]) == 1420878968
      and not re.search(r'running max of the period reaches 64', b0),
      "from 90,000 the first frozen diagonal is d*; branch 6cbda353 next freezes at 2,107,985,254, branch 93425cac at 3,340,408,059 with period <= 32 before it",
      f"first frozen after 90,000: {first_frozen_after[:1]}; branch-0 log period histogram: {re.findall(r'period histogram.*', b0)}")
v = text('results', 'fork_law', 'label_distance_period64_1e10.log')
check('d=2107985255: running max of the period reaches 64' in v,
      "on branch 6cbda353 the period first reaches 64 at d = 2,107,985,255", re.findall(r'd=\d+: running max[^\n]*', v))
def rots(x, y, p):
    m = (1 << p) - 1; return [s for s in range(p) if ((x >> s) | (x << (p - s))) & m == y]
cand = {}
for d in (53207, 58286):
    X = step(word(d - 1), 0)[0]; p = per[d - 1]
    cand[d] = (weight(X, p), p - weight(X, p), rots(X & ((1 << p) - 1), (~X) & ((1 << p) - 1), p))
cand[1420878968] = (weight(0x93425cac, 32), weight(0x6cbda353, 32), rots(0x93425cac, 0x6cbda353, 32))
cand[3340408059] = (weight(0x743a866e, 32), weight(0x8bc57991, 32), rots(0x743a866e, 0x8bc57991, 32))
check(cand[53207][:2] == (4, 12) and cand[58286][:2] == (6, 10) and cand[1420878968][:2] == (14, 18)
      and cand[3340408059][:2] == (16, 16) and all(not v[2] for v in cand.values()),
      "Corollary 3.2(iii): candidate weights 4/12, 6/10, 14/18, and 16/16 at 3,340,408,059; no rotation maps one candidate to the other",
      cand)

# ---------------------------------------------------------------------------------- Section 5: the fork law
def hist(f):
    s = text('results', 'fork_law', f)
    m = re.search(r'distance histogram at period (\d+) \((\d+) pairs\):([^\n]*)', s)
    return int(m.group(1)), int(m.group(2)), {int(k): int(x) for k, x in re.findall(r'(\d+):(\d+)', m.group(3))}
P, n, h = hist('label_distance_period32.log')
e = [n * math.comb(P, k) / 2 ** P for k in range(P + 1)]
mid = [k for k in range(P + 1) if e[k] >= 5]
chi = sum((h[k] - e[k]) ** 2 / e[k] for k in mid)
check(n == 1420788969 and abs(chi - 29.7) < 0.05 and mid == list(range(1, 32)) and h[0] == 1 and h[32] == 0 and round(e[0], 2) == 0.33,
      "Law 5.1, P = 32: 1.42*10^9 pairs, chi^2 = 29.7 on the 31 classes with expectation >= 5 (distances 1-31); distances 0 and 32: 0.33 expected, 1 and 0 seen",
      f"{n} pairs, chi^2 = {chi:.2f} on {len(mid)} classes {mid[0]}..{mid[-1]}; k=0: {h[0]} (E {e[0]:.2f}), k=32: {h[32]}")
P, n, h = hist('label_distance_period64_1e10.log')
e = [n * math.comb(P, k) / 2 ** P for k in range(P + 1)]
mid = [k for k in range(P + 1) if e[k] >= 5]; chi = sum((h[k] - e[k]) ** 2 / e[k] for k in mid)
check(n == 10 ** 10 and abs(chi - 51.7) < 0.05 and len(mid) == 47,
      "Law 5.1, P = 64: chi^2 = 51.7 on the 47 classes with expectation >= 5",
      f"chi^2 = {chi:.1f} on {len(mid)} classes; distance 8: {e[8]:.1f} expected, {h[8]} seen; 56: {e[56]:.1f}/{h[56]}; 57: {e[57]:.2f}/{h[57]}")
check(round(e[8], 1) == 2.4 and h[8] == 3 and round(e[56], 1) == 2.4 and h[56] == 1 and round(e[57], 1) == 0.3 and h[57] == 1,
      "tails at P = 64: distance 8 (2.4 expected, 3 seen), 56 and 57 (2.4 and 0.3 expected, 1 and 1 seen)",
      [(k, round(e[k], 2), h[k]) for k in (8, 56, 57)])
check(87866 - 400 + 1 == 87467 and abs(2 ** 17 - 1.3e5) < 2e3 and abs(2 ** 33 / 8.6e9 - 1) < 0.01 and abs(2 ** 65 / 3.7e19 - 1) < 0.01,
      "regime lengths: period 16 from 400 to 87,866 is 87,467 long; 2^(P+1) = 1.3e5, 8.6e9, 3.7e19",
      (87866 - 400 + 1, 2 ** 17, f"{2**33:.3g}", f"{2**65:.3g}"))

# ------------------------------------------------------------------------ Section 6: branch and advance
seed_check = []
for d in (2, 7, 28, 399, 53207, 58286, 87866):
    m = M[d - 1]
    if M[d] > m: case, row = 'advance', M[d] + 1
    elif d >= 2 and M[d - 1] > M[d - 2]: case, row = 'stall after a move', m
    else: case, row = 'stall after a stall', None
    seed_check.append((d, case, None if row is None else r(d + 1, row) == 0))
check([x[2] for x in seed_check] == [True, None, True, True, True, True, True]
      and [x[1] for x in seed_check].count('advance') == 3,
      "Check after Theorem 6.2: advances at 2, 28, 58,286 and stalls right after a move at 399, 53,207, 87,866, each with the predicted branch; d = 7 in neither case",
      seed_check)
check(tau[7] == 0, "d = 7 is white from row 0", f"tau_7 = {tau[7]}")

# phase congruence: the landing rule shifted by s
def census(d_lo, d_hi, P):
    moves = [j for j in range(d_lo, d_hi + 1) if M[j] > M[j - 1] and per[j - 1] <= P]
    hits = [0] * P
    for j in moves:
        for s in range(P):
            t = M[j - 1] + s + 1
            while r(j - 1, t) == 0: t += 1
            if t == M[j] + s: hits[s] += 1
    return len(moves), hits
n16, h16 = census(20001, 87866, 16)
n32, h32 = census(90000, 800000, 32)
check(n16 == 37086 and h16[0] == n16 and round(max(h16[1:]) / n16, 3) == 0.261 and n32 == 385784 and h32[0] == n32 and round(max(h32[1:]) / n32, 3) == 0.261,
      "Proposition 6.3: shifted landing rule holds for s = 0 at all 37,086 moves between 20,000 and 87,866 and all 385,784 moves between 90,000 and 800,000, for no other residue; best other 26.1% in both",
      f"period 16: {n16} moves, s=0 {h16[0]}, best other {max(h16[1:])/n16:.3f}; period 32: {n32} moves, s=0 {h32[0]}, best other {max(h32[1:])/n32:.3f}")

# Theorem 6.5
lab = {}
for line in open(path('results', 'fork_mechanism', 'labels_dstar_minus13_to_dstar.txt')):
    j, hx = line.split(); lab[int(j)] = int(hx, 16)
ds = 1420878968; bit = lambda j, t: (lab[j] >> (t % 128)) & 1
check(lab[ds - 1] == lab[ds - 2] and lab[ds] == 0 and bit(ds - 3, 30) == 0 and bit(ds - 3, 31) == 1 and bit(ds - 2, 31) == 1
      and bit(ds - 1, 0) == 0 and bit(ds - 1, 1) == 1 and (lab[ds - 3] & 0xffffffff) == 0xb792cb87 and (lab[ds - 2] & 0xffffffff) == 0xdae372fa
      and (0x6cbda353 >> 2) & 1 == 0 and (0x93425cac >> 2) & 1 == 1,
      "Theorem 6.5: r_{d-3} = b792cb87 with bits 30, 31 = 0, 1; r_{d-2} = r_{d-1} = dae372fa with bit 31 = 1, bits 0, 1 = 0, 1; bit 2 of 6cbda353 is 0, of 93425cac is 1",
      (hex(lab[ds - 3] & 0xffffffff), hex(lab[ds - 2] & 0xffffffff), bit(ds - 3, 30), bit(ds - 3, 31), bit(ds - 2, 31), bit(ds - 1, 0), bit(ds - 1, 1)))
tr = subprocess.run([sys.executable, path('scripts', 'stall_lookback_trace.py'), path('results', 'fork_mechanism', 'labels_dstar_minus13_to_dstar.txt'),
                     str(ds), '3', '29,31,31'], capture_output=True, text=True).stdout
check('target D_d(32) = 1  (advance)' in tr and len([l for l in tr.splitlines() if re.match(r'\s*\d+\.', l)]) == 5,
      "the propagation trace of Theorem 6.5 is the five-line derivation", tr.strip().splitlines()[0])

# the twenty runs and the eight machines
brs = sorted(glob.glob(path('results', 'edge_machine_runs', 'br_s*.txt')))
br = [re.search(r'd=1420878968: frozen diagonal.*branch (\d) taken \(M_d=(\d+)', open(f).read()) for f in brs]
check(len(brs) == 20 and all(x and x.group(1) == '1' and int(x.group(2)) % 32 == 1 for x in br),
      "Computation 6.4: twenty runs from 90,000 all take branch 6cbda353, with M_{d*} = 1 mod 32",
      sorted(set((x.group(1), int(x.group(2)) % 32) for x in br if x)))
mf = text('results', 'fork_mechanism', 'machines_across_the_fork.log')
hists = re.findall(r'M_\(d\*-12\.\.d\*\+3\) = \[([^\]]+)\]\s+mod 32: \[([^\]]+)\]', mf)
mods = set(h[1] for h in hists); Ms = [int(h[0].split(',')[12]) for h in hists]
check(len(hists) == 8 and len(mods) == 1 and max(Ms) - min(Ms) == 2336 and list(map(int, list(mods)[0].split(','))) == [19, 19, 21, 21, 24, 24, 27, 27, 29, 29, 31, 31, 1, 1, 4, 5],
      "eight machines: one front history mod 32 over d*-12..d*+3; absolute rows differ by up to 2,336", (mods, max(Ms) - min(Ms)))

# census (Remark 6.6)
if QUICK:
    check(None, "Remark 6.6 census (depth 8 at d*: 1,038 / 1,791 / 1,267 of 4,096)", "skipped (--quick)")
else:
    out = subprocess.run([sys.executable, path('scripts', 'stall_census.py'), path('results', 'fork_mechanism', 'labels_dstar_minus13_to_dstar.txt'), str(ds), '8'],
                         capture_output=True, text=True, cwd=path('scripts')).stdout
    check('4096 histories: forced advance 1038, forced stall 1791, undecided 1267' in out,
          "Remark 6.6: at depth 8, 1,038 advances, 1,791 stalls, 1,267 undecided of 4,096", out.splitlines()[0] if out else out)
cs = text('results', 'fork_mechanism', 'census_seed_frozen_K5.log')
both = re.findall(r'd=(\d+) K=5: \d+ histories: forced advance (\d+), forced stall (\d+)', cs)
check(sorted(int(d) for d, a, s in both) == [28, 399, 53207, 58286, 87866] and all(int(a) > 0 and int(s) > 0 for d, a, s in both),
      "Remark 6.6: both verdicts occur at depth 5 at every frozen diagonal from 28 to 87,866", both)

# the doubling
md = text('results', 'fork_mechanism', 'machines_across_the_doubling.log')
brd = re.findall(r'branch (\d) taken \(M_d=(\d+)', md)
hd = re.findall(r'M = \[([^\]]+)\]\s+mod 32: \[([^\]]+)\]\s+M_DD mod 64 = (\d+)', md)
stalls = all(int(x[0].split(',')[8]) == int(x[0].split(',')[7]) for x in hd)
m64 = sorted(int(x[2]) for x in hd)
check(len(brd) == 8 and sorted(x[0] for x in brd).count('0') == 4 and stalls and m64.count(8) + m64.count(40) == 7 and 26 in m64
      and 'start DD-2e6' in md,
      "Remark 6.7: eight machines started 2*10^6 before it all stall at the doubling, branches four and four, seven at front row 8 mod 32 (8 or 40 mod 64), the eighth at 26",
      (sorted(x[0] for x in brd), m64, 'all stall' if stalls else 'not all stall'))

# ------------------------------------------------------------------ Section 7: the frozen column
ft = text('results', 'frozen_column', 'frozen_column_tests.log')
dens = float(re.search(r'density of ones: ([0-9.]+)', ft).group(1))
agr = re.search(r'agreement with the real centre column on d = 90,001 \.\. 1,130,000: ([0-9.]+) \(chance ½ ± [0-9.]+; z = (-?[0-9.]+)\)', ft)
runs = [(int(k), int(o.replace(',', '')), float(x.replace(',', ''))) for k, o, x in re.findall(r'L=\s*(\d+):\s+([\d,]+)\s+expected\s+([\d,.]+)', ft)]
zruns = [(k, (o - x) / math.sqrt(x)) for k, o, x in runs if k <= 12]
blocks = [(int(b), float(c.replace(',', '')), float(m.replace(',', '')), float(s)) for b, c, m, s in
          re.findall(r'blocks of (\d+): chi-square ([\d,.]+) on [\d,]+ dof \(mean ([\d,]+), sd (\d+)\)', ft)]
zb = [(b, (c - m) / s) for b, c, m, s in blocks]
check(dens == 0.500029 and agr and agr.group(1) == '0.49973' and agr.group(2) == '-0.55' and 'none' in re.search(r'periods p <= 10,000,000[^\n]*', ft).group(0),
      "Section 7: density 0.500029; agreement 0.49973 (z = -0.55) on rows 90,001-1,130,000; no period up to 10^7 on the last 10^8 bits",
      (dens, agr.group(1) if agr else None))
check(max(abs(z) for k, z in zruns) < 1.3 and [round(z, 1) for b, z in zb] == [-1.5, -0.3, -0.2],
      "Section 7: run counts for lengths 1-12 within 1.3 sd of a coin's; block chi-squares of 8, 12, 16 at 1.5, 0.3, 0.2 sd below the mean",
      ([(k, round(z, 2)) for k, z in zruns], [(b, round(z, 2)) for b, z in zb]))
cm = text('results', 'frozen_column', 'column_memory.log')
row16 = re.search(r'\n\s*16 \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+)', cm)
check(row16 and abs(float(row16.group(1)) - float(row16.group(2))) < 1e-4 and abs(float(row16.group(3)) - float(row16.group(5))) < 1e-3,
      "Section 7: the best guess from 16 bits scores as a coin of the same length (10^9 bits; and 1,040,000 bits with the real column)",
      row16.groups() if row16 else None)
AC_CLAIM = "Section 7: autocorrelation at the 26 lags 1-16, 24, 31-33, 48, 63-65, 96, 128 within 2.1 sd except lag 64, at -6.2 sd"
def ac_ok(out):
    others = [abs(z) for k, z in out.items() if k != 64]
    return len(out) == 26 and max(others) <= 2.1 and abs(out[64] + 6.2) <= 0.05 + 1e-9      # -6.2 to one decimal
colfile = os.environ.get('FROZEN_COLUMN', path('results', 'frozen_column', 'frozen_column_90k_1e9.bin'))
if QUICK or not os.path.exists(colfile):                   # the logged values (scripts/frozen_column_autocorr.py)
    out = {int(k): float(z) for k, z in re.findall(r'lag +(\d+): r = \S+ +z = ([-+0-9.]+)', text('results', 'frozen_column', 'autocorrelation.log'))}
    check(ac_ok(out), AC_CLAIM + " [from results/frozen_column/autocorrelation.log]", out)
else:                                                      # recomputed from the column itself
    bits = np.unpackbits(np.fromfile(colfile, dtype=np.uint8), bitorder='little')[:999910000].astype(np.int8)
    x = 1 - 2 * bits; N = len(x); out = {}
    for lag in list(range(1, 17)) + [24, 31, 32, 33, 48, 63, 64, 65, 96, 128]:
        s = 0
        for i in range(0, N - lag, 50_000_000):
            j = min(i + 50_000_000, N - lag); s += int(np.dot(x[i:j].astype(np.int32), x[i + lag:j + lag].astype(np.int32)))
        out[lag] = s / (N - lag) * math.sqrt(N - lag)
    check(ac_ok(out), AC_CLAIM + " [recomputed from the column]", {k: round(z, 2) for k, z in out.items()})
dd = text('results', 'frozen_column', 'defect_diagonal.log')
pc = re.search(r'P\(c=1\|b_B=0\) = ([0-9.]+), P\(c=1\|b_B=1\) = ([0-9.]+)', dd)
phase = [float(z) for z in re.findall(r'P\(Delta=1 \| t mod \d+\): range [0-9.]+\.\.[0-9.]+, max \|dev\| from ½ = ([0-9.]+) sd', dd)]
condvals = [float(v) for v in re.findall(r'(?:wall cells[^\n]*|p_t = \d+\) = |move\) = |wait\) = |P\(Delta=1\) = )([0-9]\.[0-9]{4})', dd)]
ctx = [float(v) for v in re.findall(r'\d{3}: ([0-9.]+) \(n=', dd)]
allc = ctx + [float(v) for v in re.findall(r'p_t = \d+\) = ([0-9.]+)', dd)] + [float(v) for v in re.findall(r'(?:move|wait)\) = ([0-9.]+)', dd)] \
       + [float(v) for v in re.findall(r'(?:below the front in \[[\d,]+\]|mode \w+ at t|column): P\(Delta=1\) = ([0-9.]+)', dd)]
check('density of Delta: 0.50015' in dd and pc and abs(float(pc.group(1)) - 0.5) < 1e-3 and abs(float(pc.group(2)) - 0.5) < 1e-3
      and max(phase) <= 2.1,
      "Proposition 7.1: density 0.50015; 2x2 table flat to 10^-3; phase residues within 2.1 sd",
      (pc.groups() if pc else None, phase, (min(allc), max(allc), len(allc))))

# ------------------------------------------------------------------ Section 8 and Figure 1
sp = [float(re.search(r'diagonals \d+: speed ([0-9.]+)', open(f).read()).group(1)) for f in brs]
mean, se = np.mean(sp), np.std(sp, ddof=1) / math.sqrt(len(sp))
st32 = text('results', 'edge_machine_runs', 'stall_tail_period32.txt'); st64 = text('results', 'edge_machine_runs', 'stall_tail_period64.txt')
g = lambda s, pat: re.search(pat, s).group(1)
sd = np.std(sp, ddof=1)
check(abs(mean - 1.321294) < 1e-6 and abs(sd - 0.000019) < 1.5e-6 and abs(se - 0.000004) < 1e-6 and '28,416,200,000 diagonals' in st32,
      "Section 8: speed 1.321294, mean of twenty runs over 2.8*10^10 diagonals; spread 0.000019, error of the mean 0.000004", f"{mean:.6f}, spread {sd:.6f}, error {se:.6f}")
p64 = sorted(glob.glob(path('results', 'edge_machine_runs', 'p64_s*.txt')))
sp64 = [float(re.findall(r'diagonals \d+: speed ([0-9.]+)', open(f).read())[-1]) for f in p64]
check(abs(np.mean(sp64) - 1.321438) < 1e-6 and abs(np.std(sp64, ddof=1) - 0.000007) < 1.5e-6, "Section 8: 1.321438 on the period-64 wall, four runs, spread 0.000007",
      f"{np.mean(sp64):.6f}, spread {np.std(sp64, ddof=1):.6f} over {len(sp64)} runs")
check(g(st32, r'P\(k\+1\)/P\(k\) = ([0-9.]+)') == '0.433' and g(st64, r'P\(k\+1\)/P\(k\) = ([0-9.]+)') == '0.434'
      and g(st32, r'longest stall seen (\d+)') == '24' and g(st64, r'longest stall seen (\d+)') == '26',
      "Section 8: stall tail ratios 0.433 and 0.434, longest stalls 24 and 26",
      (g(st32, r'P\(k\+1\)/P\(k\) = ([0-9.]+)'), g(st64, r'P\(k\+1\)/P\(k\) = ([0-9.]+)'), g(st32, r'longest stall seen (\d+)'), g(st64, r'longest stall seen (\d+)')))
check(abs(float(g(st32, r'after [0-9.e+]+ cycles = ([0-9.e+]+) diagonals')) - 2.7e12) < 0.1e12 and 0.5e24 < float(g(st64, r'after [0-9.e+]+ cycles = ([0-9.e+]+) diagonals')) < 2e24,
      "Section 8: a stall of 32 once per 2.7*10^12 diagonals, of 64 once per 10^24",
      (g(st32, r'after [0-9.e+]+ cycles = ([0-9.e+]+) diagonals'), g(st64, r'after [0-9.e+]+ cycles = ([0-9.e+]+) diagonals')))
ds_ = np.arange(1000, 800001); slope = np.polyfit(np.array(M)[ds_], (ds_ - np.array(M)[ds_]), 1)[0]
check(abs(slope + 0.24) < 0.01, "Figure 1: over the first 800,000 diagonals the front's average lean is 0.24 cells per row", f"slope {slope:.4f}")

# ------------------------------------------------------------------ checks added in the second referee pass
# Introduction and Section 7: the speed of the recursion
tm = text('results', 'wall_recursion_logs', 'timing.txt')
r128 = [float(x) for pair in re.findall(r'128-bit words: [^(]*\(([\d.]+), ([\d.]+) million', tm) for x in pair]
check(len(r128) == 4 and all(25 <= x <= 35 for x in r128) and 25 <= 1e9 / 35 / 1e6 <= 35,
      "Introduction: the recursion runs at about 3*10^7 diagonals a second on one core (128-bit words), consistent with 10^9 in 35 s",
      f"128-bit rates {r128} million/s; 10^9 in 35 s = {1e9 / 35 / 1e6:.1f} million/s")

# Theorem 4.3: the chain of logs from the doubling to 10^12
W = lambda f: text('results', 'wall_recursion_logs', f)
L1, L2, L3, L4, CK, B1 = (W('frontier_2.107e9_to_8.8e9.log'), W('frontier_8.6e9_to_1e11.log'), W('frontier_1e11_to_1.508e11.log'),
                          W('frontier_second_machine.txt'), W('frontier_checkpoint.txt'), W('wall_recursion_b1.txt'))
cps = lambda s_, pre='': {int(d): (a, b) for d, a, b in re.findall(r'^' + pre + r'checkpoint d=(\d+): a=([0-9a-f]+) b=([0-9a-f]+)', s_, re.M)}
start = lambda s_: re.search(r'start d=(\d+): a=([0-9a-f]+) \(p=\d+\) b=([0-9a-f]+)', s_)
c1, c3, nb, rc = cps(L1), cps(L3), cps(L4), cps(L4, r'#\s+')
s1, s2, s3 = start(L1), start(L2), start(L3)
fb1 = re.search(r'final labels: a=([0-9a-f]+) b=([0-9a-f]+)', B1); f2 = re.search(r'final labels: a=([0-9a-f]+) b=([0-9a-f]+)', L2)
h2 = re.search(r'period histogram \(p: count\): ([^\n]*)', L2).group(1); h2d = {int(k): int(v) for k, v in re.findall(r'(\d+): (\d+)', h2)}
ck = re.search(r'd=(\d+)\s+a=([0-9a-f]+)\s+b=([0-9a-f]+)', CK)
nb_body = '\n'.join(l for l in L4.splitlines() if not l.startswith('#'))
ok = (s1.group(1) == '2107985254' and fb1 and s1.group(2)[:16] == fb1.group(1) and int(s1.group(3), 16) == 0      # joins the branch-1 run
      and L1.count('frozen diagonal') == 1 and 'd=2107985255: running max of the period reaches 64' in L1
      and c1.get(int(s2.group(1))) == (s2.group(2), s2.group(3))                                                     # joins at 8.6e9
      and 1 not in h2d and 128 not in h2d and sum(h2d.values()) == 100000000000 - 8600000000 and 'frozen' not in L2
      and f2 and (f2.group(1), f2.group(2)) == (s3.group(2), s3.group(3)) and s3.group(1) == '100000000000'          # joins at 1e11
      and 'frozen' not in L3 and max(c3) == 150800000000 and c3[150800000000] == (ck.group(2), ck.group(3))         # the checkpoint
      and 'frozen' not in nb_body and 'running max' not in nb_body and 'labels for the machine' not in nb_body)
check(ok, "Theorem 4.3: logs join from the doubling to 1.508*10^11 with no frozen diagonal and no period 128, and the second machine's lines report no event",
      f"segments start at {s1.group(1)}, {s2.group(1)}, {s3.group(1)}; last cloud checkpoint {max(c3)}")
# the end of the second machine's log: the histogram's total is every diagonal of the run, so it never stopped early
hn = {int(k): int(v) for k, v in re.findall(r'(\d+): (\d+)', re.search(r'^period histogram \(p: count\):(.*)$', nb_body, re.M).group(1))}
fn = re.search(r'^final labels: a=([0-9a-f]+) b=([0-9a-f]+)', nb_body, re.M)
n_run = 10**12 - 150800000000; exp32 = n_run / 2**32
check(max(nb) == 10**12 and nb[10**12] == (fn.group(1), fn.group(2)) and sum(hn.values()) == n_run and set(hn) == {32, 64}
      and int(fn.group(2), 16) != 0 and abs(hn[32] - exp32) < 3 * math.sqrt(exp32) and n_run // 10**8 == 8492,
      "Theorem 4.3: the second machine's run ends at 10^12 having computed all 8,492 steps (849,200,000,000 labels): no frozen diagonal, no period 128; 188 labels of period 32 by chance against 197.7 expected",
      f"last checkpoint {max(nb)}; histogram {hn}; expected period-32 labels {exp32:.1f} +- {math.sqrt(exp32):.1f}")
# thirteen steps of 10^8 diagonals at six places, recomputed on the first machine (recorded in the file)
ST = [(d, d + 10**8) for d0, k in ((150800000000, 3), (304600000000, 2), (430400000000, 2), (543700000000, 2),
                                   (567600000000, 2), (650100000000, 2)) for d in range(d0, d0 + k * 10**8, 10**8)]
lab = lambda d: (ck.group(2), ck.group(3)) if d == int(ck.group(1)) else nb[d]
check(len(ST) == 13 and len({d0 // 10**10 for d0, _ in ST}) == 6 and max(d1 for _, d1 in ST) == 650300000000
      and all(nb[d1] == rc[d1] for _, d1 in ST),
      "Theorem 4.3: thirteen of the second machine's steps, at six places up to 6.503*10^11, equal the recomputation recorded on the first machine",
      [d1 for _, d1 in ST if nb[d1] == rc[d1]])
cc = shutil.which('gcc') or shutil.which('cc')
if QUICK or not cc:
    check(None, "Theorem 4.3: recompute the thirteen steps here with the 64-bit and the 128-bit program", "skipped (--quick or no C compiler)")
else:
    with tempfile.TemporaryDirectory() as tmp:
        ex64, ex128 = os.path.join(tmp, 'wr64'), os.path.join(tmp, 'wr128')
        subprocess.run([cc, '-O2', '-o', ex64, path('scripts', 'wall_recursion.c')], check=True)
        subprocess.run([cc, '-O2', '-o', ex128, path('scripts', 'wall_recursion128.c')], check=True)
        runs = []
        for d0, d1 in ST:                     # the 64-bit program reads the halves of the 128-bit words (period <= 64)
            a, b = lab(d0)
            assert all(w[:16] == w[16:] for w in (a, b, *nb[d1]))
            runs.append((d1, 16, subprocess.Popen([ex64, a[:16], b[:16], str(d0), str(d1), ''], stdout=subprocess.PIPE, text=True)))
            runs.append((d1, 32, subprocess.Popen([ex128, a, b, str(d0), str(d1), '', '128', '-'], stdout=subprocess.PIPE, text=True)))
        got = [(d1, n, re.search(r'final labels: a=([0-9a-f]+) b=([0-9a-f]+)', p.communicate()[0])) for d1, n, p in runs]
    agree = [d1 for d1, n, m in got if m and (m.group(1), m.group(2)) == (nb[d1][0][:n], nb[d1][1][:n])]
    check(len(agree) == 26, "Theorem 4.3: the thirteen steps recomputed here by both programs reach the second machine's checkpoints",
          f"{len(agree)} of 26 runs agree")

# Proposition 6.3 context and Computation 6.4: machine starts and the validations recorded in the ledger
led = lambda n: open(glob.glob(path('ledger', f'entry_10_{n}d_HELD.md'))[0]).read()
check(all(open(f).read().startswith('start d=90000 M=120085') for f in brs) and int(M[90000]) == 120085
      and 'M0 = 1.3213*D0 + 7*seed' in mf,
      "Computation 6.4: the twenty runs start at the seed's row M_90000 = 120,085; the eight start with rows 7 apart",
      (int(M[90000]), sorted(set(open(f).read()[:30] for f in brs))[:1]))
l31, l33 = led('231'), led('233')
check(re.search(r'zero\s+deviations over 5,000 diagonals at W = 32, and over 3,000 at W = 90 and at W = 4', l31) is not None
      and re.search(r"reproduces the Python machine's front with 0\s+mismatches over 200,000 diagonals", l33) is not None
      and re.search(r'in all twelve cases s_b ≡ s_w modulo the period', l31) is not None
      and '((01)^∞, (0011)^∞, coin, b ≡ 1, random states of\nwidth 40 and 100)' in l31,
      "Computation 6.4 and Proposition 6.3 [from the ledger]: true-bit runs without deviation over 3,000-5,000 diagonals at three depths; C = Python over 200,000 with constant inputs; six machines, twelve congruences",
      'found in 10.231d and 10.233d')

# Remark 6.6: the seed's own advances at 28 and 58,286 are not forced at any depth up to 32
lk = subprocess.run([sys.executable, path('scripts', 'stall_lookback.py')], capture_output=True, text=True, cwd=path('scripts')).stdout
und = {int(d): v for d, v in re.findall(r'd=(\d+): observed \w+; verdict by lookback depth K: (\[[^\n]*\])', lk)}
check(set(und) >= {28, 58286} and all(re.findall(r"'([^']+)'\)", und[d]) == ['?'] * 9 for d in (28, 58286)),
      "Remark 6.6: at 28 and 58,286 the seed's advance is not forced by the propagation at any depth up to 32",
      {d: und.get(d) for d in (28, 58286)})

# Section 7, Proposition 7.1 and Figure 1: tau and the front
bad_tau = [t for t in range(1, 800001) if tau[t] <= t]
check(max(bad_tau) == 19 and tau[2] == 2 and tau[7] == 0 and tau[28] == 31,
      "Section 7: tau_t > t for every t from 20 to 800,000; Figure 1: tau_2 = 2, tau_7 = 0, tau_28 = 31",
      f"t <= 800,000 with tau_t <= t: {bad_tau}; tau_2, tau_7, tau_28 = {tau[2]}, {tau[7]}, {tau[28]}")
fd = [d for d in range(12, 400) if 0 < M[d] < 360]
lean = np.polyfit([M[d] for d in fd], [d - M[d] for d in fd], 1)[0]
check(abs(lean + 0.2) < 0.02, "Figure 1: in rows 0-359 the front leans left by about 0.2 cells per row", f"slope {lean:.3f}")
zc = []
for cm_ in re.finditer(r'(\d{3}): ([0-9.]+) \(n=(\d+)\)', dd): zc.append((float(cm_.group(2)) - 0.5) / (0.5 / math.sqrt(int(cm_.group(3)))))
for cm_ in re.finditer(r'P\(Delta=1(?: \| (?:p_t = \d+|move|wait))?\) = ([0-9.]+) \(n=(\d+)\)', dd):
    zc.append((float(cm_.group(1)) - 0.5) / (0.5 / math.sqrt(int(cm_.group(2)))))
check(len(zc) >= 31 and max(abs(z) for z in zc) <= 2.1 and max(phase) <= 2.1,
      "Proposition 7.1: every conditional frequency of the transient bit within 2.1 sd of 1/2 (phase residues, wall cells, next black cell, period, move or wait, distance from the front, mode)",
      f"{len(zc)} logged conditionals, max |z| {max(abs(z) for z in zc):.2f}; phase residues max {max(phase)} sd")

# Section 8: coins against stretches of wall; the start-up stall excluded from the tails
Sb = np.array([[float(l.split()[4]) for l in open(f) if not l.startswith('#')][:14]
               for f in sorted(glob.glob(path('results', 'edge_machine_runs', 'br_s*.blocks')))])
coin_sd = Sb.std(axis=0, ddof=1); bm = Sb.mean(axis=0); ratio = bm.std(ddof=1) / (coin_sd.mean() / math.sqrt(Sb.shape[0]))
check(Sb.shape == (20, 14) and coin_sd.max() < 0.00011 and abs(bm.std(ddof=1) - 0.00013) < 0.000005 and 6.5 < ratio < 7.5,
      "Section 8: per block of 10^8 the twenty coins agree to 0.0001 (sd); block means vary with sd 0.00013, seven times the coins' share",
      f"max coin sd {coin_sd.max():.6f}; sd of block means {bm.std(ddof=1):.6f}; ratio {ratio:.1f}")
check('ignition delays excluded' in st32 and 'ignition delays excluded' in st64,
      "Section 8: the one stall at the start of each run is left out of the tails", re.findall(r'ignition delays excluded: \[[^\]]*\]', st32 + st64))

print(f"\n{results.count('PASS')} PASS, {results.count('FAIL')} FAIL, {results.count('SKIP')} SKIP")
