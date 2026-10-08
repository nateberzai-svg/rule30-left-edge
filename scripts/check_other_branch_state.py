"""Does a given finite starting row take the other branch at the seed's fork after diagonal 53,207?

Usage:  python scripts/check_other_branch_state.py STATE.txt
STATE.txt is a row of 0s and 1s whose first cell is its leftmost black cell, for example the starting state published by
szymon-lania, https://gist.github.com/szymon-lania/2535d62468a7912137b5002d900ebdd6 (not redistributed here; its sha256
begins 54eb2f35abb9e3d2). About a minute and 120 MB. For that state it prints: both evolutions periodic on every
diagonal checked; frozen diagonal 53,207 in both; the walls equal up to a shift of rows on [52800, 53199]; 292 of the
301 diagonals in [53200, 53500] differ; no shift of rows matches the walls on [53300, 53500].

Diagonal frame from the left edge: D_t[d] = u^t_{d-t}, d = 0 is the left edge (always black).
Rule 30 u^{t+1}_x = u^t_{x-1} xor (u^t_x or u^t_{x+1}) reads D_{t+1}[d] = D_t[d-2] xor (D_t[d-1] or D_t[d]).
Diagonals d <= DMAX need only d <= DMAX, so a fixed array of DMAX+1 cells is exact for them.
"""
import sys, numpy as np

DMAX, T, KEEP = 53_600, 110_000, 1024

def run(row0):
    D = np.zeros(DMAX + 1, dtype=np.uint8)
    n = min(len(row0), DMAX + 1); D[:n] = row0[:n]
    tail = np.zeros((KEEP, DMAX + 1), dtype=np.uint8)
    for t in range(T):
        if t >= T - KEEP: tail[t - (T - KEEP)] = D
        N = np.empty_like(D)
        N[0] = D[0]                                  # D_t[-2] = D_t[-1] = 0
        N[1] = D[1] | D[0]                            # D_t[-1] xor (D_t[0] or D_t[1]) with D_t[-1] = 0
        N[2:] = D[:-2] ^ (D[1:-1] | D[2:])
        D = N
    return tail                                       # rows T-KEEP .. T-1

def labels(tail, lo, hi):
    """eventual pattern of each diagonal: (period, 32-row word starting at row T-32), or None if not periodic in the tail"""
    out = {}
    for d in range(lo, hi + 1):
        col = tail[:, d]
        for p in (1, 2, 4, 8, 16, 32):
            if np.array_equal(col[p:], col[:-p]): out[d] = (p, col[-32:].copy()); break
        else: out[d] = None
    return out

state = np.frombuffer(open(sys.argv[1], 'rb').read().strip(), dtype=np.uint8) - ord('0')
assert set(np.unique(state)) <= {0, 1} and state[0] == 1, 'expected a 0/1 row whose first cell is its leftmost black cell'
seed = np.array([1], dtype=np.uint8)

LO, HI = 52_800, 53_500
Ls, Lg = labels(run(seed), LO, HI), labels(run(state), LO, HI)
print('non-periodic in tail: seed', sum(v is None for v in Ls.values()), ' gist', sum(v is None for v in Lg.values()))
frozen_s = [d for d, v in Ls.items() if v and not v[1].any()]
frozen_g = [d for d, v in Lg.items() if v and not v[1].any()]
print('frozen diagonals in [%d,%d]: seed %s  gist %s' % (LO, HI, frozen_s, frozen_g))

def shift_ok(d, s):
    a, b = Ls[d], Lg[d]
    return a is not None and b is not None and a[0] == b[0] and np.array_equal(np.roll(a[1], s), b[1])

# the row translation that matches the two walls BEFORE the fork, then whether it still holds after it
pre = range(LO, 53_200)
good = [s for s in range(32) if all(shift_ok(d, s) for d in pre)]
print('translations matching every diagonal in [%d, 53199]: %s' % (LO, good))
if good:
    s = good[0]
    bad = [d for d in range(53_200, HI + 1) if not shift_ok(d, s)]
    print('with that translation, diagonals in [53200, %d] that DIFFER: %d; first few %s' % (HI, len(bad), bad[:12]))
    for d in range(53_205, 53_212):
        a, b = Ls[d], Lg[d]
        print('  d=%d  seed p=%s %s   gist p=%s %s' % (d, a and a[0], a and ''.join(map(str, a[1][:16])),
                                                    b and b[0], b and ''.join(map(str, np.roll(b[1], -s)[:16]))))

# can ANY row translation re-match the walls after the fork?
post = range(53_300, HI + 1)
any_s = [s for s in range(32) if all(shift_ok(d, s) for d in post)]
print('translations matching every diagonal in [53300, %d]: %s' % (HI, any_s))
