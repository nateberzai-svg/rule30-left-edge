#!/usr/bin/env python3
"""Tests on the frozen universe's centre column b_B(d) = r_d(d) (rulebook B8), as written by
wall_recursion128 ... COLFILE (bit i = r_d(d) for d = D0 + 1 + i).
Usage: frozen_column.py COLFILE D0 [CENTRE.npz]   (CENTRE: the real centre column c(t), for the match test)"""
import sys, numpy as np

def main():
    col = np.unpackbits(np.fromfile(sys.argv[1], dtype=np.uint8), bitorder='little'); D0 = int(sys.argv[2]) + 1
    N = len(col); print(f"frozen column b_B(d) = r_d(d) for d = {D0:,} .. {D0 + N - 1:,} ({N:,} bits)")
    print(f"density of ones: {col.mean():.6f}; by block of 1e8: {[round(float(col[i:i+10**8].mean()), 5) for i in range(0, N, 10**8)]}")
    # match with the real centre column
    if len(sys.argv) > 3:
        c = np.load(sys.argv[3])['c'].astype(np.uint8); n = min(len(c) - D0, N)
        agree = float(np.mean(col[:n] == c[D0:D0 + n]))
        print(f"agreement with the real centre column on d = {D0:,} .. {D0 + n - 1:,}: {agree:.5f} (chance ½ ± {0.5 / np.sqrt(n):.5f}; z = {(agree - 0.5) * 2 * np.sqrt(n):.2f})")
    # eventual periodicity: no period p <= PMAX with preperiod before the tail window
    W = col[-10**8:].copy(); PMAX = 10**7
    Wp = np.packbits(W, bitorder='little').view(np.uint64) if len(W) % 64 == 0 else None
    found = []
    for p in range(1, PMAX + 1):
        # compare 64 bits at offset 0 and p (early exit), then the whole window if they agree
        if np.array_equal(W[:256], W[p:p + 256]):
            if np.array_equal(W[:len(W) - p], W[p:]): found.append(p)
    print(f"periods p <= {PMAX:,} on the last {len(W):,} bits: {found if found else 'none'}")
    # runs
    d = np.diff(col.astype(np.int8)); starts = np.flatnonzero(d != 0) + 1
    lens = np.diff(np.concatenate(([0], starts, [N]))); h = np.bincount(lens)
    print("run-length counts (both colours) vs geometric(½) expectation, L = 1..30:")
    tot = len(lens)
    for L in range(1, 31):
        exp = tot / 2 ** L
        if L <= 24 or h[L]: print(f"   L={L:2d}: {h[L]:>12,d}  expected {exp:>14,.1f}  ratio {h[L] / exp if exp else 0:.4f}")
    print(f"   longest run {len(h) - 1}; runs total {tot:,}")
    # block frequencies for k = 8, 12, 16 (contexts built by shifts on chunks)
    def contexts(k, chunk=10**7):
        out = []
        for i in range(0, N - k, chunk):
            seg = col[i:i + chunk + k].astype(np.uint32); m = len(seg) - k
            ctx = np.zeros(m, dtype=np.uint32)
            for j in range(k): ctx = (ctx << 1) | seg[j:j + m]
            out.append(np.bincount(ctx, minlength=2 ** k))
        return np.sum(out, axis=0)
    for k in (8, 12, 16):
        cnt = contexts(k); exp = cnt.sum() / 2 ** k; chi = float(((cnt - exp) ** 2 / exp).sum()); dof = 2 ** k - 1
        print(f"blocks of {k}: chi-square {chi:,.1f} on {dof:,} dof (mean {dof:,}, sd {np.sqrt(2 * dof):,.0f}); max |dev|/sd = {np.abs(cnt - exp).max() / np.sqrt(exp):.2f}")
    # predictor from the last k bits (the rulebook's B10 test)
    print("best-guess accuracy of the next bit from the last k bits (½ = no memory):")
    for k in (1, 2, 4, 8, 12, 16):
        cnt = np.zeros((2 ** k, 2), dtype=np.int64); tot_n = 0      # counts pooled over the whole sequence
        for i in range(0, N - k - 1, 10**7):
            seg = col[i:i + 10**7 + k + 1].astype(np.uint32); m = len(seg) - k - 1
            if m <= 0: break
            ctx = np.zeros(m, dtype=np.uint32)
            for j in range(k): ctx = (ctx << 1) | seg[j:j + m]
            nxt = seg[k:k + m]; cnt += np.bincount(ctx * 2 + nxt, minlength=2 ** (k + 1)).reshape(-1, 2); tot_n += m
        n_per = tot_n / 2 ** k
        print(f"   k={k:2d}: {int(cnt.max(axis=1).sum()) / tot_n:.6f}   (a fair coin scores about {0.5 + 0.399 / np.sqrt(n_per):.6f} in sample at {n_per:,.0f} samples per context)")
    # autocorrelation at lags 1..64
    def agree(l, chunk=2 * 10**8):
        s = 0
        for i in range(0, N - l, chunk):
            a = col[i:i + chunk]; b = col[i + l:i + l + len(a)]; m = min(len(a), len(b)); s += int(np.count_nonzero(a[:m] == b[:m]))
        return s / (N - l)
    ac = [2 * agree(l) - 1 for l in range(1, 129)]
    print(f"autocorrelation r(l) = 2 P(agree) - 1, lags 1..128: max |r| = {max(abs(a) for a in ac):.2e} (sd of a coin {1 / np.sqrt(N):.1e}); lags with |r| > 4 sd: {[(l + 1, round(a / (1 / np.sqrt(N)), 1)) for l, a in enumerate(ac) if abs(a) > 4 / np.sqrt(N)]}")

if __name__ == '__main__':
    main()
