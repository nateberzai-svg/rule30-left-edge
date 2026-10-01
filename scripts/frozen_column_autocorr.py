#!/usr/bin/env python3
"""Autocorrelation of the frozen universe's centre column b_B(d) = r_d(d), as written by
wall_recursion128 ... COLFILE (bit i = r_d(d) for d = D0 + 1 + i).
For each lag k it prints r_k = mean of x_i x_{i+k} with x = 1 - 2 b, and z = r_k sqrt(N - k),
the number of standard deviations from zero for a fair coin.
Usage: frozen_column_autocorr.py COLFILE D0 [NBITS]"""
import math, sys, numpy as np

LAGS = list(range(1, 17)) + [24, 31, 32, 33, 48, 63, 64, 65, 96, 128]

def main():
    bits = np.unpackbits(np.fromfile(sys.argv[1], dtype=np.uint8), bitorder='little')
    D0 = int(sys.argv[2]) + 1
    if len(sys.argv) > 3: bits = bits[:int(sys.argv[3])]
    x = (1 - 2 * bits.astype(np.int8)).astype(np.int8); N = len(x)
    print(f"frozen column b_B(d) = r_d(d) for d = {D0:,} .. {D0 + N - 1:,} ({N:,} bits)")
    print(f"autocorrelation at {len(LAGS)} lags: r_k = mean x_i x_(i+k), x = 1 - 2 b; z = r_k sqrt(N - k)")
    for k in LAGS:
        s = 0
        for i in range(0, N - k, 50_000_000):
            j = min(i + 50_000_000, N - k)
            s += int(np.dot(x[i:j].astype(np.int32), x[i + k:j + k].astype(np.int32)))
        r = s / (N - k)
        print(f"  lag {k:3d}: r = {r:+.3e}  z = {r * math.sqrt(N - k):+.2f}")

if __name__ == '__main__':
    main()
