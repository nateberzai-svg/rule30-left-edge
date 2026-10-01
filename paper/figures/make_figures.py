#!/usr/bin/env python3
"""Figures of the Wolfram-number note (paper/wolfram_number_note.tex).

fig_left_edge.pdf   the left edge of Rule 30 from one cell: stripes, frozen diagonals, the front M_d
fig_fork_law.pdf    Hamming distance of consecutive labels against binomial(P, 1/2), P = 32 and 64
fig_fork_path.pdf   the front's path through the labels at the eighth frozen diagonal (1,420,878,968)
Data: results/front_800k_derived.npz, results/fork_law/*.log, results/fork_mechanism/*.
Palette: slots 1-2 of the validated categorical palette (blue #2a78d6, orange #eb6834); ink #0b0b0b / #52514e.
"""
import logging, os, re, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
logging.getLogger('fontTools').setLevel(logging.ERROR)     # quiet the TrueType-embedding timestamp notices
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
BLUE, ORANGE, INK, INK2, GRID = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#d9d8d4'
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 9, 'axes.edgecolor': INK2,
                     'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.6,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'legend.frameon': False,
                     'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})
CELLS = ListedColormap(['#ffffff', '#2b2b2b'])

def fig_left_edge(T=360, XL=-360, XR=24):
    W = 2 * T + 3; row = np.zeros(W, np.uint8); row[T + 1] = 1; img = np.zeros((T, XR - XL), np.uint8)
    for t in range(T):
        img[t] = row[T + 1 + XL: T + 1 + XR]
        row = np.roll(row, 1) ^ (row | np.roll(row, -1))
    M = np.load(os.path.join(ROOT, 'results', 'front_800k_derived.npz'))['M']
    ds = [d for d in range(12, 400) if 0 < M[d] < T]
    fx = [d - M[d] for d in ds]; ft = [M[d] for d in ds]
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.imshow(img, cmap=CELLS, interpolation='nearest', extent=(XL - 0.5, XR - 0.5, T - 0.5, -0.5), aspect='equal')
    ax.plot(fx, ft, color=BLUE, lw=1.6, label=r'front: the cell of diagonal $d$ on row $M_d$')
    for d, tau, ty, lx in ((7, 0, 150, -330), (28, 31, 230, -330)):   # frozen diagonals: all white from row tau on, along x = d - t
        ax.plot([d - tau, d - (T - 1)], [tau, T - 1], color=ORANGE, lw=0.8, ls=(0, (3, 2)))
        ax.annotate(rf'frozen diagonal $d={d}$', xy=(d - ty, ty), xytext=(lx, ty - 55), color=INK, fontsize=7.5,
                    arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5))
    ax.plot([], [], color=ORANGE, lw=0.8, ls=(0, (3, 2)), label='frozen diagonals (eventually all white)')
    ax.set_xlabel(r'cell $x$ (the seed is at $x=0$)'); ax.set_ylabel(r'row $t$')
    ax.set_xlim(XL - 0.5, XR - 0.5); ax.set_ylim(T - 0.5, -0.5)
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), fontsize=7.5, handlelength=2.2)
    fig.savefig(os.path.join(HERE, 'fig_left_edge.pdf'), metadata={'CreationDate': None}); fig.savefig(os.path.join(HERE, 'fig_left_edge.png'), dpi=160)
    plt.close(fig)

def read_hist(path):
    s = open(path).read(); m = re.search(r'distance histogram at period (\d+) \((\d+) pairs\):([^\n]*)', s)
    P, n = int(m.group(1)), int(m.group(2))
    h = {int(k): int(v) for k, v in re.findall(r'(\d+):(\d+)', m.group(3))}
    return P, n, h

def fig_fork_law():
    files = [os.path.join(ROOT, 'results', 'fork_law', f) for f in ('label_distance_period32.log', 'label_distance_period64_1e10.log')]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.7), sharey=True)
    for ax, f in zip(axes, files):
        P, n, h = read_hist(f)
        k = np.arange(P + 1); exp = np.array([n * math.comb(P, int(x)) / 2 ** P for x in k])
        ax.plot(k, exp, color=ORANGE, lw=1.2, label=r'$n\binom{P}{k}2^{-P}$', zorder=2)
        ko = [x for x in k if h.get(int(x), 0) > 0]; vo = [h[int(x)] for x in ko]
        ax.plot(ko, vo, ls='none', marker='o', ms=3.4, mfc=BLUE, mec='white', mew=0.5, color=BLUE, label='observed', zorder=3)
        ax.set_yscale('log'); ax.set_ylim(0.08, 3e10); ax.set_xlim(-1, P + 1)
        ax.grid(True, axis='y', color=GRID, lw=0.4); ax.set_axisbelow(True)
        ax.set_xlabel(r'Hamming distance $k$ between $r_{d}$ and $r_{d+1}$')
        ax.set_title(rf'$P={P}$, ' + (r'$n=1.42\cdot 10^{9}$ pairs' if P == 32 else r'$n=10^{10}$ pairs'), fontsize=8.5, color=INK)
        if P == 32:
            ax.annotate('the fork\n$d=1{,}420{,}878{,}968$', xy=(0, 1), xytext=(3.5, 30), fontsize=7.5, color=INK,
                        arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5))
    axes[0].set_ylabel('number of pairs'); axes[1].legend(loc='upper right', fontsize=7.5)
    fig.savefig(os.path.join(HERE, 'fig_fork_law.pdf'), metadata={'CreationDate': None}); fig.savefig(os.path.join(HERE, 'fig_fork_law.png'), dpi=160)
    plt.close(fig)

MASK = (1 << 128) - 1
def rotl(x, k): return ((x << k) | (x >> (128 - k))) & MASK
def step(a, b):
    if b == 0:
        X = 0; x = 0
        for t in range(128):
            if x: X |= 1 << t
            x ^= (a >> t) & 1
        return (X, ~X & MASK)
    z = ~b & MASK; L = 0
    while z: z &= rotl(z, 1); L += 1
    c = 0
    for _ in range(L + 1): c = rotl(a ^ b ^ (~b & MASK & c), 1)
    return (c,)

def fig_fork_path():
    dstar = 1420878968
    lab = {}
    for line in open(os.path.join(ROOT, 'results', 'fork_mechanism', 'labels_dstar_minus13_to_dstar.txt')):
        j, hx = line.split(); lab[int(j)] = int(hx, 16)
    X, nX = step(lab[dstar - 1], lab[dstar]); lab[dstar + 1] = nX          # branch 1: 6cbda353...
    a, b = lab[dstar], lab[dstar + 1]
    for j in range(dstar + 1, dstar + 3):
        c = step(a, b)[0]; lab[j + 1] = c; a, b = b, c
    phases = [19, 19, 21, 21, 24, 24, 27, 27, 29, 29, 31, 31, 1, 1, 4, 5]     # M_j mod 32, j = d*-12 .. d*+3 (8 of 8 machines)
    js = list(range(dstar - 12, dstar + 4)); un = []; off = 0
    for p in phases:                      # the front never moves back: unwrap to a non-decreasing row
        if un and p + off < un[-1]: off += 32
        un.append(p + off)
    lo, hi = 14, 42
    rows = list(range(dstar - 12, dstar + 4))
    img = np.array([[(lab[j] >> (t % 128)) & 1 for t in range(lo, hi)] for j in rows], dtype=np.uint8)
    fig, ax = plt.subplots(figsize=(6.2, 3.5))
    ax.imshow(img, cmap=ListedColormap(['#ffffff', '#b9b8b3']), interpolation='nearest',
              extent=(lo - 0.5, hi - 0.5, len(rows) - 0.5, -0.5), aspect='auto')
    for t in range(lo, hi + 1): ax.axvline(t - 0.5, color='#ecebe7', lw=0.3)
    for i in range(len(rows) + 1): ax.axhline(i - 0.5, color='#ecebe7', lw=0.3)
    y = [rows.index(j) for j in js]
    ax.plot(un, y, color=BLUE, lw=1.4, marker='o', ms=4.2, mfc=BLUE, mec='white', mew=0.6, label=r'front $M_j$ mod 32 (all eight machines)', zorder=3)
    i0 = rows.index(dstar)
    ax.plot([un[js.index(dstar)]], [i0], marker='o', ms=7, mfc='none', mec=ORANGE, mew=1.3, zorder=4)
    ax.add_patch(plt.Rectangle((un[js.index(dstar)] + 1 - 0.5, i0 + 1 - 0.5), 1, 1, fill=False, ec=ORANGE, lw=1.3, zorder=4))
    ax.annotate('frozen diagonal $d^*$ (all white):\nthe front advances two rows', xy=(un[js.index(dstar)] + 0.35, i0),
                xytext=(hi + 0.3, i0 - 2.5), fontsize=7.5, color=INK, annotation_clip=False,
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5))
    ax.annotate('branch taken: the label\nwith $r_{d^*+1}(M_{d^*}+1)=0$', xy=(un[js.index(dstar)] + 1.5, i0 + 1),
                xytext=(hi + 0.3, i0 + 2.2), fontsize=7.5, color=INK, annotation_clip=False,
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5))
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([rf'$d^*{j - dstar:+d}$' if j != dstar else r'$d^*$' for j in rows], fontsize=7)
    xt = list(range(16, hi, 4)); ax.set_xticks(xt); ax.set_xticklabels([str(t % 32) for t in xt])
    ax.set_xlabel(r'row $t$ modulo 32 (grey: $r_j(t)=1$)'); ax.set_ylabel(r'diagonal $j$')
    ax.legend(loc='upper left', fontsize=7.5, bbox_to_anchor=(0.0, 1.13), ncol=1)
    fig.savefig(os.path.join(HERE, 'fig_fork_path.pdf'), metadata={'CreationDate': None}); fig.savefig(os.path.join(HERE, 'fig_fork_path.png'), dpi=160)
    plt.close(fig)
    return un

if __name__ == '__main__':
    fig_left_edge(); fig_fork_law(); un = fig_fork_path(); print('unwrapped front phases', un)
