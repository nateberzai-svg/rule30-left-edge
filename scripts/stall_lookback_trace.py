#!/usr/bin/env python3
"""Proof trace of the lookback verdict at a frozen diagonal (10.237d D; the paper's Theorem 4.5).

Same propagation as stall_lookback.py, but every inferred cell records the constraint instance and the cells
it used; the dependency tree of the target cell D_d(M_{d-1}+1) is printed bottom-up as a numbered derivation.
Usage: python stall_lookback_trace.py LABELFILE d K M_{d-K},...,M_{d-1}   (front rows; any common offset works)
"""
import sys
def trace(r, M, d, K):
    js = list(range(d - K, d + 1)); t_lo = min(M[j] for j in js if j < d) - 1; t_hi = M[d - 1] + 1
    val, why = {}, {}
    for j in js:
        for t in range(t_lo, t_hi + 1):
            if j < d and t >= M[j] + 1:
                val[(j, t)] = r(j, t); why[(j, t)] = ('label', ())
            elif j < d and t == M[j] and (j - 1) in M and M[j] > M[j - 1]:
                val[(j, t)] = 1 - r(j, t); why[(j, t)] = ('defect', ())
    cons = [(j, t) for j in js[2:] for t in range(t_lo, t_hi)]
    changed = True
    while changed:
        changed = False
        for (j, t) in cons:
            Y, X, A, B = (j, t + 1), (j - 2, t), (j - 1, t), (j, t)
            y, x, a, b = (val.get(k) for k in (Y, X, A, B))
            new = []
            o = 1 if (a == 1 or b == 1) else (0 if (a == 0 and b == 0) else None)
            src_o = [k for k in (A, B) if val.get(k) == 1][:1] if o == 1 else [A, B]
            if o is not None:
                if y is None and x is not None: new.append((Y, x ^ o, [X] + src_o))
                if x is None and y is not None: new.append((X, y ^ o, [Y] + src_o))
            if y is not None and x is not None:
                oo = y ^ x
                if oo == 0:
                    if a is None: new.append((A, 0, [Y, X]))
                    if b is None: new.append((B, 0, [Y, X]))
                else:
                    if a == 0 and b is None: new.append((B, 1, [Y, X, A]))
                    if b == 0 and a is None: new.append((A, 1, [Y, X, B]))
            for k, v, src in new:
                if k not in val:
                    val[k] = v; why[k] = (f'rule at ({j},{t})', tuple(src)); changed = True
    return val, why
def show(val, why, target, d, r):
    order, seen = [], set()
    def visit(k):
        if k in seen: return
        seen.add(k)
        for s in why[k][1]: visit(s)
        order.append(k)
    visit(target)
    name = lambda k: f"D_{{d{k[0]-d:+d}}}({k[1]})" if k[0] != d else f"D_{{d}}({k[1]})"
    for n, k in enumerate(order, 1):
        kind, src = why[k]
        if kind == 'label': txt = f"= r_{{d{k[0]-d:+d}}}({k[1]}) (row above the front)"
        elif kind == 'defect': txt = f"= 1 - r_{{d{k[0]-d:+d}}}({k[1]}) (defect at the landing row)"
        else: txt = "from " + ', '.join(f"{name(s)}={val[s]}" for s in src) + f" by D_j(t+1) = D_(j-2)(t) xor (D_(j-1)(t) or D_j(t))"
        print(f"{n:2d}. {name(k)} = {val[k]}   {txt}")
if __name__ == '__main__':
    lab = {}
    for line in open(sys.argv[1]):
        j, h = line.split(); lab[int(j)] = int(h, 16)
    d = int(sys.argv[2]); K = int(sys.argv[3]); rows = [int(x) for x in sys.argv[4].split(',')]
    M = {d - K + i: rows[i] for i in range(K)}
    r = lambda j, t: (lab[j] >> (t % 128)) & 1
    val, why = trace(r, M, d, K)
    tgt = (d, M[d - 1] + 1)
    print(f"target D_d({M[d-1]+1}) = {val.get(tgt)}  ({'advance' if val.get(tgt)==1 else 'stall' if val.get(tgt)==0 else 'undecided'})")
    if tgt in val: show(val, why, tgt, d, r)
