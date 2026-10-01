#!/usr/bin/env python3
"""Lookback by constraint propagation: is the stall or advance of a frozen diagonal forced by the wall and
the front's recent history alone?  (10.237d)

Frozen diagonal d (r_d = 0, hence r_{d-1} = r_{d-2}) after a move at d-1.  The front stalls at d iff
D_d(M_{d-1}) = 1 (then the label with r_{d+1}(M_{d-1}) = 0 is taken); else it advances to the next 1 of r_{d-1}
(and the label with r_{d+1}(M_d) = 1 is taken).  Known cells: D_j(t) = r_j(t) for t >= M_j + 1 (j < d), and
D_j(M_j) = 1 - r_j(M_j) where the front moved at j.  The frame recursion D_j(t+1) = D_{j-2}(t) xor (D_{j-1}(t)
or D_j(t)) is propagated as a constraint (arc consistency, both directions) over the cells of diagonals
d-K..d on the rows between M_{d-K} and M_{d-1}; the answer is D_d(M_{d-1}) if it becomes determined.
API: lookback(r, M, d, K) -> D_d(M_{d-1}+1): 1 (advance forced), 0 (stall forced), None (not decided by K fronts).
The front advances at d iff D_d(M_{d-1}+1) = 1 (for t >= M_{d-1}+1 the frozen diagonal keeps its value until the
next 1 of r_{d-1} and is 0 from there on); with a move at d-1 this is the same as D_d(M_{d-1}) = 0.
"""
def lookback(r, M, d, K):
    js=list(range(d-K, d+1)); t_lo=min(M[j] for j in js if j<d)-1; t_hi=M[d-1]+1
    def has(j):
        try: return M[j] is not None
        except (KeyError, IndexError): return False
    beta={j:(has(j-1) and M[j]>M[j-1]) for j in js if j<d}   # unknown history below d-K: no defect claimed there
    val={}
    for j in js:
        for t in range(t_lo, t_hi+1):
            if j<d and t>=M[j]+1: val[(j,t)]=r(j,t)
            elif j<d and t==M[j] and beta[j]: val[(j,t)]=1-r(j,t)
    cons=[(j,t) for j in js[2:] for t in range(t_lo, t_hi)]     # y=D_j(t+1), x=D_{j-2}(t), a=D_{j-1}(t), b=D_j(t)
    changed=True
    while changed:
        changed=False
        for (j,t) in cons:
            Y,Xk,A,B=(j,t+1),(j-2,t),(j-1,t),(j,t)
            y,x,a,b=val.get(Y),val.get(Xk),val.get(A),val.get(B)
            new={}
            o = 1 if (a==1 or b==1) else (0 if (a==0 and b==0) else None)   # a or b
            if o is not None:
                if y is None and x is not None: new[Y]=x^o
                if x is None and y is not None: new[Xk]=y^o
            if y is not None and x is not None:
                oo=y^x                       # forced value of (a or b)
                if oo==0:
                    if a is None: new[A]=0
                    if b is None: new[B]=0
                else:
                    if a==0 and b is None: new[B]=1
                    if b==0 and a is None: new[A]=1
            for k,v in new.items():
                if k in val and val[k]!=v: raise ValueError(f"inconsistent at {k}")
                if k not in val: val[k]=v; changed=True
    return val.get((d,M[d-1]+1))
if __name__=='__main__':
    import numpy as np
    import os
    R=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    z=np.load(os.path.join(R,'paper','data','labels_800k.npz')); f=np.load(os.path.join(R,'results','front_800k_derived.npz'))
    per=[int(x) for x in z['per']]; cw=[int(x) for x in z['colword']]; Mf=[int(x) for x in f['M']]
    rr=lambda j,t: (cw[j]>>(t%per[j]))&1
    for d in (28,399,53207,58286,87866):
        obs='stall' if Mf[d]==Mf[d-1] else 'advance'
        res=[(K, {0:'stall',1:'advance',None:'?'}[lookback(rr, Mf, d, K)]) for K in (2,3,4,6,8,12,16,24,32)]
        print(f"d={d}: observed {obs}; verdict by lookback depth K: {res}")
