#!/usr/bin/env python3
"""Census of the lookback verdict over every front history consistent with the landing rule (10.237d).

For a frozen diagonal d and depth K: every starting row M_{d-K} (all phases of the period) and every pattern of
moves and stalls on d-K+1..d-1; a move at j lands on min{t >= M_{j-1}+1 : r_{j-1}(t) = 1} (Theorem B(c)).  For
each history the lookback says whether D_d(M_{d-1}+1) is forced (1: the frozen diagonal advances, 0: it stalls).
Usage: python stall_census.py LABELFILE d K    with LABELFILE holding lines "j hexlabel" for j = d-K-1..d-1
   or: python stall_census.py --seed d K       (labels from the 800,000 frame)
"""
import sys, itertools
from stall_lookback import lookback
def census(rr, d, K, P):
    out={'advance':0,'stall':0,'?':0}; examples={'stall':[], '?':[]}
    for p in range(P):
        for bet in itertools.product((0,1), repeat=K-1):
            M={d-K:p}
            for i,j in enumerate(range(d-K+1,d)):
                if bet[i]==0: M[j]=M[j-1]
                else:
                    t=M[j-1]+1
                    while rr(j-1,t)==0: t+=1
                    M[j]=t
            v=lookback(rr,M,d,K); key={1:'advance',0:'stall',None:'?'}[v]; out[key]+=1
            if key!='advance' and len(examples[key])<6: examples[key].append((p,bet,[M[j]%P for j in range(d-K,d)]))
    return out, examples
if __name__=='__main__':
    if sys.argv[1]=='--seed':
        import numpy as np
        import os
        R=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        z=np.load(os.path.join(R,'paper','data','labels_800k.npz')); per=[int(x) for x in z['per']]; cw=[int(x) for x in z['colword']]
        rr=lambda j,t: (cw[j]>>(t%per[j]))&1; d=int(sys.argv[2]); K=int(sys.argv[3]); P=per[d-1]
    else:
        lab={}
        for line in open(sys.argv[1]):
            j,h=line.split(); lab[int(j)]=int(h,16)
        rr=lambda j,t: (lab[j]>>(t%128))&1; d=int(sys.argv[2]); K=int(sys.argv[3]); P=32
    out,ex=census(rr,d,K,P)
    tot=sum(out.values())
    print(f"d={d} K={K}: {tot} histories: forced advance {out['advance']}, forced stall {out['stall']}, undecided {out['?']}")
    for key in ('stall','?'):
        for p,bet,phases in ex[key]: print(f"   {key}: start phase {p}, moves {bet}, fronts mod P {phases}")
