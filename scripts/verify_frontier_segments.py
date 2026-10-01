#!/usr/bin/env python3
"""verify_frontier_segments.py: recompute a wall_recursion128 run step by step, with a second program, in parallel.

wall_recursion128 prints a checkpoint line every 10^8 diagonals. Every step between two consecutive checkpoint
lines of its log (and from the start line to the first) is recomputed here, starting from the log's own labels
at the first line, by the separately written 64-bit program wall_recursion.c; it must end at the log's labels at
the second. If every step agrees, the whole run is confirmed by a second implementation and a second computation.
The 64-bit program holds periods up to 64: a step whose words have period 128 (their halves differ) is recomputed
by wall_recursion128.c instead, and --program 128 uses that program for every step. A step also fails if the
recomputation stops early (a frozen diagonal: the log's run, given no branch bits, would have stopped there too).

Usage, from the repository root (needs python3 and gcc, nothing else):
  python3 scripts/verify_frontier_segments.py LOG [--out FILE] [--jobs N] [--program 64|128]
    LOG        the run's log, e.g. ~/rule30run/wall128_to_1e12.log
    --out      one line per step (default LOG.steps.txt). A rerun skips the steps already there, so the job can be
               stopped (Ctrl-C) and started again; delete a line to have that step run again. When every step is
               done the file is rewritten in order, with a summary at the top.
    --jobs     processes at once (default: the number of CPUs)
One step takes about 7 s with the 64-bit program and 3.5 s with the 128-bit one on one core of a 2.8 GHz Xeon, so
the 8,492 steps of the run to 10^12 need about 17 core-hours (64-bit) or 8 (128-bit).
Exit code 0 when every step of the log is done and agrees and the steps' period histograms add up to the log's.
"""
import argparse, concurrent.futures as cf, os, re, shutil, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
LINE = re.compile(r'^(\d+) (\d+) (agree|DIFFER) ')
HIST = re.compile(r'^period histogram \(p: count\):(.*)$', re.M)
FINAL = re.compile(r'^final labels: a=([0-9a-f]+) b=([0-9a-f]+)', re.M)


def read_log(log):
    s = open(log).read()
    pts = {int(d): (a, b) for d, a, b in re.findall(r'^checkpoint d=(\d+): a=([0-9a-f]+) b=([0-9a-f]+)', s, re.M)}
    st = re.search(r'^start d=(\d+): a=([0-9a-f]+) \(p=\d+\) b=([0-9a-f]+)', s, re.M)
    if st: pts.setdefault(int(st.group(1)), (st.group(2), st.group(3)))
    h = HIST.search(s)
    hist = {int(k): int(v) for k, v in re.findall(r'(\d+): (\d+)', h.group(1))} if h else None
    return sorted(pts.items()), hist, re.findall(r'^d=\d+: .*$', s, re.M)


def run_step(exe64, exe128, prog, d0, lab0, d1, lab1):
    words = (*lab0, *lab1)
    use64 = prog == 64 and all(len(w) == 16 or (len(w) == 32 and w[:16] == w[16:]) for w in words)
    if use64:
        cmd = [exe64, lab0[0][:16], lab0[1][:16], str(d0), str(d1), '']
        want = (lab1[0][:16], lab1[1][:16])
    else:                                    # wall_recursion128 duplicates a 16-digit word to 128 bits itself
        cmd = [exe128, lab0[0], lab0[1], str(d0), str(d1), '', '128', '-']
        want = tuple(w if len(w) == 32 else w + w for w in lab1)
    t = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode < 0:                     # killed (Ctrl-C): not a result, the step stays to do
        raise RuntimeError('interrupted')
    o = p.stdout
    fin, hs = FINAL.search(o), HIST.search(o)
    hist = {int(k): int(v) for k, v in re.findall(r'(\d+): (\d+)', hs.group(1))} if hs else {}
    got = (fin.group(1), fin.group(2)) if fin else None
    ok = got == want and sum(hist.values()) == d1 - d0
    line = (f"{d0} {d1} {'agree' if ok else 'DIFFER'} {64 if use64 else 128} "
            f"{','.join(f'{k}:{v}' for k, v in sorted(hist.items())) or '-'} {time.time() - t:.1f}s")
    if not ok:
        line += f" | got {got} want {want} events {[l for l in o.splitlines() if l.startswith('d=')]}"
    return line


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('log')
    ap.add_argument('--out')
    ap.add_argument('--jobs', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--program', type=int, choices=(64, 128), default=64)
    a = ap.parse_args()
    out = a.out or a.log + '.steps.txt'
    pts, hist, events = read_log(a.log)
    if len(pts) < 2: sys.exit(f'{a.log}: fewer than two checkpoint lines')
    steps = [(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1]) for i in range(len(pts) - 1)]
    done = {}
    if os.path.exists(out):
        for l in open(out):
            m = LINE.match(l)
            if m: done[(int(m.group(1)), int(m.group(2)))] = l.rstrip('\n')
    todo = [s for s in steps if (s[0], s[2]) not in done]
    cc = shutil.which('gcc') or shutil.which('cc')
    if not cc:
        sys.exit('no C compiler found: on Windows, run this in the Ubuntu (WSL) window, not in PowerShell'
                 if os.name == 'nt' else 'no C compiler found: install it with  sudo apt install gcc')
    tmp = tempfile.mkdtemp(prefix='verify_frontier_')
    exe64, exe128 = os.path.join(tmp, 'wall_recursion'), os.path.join(tmp, 'wall_recursion128')
    for exe, src in ((exe64, 'wall_recursion.c'), (exe128, 'wall_recursion128.c')):
        subprocess.run([cc, '-O2', '-o', exe, os.path.join(HERE, src)], check=True)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    print(f'{a.log}: {len(steps)} steps from d={steps[0][0]} to d={steps[-1][2]}; {len(done)} already in {out}; '
          f'{len(todo)} to run, {a.jobs} at a time, {a.program}-bit program', flush=True)
    uneven = [(s[0], s[2]) for s in steps if s[2] - s[0] != 10**8]
    if uneven: print(f'note: {len(uneven)} steps are not 10^8 long (missing checkpoint lines?), e.g. {uneven[:3]}')
    if events: print(f'note: the log reports {len(events)} events, e.g. {events[:3]}')
    t0 = time.time(); n = 0
    ex = cf.ThreadPoolExecutor(a.jobs)
    futs = [ex.submit(run_step, exe64, exe128, a.program, *s) for s in todo]
    try:
        with open(out, 'a') as fh:
            for f in cf.as_completed(futs):
                try: line = f.result()
                except RuntimeError: continue
                fh.write(line + '\n'); fh.flush(); n += 1
                m = LINE.match(line); done[(int(m.group(1)), int(m.group(2)))] = line
                if m.group(3) == 'DIFFER': print(line, flush=True)
                if n % 25 == 0 or n == len(todo):
                    el = time.time() - t0
                    print(f'  {n}/{len(todo)} steps, {el / 60:.0f} min so far, about {el / n * (len(todo) - n) / 60:.0f} min to go',
                          flush=True)
    except KeyboardInterrupt:
        for f in futs: f.cancel()
        ex.shutdown(wait=False)
        print(f'\nstopped after {n} steps; they are in {out}. Run the same command again to continue.')
        sys.exit(130)
    ex.shutdown()
    shutil.rmtree(tmp, ignore_errors=True)
    res = [done.get((s[0], s[2])) for s in steps]
    missing = sum(r is None for r in res)
    if missing:
        print(f'{missing} of {len(steps)} steps not done (interrupted?); run the same command again'); sys.exit(1)
    bad = [r for r in res if LINE.match(r).group(3) == 'DIFFER']
    tot = {}
    for r in res:
        for k, v in re.findall(r'(\d+):(\d+)', r.split()[4]): tot[int(k)] = tot.get(int(k), 0) + int(v)
    fmt = lambda h: ', '.join(f'{k}: {v}' for k, v in sorted(h.items()))
    same = hist is None or tot == hist
    used = sorted({r.split()[3] for r in res})
    head = [f'# {len(steps)} steps of {os.path.basename(a.log)}, d = {steps[0][0]} to {steps[-1][2]}, each recomputed from '
            f"the log's labels at its first line by scripts/verify_frontier_segments.py ({' and '.join(used)}-bit "
            f"program), {time.strftime('%Y-%m-%d')}",
            f'# agree: {len(steps) - len(bad)}   differ: {len(bad)}',
            f'# period histogram of the recomputation: {fmt(tot)}; of the log: '
            f"{fmt(hist) if hist else 'none'}{'' if hist is None else ' (equal)' if same else ' (DIFFERENT)'}",
            '# columns: first diagonal, last diagonal, verdict, program (bits), the step\'s period histogram, seconds']
    with open(out, 'w') as fh: fh.write('\n'.join(head + res) + '\n')
    print('\n'.join(head))
    for r in bad: print(r)
    sys.exit(0 if not bad and same else 1)


if __name__ == '__main__':
    main()
