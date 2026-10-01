/* wall_recursion.c: the wall (backward orbit of pi, Lemma 18.2 of the paper) continued by the label
   recursion alone, without any Rule 30 simulation.
   A label r_d is stored as a 64-bit word holding its row-indexed pattern repeated to 64 bits
   (bit t = r_d(t mod 64)); periods are powers of two <= 64.  Given (a, b) = (r_{d-1}, r_d), the next
   label c = r_{d+1} solves c(t+1) = a(t) + b(t) + (1 + b(t)) c(t).  For b != 0 the periodic solution is
   the fixed point of c -> rotl(a ^ b ^ (~b & c), 1), reached from any start within 64 iterations (the
   dependence on the start dies at every 1 of b).  For b == 0 (a frozen diagonal) the two solutions are
   X and ~X with X(t) = xor_{s<t} a(s), of period p(a) or 2 p(a) (Rowland's parity law); the branch is
   taken from BRANCHBITS ('0': X, '1': ~X), and the run stops at a frozen diagonal with no bit left.
   Usage: wall_recursion A_HEX B_HEX D0 DMAX [BRANCHBITS] [PERFILE]
     a = r_{D0-1}, b = r_{D0}; computes r_{d+1} for d = D0 .. DMAX-1; PERFILE receives log2(p_{d+1}) as
     one byte per diagonal (for validation against a frame).  Prints running-max increases of the period,
     frozen diagonals, checkpoints every 10^8 diagonals, the period histogram and the final labels. */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static inline uint64_t rotl(uint64_t x, int k) { return k ? (x << k) | (x >> (64 - k)) : x; }
static int period(uint64_t x) { int p = 64; while (p > 1 && rotl(x, p / 2) == x) p /= 2; return p; }
static int lg(int p) { int k = 0; while ((1 << k) < p) k++; return k; }

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: wall_recursion A_HEX B_HEX D0 DMAX [BRANCHBITS] [PERFILE]\n"); return 2; }
    uint64_t a = strtoull(argv[1], 0, 16), b = strtoull(argv[2], 0, 16);
    long long D0 = atoll(argv[3]), DMAX = atoll(argv[4]);
    const char *bits = argc > 5 ? argv[5] : ""; int nb = 0;
    FILE *pf = argc > 6 ? fopen(argv[6], "wb") : NULL;
    long long hist[8] = {0};
    int pmax = period(a) > period(b) ? period(a) : period(b);
    printf("start d=%lld: a=%016llx (p=%d) b=%016llx (p=%d)\n", D0, (unsigned long long)a, period(a), (unsigned long long)b, period(b));
    for (long long d = D0; d < DMAX; d++) {
        uint64_t c;
        if (b == 0) {
            int pa = period(a);
            uint64_t X = 0; int x = 0;
            for (int t = 0; t < 64; t++) { if (x) X |= 1ULL << t; x ^= (int)((a >> t) & 1); }
            if (x) { printf("d=%lld: frozen diagonal with p_(d-1)=64 of odd weight: period 128 needed, stop\n", d); break; }
            printf("d=%lld: frozen diagonal (r_d = 0), p_(d-1)=%d, weight %d, candidates %016llx / %016llx (period %d)\n",
                   d, pa, __builtin_popcountll(a) / (64 / pa), (unsigned long long)X, (unsigned long long)~X, period(X));
            if (nb >= (int)strlen(bits)) { printf("  no branch bit given: stop\n"); break; }
            c = bits[nb] == '1' ? ~X : X; printf("  branch %c taken\n", bits[nb]); nb++;
        } else {
            c = 0;
            for (int k = 0; k < 64; k++) c = rotl(a ^ b ^ (~b & c), 1);
            if (rotl(a ^ b ^ (~b & c), 1) != c) { printf("d=%lld: no fixed point\n", d); return 1; }
        }
        int pc = period(c); hist[lg(pc)]++;
        if (pf) fputc(lg(pc), pf);
        if (pc > pmax) { pmax = pc; printf("d=%lld: running max of the period reaches %d\n", d + 1, pc); fflush(stdout); }
        if ((d + 1) % 100000000 == 0) { printf("checkpoint d=%lld: a=%016llx b=%016llx\n", d + 1, (unsigned long long)b, (unsigned long long)c); fflush(stdout); }
        a = b; b = c;
    }
    printf("period histogram (p: count):"); for (int k = 0; k < 8; k++) if (hist[k]) printf(" %d: %lld", 1 << k, hist[k]); printf("\n");
    printf("final labels: a=%016llx b=%016llx\n", (unsigned long long)a, (unsigned long long)b);
    if (pf) fclose(pf);
    return 0;
}
