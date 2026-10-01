/* label_distance.c (derived from wall_recursion128.c): the same recursion, recording the Hamming distance
   between consecutive labels over one period, per period regime.  Output: for each period P, the counts of
   distance k = 0..P.  Usage: label_distance A_HEX B_HEX D0 DMAX [BRANCHBITS]
   original header: the label recursion of wall_recursion.c with 128-bit labels, for walls of period
   up to 128.  A label is a 128-bit word holding its row-indexed pattern repeated to 128 bits
   (bit t = r(t mod 128)).  Given (a, b) = (r_{d-1}, r_d), r_{d+1} is the fixed point of
   c -> rotl(a ^ b ^ (~b & c), 1), reached from 0 in P = max(period a, period b) iterations when b != 0
   (the dependence on the start dies at every 1 of b); when b == 0 the two solutions are X and ~X with
   X(t) = xor_{s<t} a(s) (period p(a) or 2 p(a)).  Branches at frozen diagonals come from BRANCHBITS
   ('0': X, '1': ~X); the run stops at a frozen diagonal with no bit left, at a period-256 onset, or when
   the running maximum of the period reaches STOPP, printing the label pair at that point.
   Usage: wall_recursion128 A_HEX B_HEX D0 DMAX [BRANCHBITS] [STOPP] [PERFILE|-] [COLFILE]
     COLFILE receives the frozen universe's centre column r_d(d) for d = D0+1 .. DMAX, packed 8 bits per byte (bit d - D0 - 1)
     A_HEX, B_HEX: 16 hex digits (a 64-bit word as wall_recursion.c prints, duplicated to 128 bits) or
     32 hex digits (a 128-bit word). PERFILE receives log2(p_{d+1}) per diagonal. */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned __int128 u128;
static inline u128 rotl(u128 x, int k) { return k ? (x << k) | (x >> (128 - k)) : x; }
static int period(u128 x) { int p = 128; while (p > 1 && rotl(x, p / 2) == x) p /= 2; return p; }
static int lg(int p) { int k = 0; while ((1 << k) < p) k++; return k; }
static int popcount(u128 x) { return __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x >> 64)); }
static u128 parse_hex(const char *s) {
    u128 v = 0; size_t n = strlen(s);
    for (; *s; s++) { int c = *s | 32; int d = (c >= '0' && c <= '9') ? c - '0' : c - 'a' + 10; v = (v << 4) | (u128)d; }
    if (n <= 16) v |= v << 64;                       /* a 64-bit periodic word, duplicated */
    return v;
}
static void print128(u128 x) { printf("%016llx%016llx", (unsigned long long)(x >> 64), (unsigned long long)x); }

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: wall_recursion128 A_HEX B_HEX D0 DMAX [BRANCHBITS] [STOPP] [PERFILE]\n"); return 2; }
    u128 a = parse_hex(argv[1]), b = parse_hex(argv[2]);
    long long D0 = atoll(argv[3]), DMAX = atoll(argv[4]);
    const char *bits = argc > 5 ? argv[5] : ""; int nb = 0;
    int STOPP = argc > 6 ? atoi(argv[6]) : 0;
    FILE *pf = (argc > 7 && strcmp(argv[7], "-")) ? fopen(argv[7], "wb") : NULL;
    FILE *cf = argc > 8 ? fopen(argv[8], "wb") : NULL;     /* COLFILE: bit d of the column = r_d(d), packed 8 per byte */
    unsigned char cbyte = 0; int cbits = 0;
    long long hist[9] = {0}; static long long dist[8][129];
    int pmax = period(a) > period(b) ? period(a) : period(b);
    printf("start d=%lld: a=", D0); print128(a); printf(" (p=%d) b=", period(a)); print128(b); printf(" (p=%d)\n", period(b));
    for (long long d = D0; d < DMAX; d++) {
        u128 c;
        if (b == 0) {
            int pa = period(a);
            u128 X = 0; int x = 0;
            for (int t = 0; t < 128; t++) { if (x) X |= (u128)1 << t; x ^= (int)((a >> t) & 1); }
            if (x) { printf("d=%lld: frozen diagonal with p_(d-1)=128 of odd weight: period 256 needed, stop\n", d); break; }
            printf("d=%lld: frozen diagonal (r_d = 0), p_(d-1)=%d, weight %d, candidates ", d, pa, popcount(a) / (128 / pa));
            print128(X); printf(" / "); print128(~X); printf(" (period %d)\n", period(X));
            if (nb >= (int)strlen(bits)) { printf("  no branch bit given: stop\n"); break; }
            c = bits[nb] == '1' ? ~X : X; printf("  branch %c taken\n", bits[nb]); nb++; fflush(stdout);
        } else {
            /* the dependence on the start survives k iterations only across a zero run of b of length k,
               so (longest cyclic zero run of b) + 1 iterations reach the fixed point */
            u128 z = ~b; int L = 0;
            while (z) { z &= rotl(z, 1); L++; }
            c = 0;
            for (int k = 0; k <= L; k++) c = rotl(a ^ b ^ (~b & c), 1);
            if (rotl(a ^ b ^ (~b & c), 1) != c) { printf("d=%lld: no fixed point\n", d); return 1; }
        }
        int pc = period(c); hist[lg(pc)]++;
        { int P = period(b) > pc ? period(b) : pc; int k = popcount(b ^ c) / (128 / P); dist[lg(P)][k]++; }
        if (pf) fputc(lg(pc), pf);
        if (cf) { cbyte |= (unsigned char)(((c >> ((d + 1) & 127)) & 1) << cbits); if (++cbits == 8) { fputc(cbyte, cf); cbyte = 0; cbits = 0; } }
        if (pc > pmax) {
            pmax = pc; printf("d=%lld: running max of the period reaches %d\n", d + 1, pc); fflush(stdout);
            if (STOPP && pc >= STOPP) {
                printf("labels for the machine at d=%lld: a=", d + 1); print128(b); printf(" b="); print128(c); printf("\n");
                a = b; b = c; break;
            }
        }
        if ((d + 1) % 100000000 == 0) { printf("checkpoint d=%lld: a=", d + 1); print128(b); printf(" b="); print128(c); printf("\n"); fflush(stdout); }
        a = b; b = c;
    }
    for (int j = 0; j < 8; j++) { long long tot = 0; for (int k = 0; k <= 128; k++) tot += dist[j][k];
        if (tot) { printf("distance histogram at period %d (%lld pairs):", 1 << j, tot); for (int k = 0; k <= (1 << j); k++) printf(" %d:%lld", k, dist[j][k]); printf("\n"); } }
    printf("period histogram (p: count):"); for (int k = 0; k < 9; k++) if (hist[k]) printf(" %d: %lld", 1 << k, hist[k]); printf("\n");
    printf("final labels: a="); print128(a); printf(" b="); print128(b); printf("\n");
    if (pf) fclose(pf);
    if (cf) { if (cbits) fputc(cbyte, cf); fclose(cf); }
    return 0;
}
