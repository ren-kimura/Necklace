#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "ds.h"
#include "utils.h"
#include "write.h"
#include "veri.h"
#include "eutils.h"

static void usage(const char *s) {
    fprintf(stderr,
            "Usage:\n"
            "\tgenerate: %s -i [in.fa] -k [k]\n"
            "\tverify:   %s -m 1 -i [in.fa] -k [k] -f [target.fa]\n",
            s, s);
    exit(EXIT_FAILURE);
}

static int parse_int(const char *s) {
    char *e;
    errno = 0;
    int v = strtol(s, &e, 10);
    if (errno || *e) usage("./nkl");
    return v;
}

int main (int argc, char *argv[]) {
    const char *infile = NULL;
    const char *outfile = NULL;
    int k, m, opt;
    k = -1; m = 0;

    while ((opt = getopt(argc, argv, "i:k:m:f:")) != -1) {
        switch (opt) {
            case 'i':
                if (infile) usage(argv[0]);
                infile = optarg;
                break;
            case 'k':
                if (k != -1) usage(argv[0]);
                k = parse_int(optarg);
                if (k < 2 || k > 31) usage(argv[0]);
                break;
            case 'm':
                m = parse_int(optarg);
                if (m != 0 && m != 1) usage(argv[0]);
                break;
            case 'f':
                outfile = optarg;
                break;
            default:
                usage(argv[0]);
        }
    }
    if (m == 1) {
        if (!infile || !outfile || k == -1) {
            fprintf(stderr, "Error: provide -i -k\n");
            usage(argv[0]);
        }
        return veri_fa(infile, outfile, k);
    } else {
        if (!infile || k == -1) {
            usage(argv[0]);
        }

        printf("infile = %s\n", infile);
        printf("k = %d\n", k);

        Hm *km = NULL; u64 *ka = NULL; u64 N;
        N = extract(infile, k, &km, &ka, 0);

        Node *g = build(km, k, N);
        balance(&g);

        VV tt; init_vv(&tt);
        etigs(&g, &tt, k);

        size_t ns = 0;
        char **ss = spell(&tt, k, &ns);
        char *b = rm_ext(infile);
        if (ss) {
            wrt_fa(b, k, ss, ns);
            free_ss(ss, ns);
        }
        free(b);
        free(ka); free_hm(&km); free_g(&g); free_vv(&tt);
        return 0;
    }
}