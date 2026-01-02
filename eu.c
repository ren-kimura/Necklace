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
#include "cov.h"
#include "out.h"

static void usage(const char *s) {
    fprintf(stderr,
            "Usage:\n"
            "\tgenerate: %s -i [in.fa] -k [k] -o [o]\n"
            "\tverify:   %s -m 1 -i [in.fa] -k [k] -d [d] -f [target.fa]\n"
            "Modes:\n"
            "\t-m MODE\t0:generate(default) 1:verify\n\n"
            "Options:\n"
	        "\t-i FILE\t\tinput FASTA file\n"
	        "\t-k INT\t\tk-mer length (>=2 && <=31)\n"
            "\t-d GRAPH TYPE\t0:unidirected 1:bidirected\n"
	        "\t-o OPTION\t0:flat (-d=0 only) 1:pointer (-d=0 only) 2:bp\n"
            "\t-f TARGET FILE\t target file of verification (***.fa)\n",
            s, s);
    exit(EXIT_FAILURE);
}

static int parse_int(const char *s) {
    char *e;
    errno = 0;
    int v = strtol(s, &e, 10);
    if (errno || *e) usage("./eu");
    return v;
}

int main (int argc, char *argv[]) {
    const char *infile = NULL;
    const char *outfile = NULL;
    int k, m, opt, di, out;
    k = -1; m = 0; di = 0; out = 0;

    while ((opt = getopt(argc, argv, "i:k:m:f:d:o:")) != -1) {
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
            case 'd':
                di = parse_int(optarg);
                if (di != 0 && di != 1) usage(argv[0]);
                break;
            case 'o':
                out = parse_int(optarg);
                if (out != 0 && out != 1 && out != 2) usage(argv[0]);
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
        return veri_fa(infile, outfile, k, di);
    } else {
        if (!infile || k == -1) {
            usage(argv[0]);
        }

        printf("infile = %s\n", infile);
        printf("k = %d\n", k);

        Hm *km = NULL; u64 *ka = NULL; u64 N;
        N = extract(infile, k, &km, &ka, 0);
        Rep r; init_rep(&r);
        char *b = rm_ext(infile);

        if (di == 1) {
            VV cc, pp; init_vv(&cc); init_vv(&pp);
            VVb ccb, ppb; init_vvb(&ccb); init_vvb(&ppb);

            // extract canonical k-mers from fasta input
            // construct closed and open paths with strands at the same time
            N = bdextract(infile, k, &km, &ka, &cc, &pp, &ccb, &ppb);

            if (out == 2) {
                r = bbp(km, ka, &cc, &pp, &ccb, &ppb, k);
                u64 np = pp.size;
                wrt(b, &r, k, di, 9, out, np); // cov = 9 (Eulertigs) for convenience
            } else {
                fprintf(stderr, "Error: out=%d is not supported for di=1 in eu\n", out);
            }
            free_vv(&cc); free_vv(&pp); free_vvb(&ccb); free_vvb(&ppb);
        } else {
            Node *g = build(km, k, N);
            balance(&g);

            VV tt; init_vv(&tt);
            etigs(&g, &tt, k);

            if (out == 0) { // yield FASTA
                size_t ns = 0;
                char **ss = spell(&tt, k, &ns);
                char *b = rm_ext(infile);
                if (ss) {
                    wrt_fa(b, k, ss, ns);
                    free_ss(ss, ns);
                }
            } else { // yield .str (possibly plus .ptr)
                VV cc, pp;
                init_vv(&cc); init_vv(&pp);

                tt_to_cc_and_pp(&tt, km, &cc, &pp);

                Rep r; init_rep(&r);
                if (out == 1) {
                    r = ptr(km, ka, &cc, &pp, k);
                } else {
                    r = rbp(km, ka, &cc, &pp, k);
                }

                u64 np = pp.size;
                char *b = rm_ext(infile);
                wrt(b, &r, k, 0, 9, out, np);

                free_vv(&cc); free_vv(&pp);
            }
            free_vv(&tt); free_g(&g);
        }

        free(ka); free_hm(&km); free(b); free_rep(&r);
        return 0;
    }
}
