#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "ds.h"
#include "utils.h"
#include "eutils.h"
#include "cov.h"
#include "out.h"
#include "stat.h"
#include "write.h"
#include "veri.h"

static void usage(const char *s) {
    fprintf(stderr,
            "Usage:\n"
	        "\tgenerate: %s -i [in.fa] -k [k] -a [a] (-p) (-u)\n"
            "\tverify:   %s -i [in.fa] -k [k] (-p) (-u) -v [target.str]\n"
            "\tsparsity: %s -i [in.fa] -k [k] (-u) -s\n\n"
	        "\t-i: set a path to an input FASTA file\n"
            "\t-o: specify output str filename (optional, default:in_k.str)\n"
	        "\t-k: choose k-mer length k s.t. 1 < k < 32\n"
            "\t-a: specify an algorithm to be run (eu:Eulertigs(default) fg:FullGreedy gb:GreedyBaseline ba:BaselineA gc:GreedyCover)\n"
            "\t-p: parenthesis representation (optional except -a fg and gb)\n"
            "\t-u: distinguish a k-mer and its reverse complement (optional)\n\n"
            "\t-v: verify the specified output\n"
            "\t-s: measure the dBG sparsity as average outdegree\n"
            "\t-h: print help\n",
	        s, s, s);
    exit(EXIT_FAILURE);
}

static int parse_int(const char *s) {
    char *e;
    errno = 0;
    int v = strtol(s, &e, 10);
    if (errno || *e) usage("./nkl");
    return v;
}

int main(int argc, char *argv[]) {
    const char *infile = NULL;
    const char *target_file = NULL;
    char *outfile = NULL;
    const char *algo = "eu"; // default
    int k = -1;
    bool p_flg = false, u_flg = false, s_flg = false;
    int opt;

	while ((opt = getopt(argc, argv, "i:o:k:a:puv:sh")) != -1) {
	    switch (opt) {
		    case 'i': infile = optarg; break;
            case 'k':
                k = parse_int(optarg);
                if (k < 2 || k > 31) usage(argv[0]);
                break;
            case 'a': algo = optarg; break;
            case 'p': p_flg = true; break;
            case 'u': u_flg = true; break;
            case 'v': target_file = optarg; break;
            case 's': s_flg = true; break;
            case 'o': outfile = strdup(optarg); break;
            case 'h': usage(argv[0]); break;
            default: usage(argv[0]);
		}
	}

    if (!infile || k == -1) usage(argv[0]);

    if (outfile == NULL) {
        // infile = "path/to/c.fa", k = 11 -> "c_11.str"
        char *base = rm_ext(infile);
        size_t out_len = strlen(base) + 12;
        outfile = (char *)malloc(out_len);
        if (!outfile) { perror("malloc"); exit(EXIT_FAILURE); }

        sprintf(outfile, "%s_%d.str", base, k);
        free(base);
    }

    if (strcmp(algo, "fg") == 0 || strcmp(algo, "gb") == 0) {
        if (p_flg == false) {
            fprintf(stderr, "Warning: this algorithm requires p_flg to be true. Overwritten\n");
            p_flg = true;
        }
    }

    if (target_file) { // verify
        if (!infile || !target_file ||  k == -1) usage(argv[0]);
        return veri(infile, target_file, k, u_flg);
    } else if (s_flg) { // graph stat
        if (!infile || k == -1) usage(argv[0]);
        printf("sparsity analysis\n");
        printf("infile = %s, k = %d, %s\n", infile, k, u_flg ? "uni-directional" : "bi-directional");

        Hs *ks = NULL;
        u64 N = extract_hs(infile, k, &ks, u_flg);
        if (N == 0) {
            fprintf(stderr, "Error: no k-mers found\n");
            return -1;
        }

        u64 M = 0; Hs *cur, *tmp; u64 z = 0;
        
        HASH_ITER (hh, ks, cur, tmp) {
            prog(++z, N, "measuring sparsity");
            u64 h = cur->key;

            if (u_flg) {
                for (int i = 0; i < 4; i++) {
                    if (step_hs(ks, k, h, B[i], true) != INF) {
                        M++;
                    }
                }
            } else {
                for (int from = 0; from < 2; from++) {
                    for (int i = 0; i < 4; i++) {
                        for (int to = 0; to < 2; to++) {
                            if (bstep_hs(ks, k, h, B[i], true, from, to) != INF) {
                                M++;
                            }
                        }
                    }
                }
            }
        }
        fin("measuring sparsity");

        printf("|V| (unique k-mers): %lu\n", N);
        printf("|E| (total edges)  : %lu\n", M);
        if (N > 0) {
            printf("Sparsity (|E|/|V|) : %.4f\n", (double)M/N);
        }

        free_hs(&ks);
    } else {
        if (!infile || k == -1) usage(argv[0]);

        printf("infile = %s\n", infile);
        printf("k = %d\n", k);
        printf("algorithm = %s\n", algo);
        printf("%s\n", u_flg ? "uni-directional" : "bi-directional");
        if (p_flg) printf("parenthesis represenation\n");

        Hm *km = NULL; u64 *ka = NULL; u64 N;
        char* r = NULL;

        if (u_flg) {
            if (strcmp(algo, "fg") == 0 || strcmp(algo, "gb") == 0) {
                N = extract(infile, k, &km, &ka, u_flg);
                if (strcmp(algo, "fg") == 0) {
                    r = full_greedy(km, ka, k);
                } else {
                    r = greedy_baseline_close(km, ka, k);
                }
            } else {
                VV cc, pp; init_vv(&cc); init_vv(&pp);
                if (strcmp(algo, "ba") == 0) {
                    N = baseline_a(infile, k, &km, &ka, &cc, &pp);
                } else {
                    N = extract(infile, k, &km, &ka, u_flg);
                    if (strcmp(algo, "eu") == 0) {
                        Node *g = build(km, k, N);
                        balance(&g);
                        VV tt; init_vv(&tt);
                        etigs(&g, &tt, k);
                        tt_to_cc_and_pp(&tt, km, &cc, &pp);
                        free_vv(&tt);
                    } else if (strcmp(algo, "gc") == 0) {
                        greedy_cover(km, ka, &cc, &pp, k);
                    }
                }
                if (p_flg) {
                    r = necklace_cover2(km, ka, &cc, &pp, k);
                } else {
                    r = pccover_to_cspss(ka, &cc, &pp, k);
                }
                free_vv(&cc); free_vv(&pp);
            }
        } else {
            if (strcmp(algo, "gb") == 0) {
                N = extract(infile, k, &km, &ka, u_flg);
                r = bi_greedy_baseline(km, ka, k);
            } else {
                VV cc, pp; init_vv(&cc); init_vv(&pp);
                VVb ccb, ppb; init_vvb(&ccb); init_vvb(&ppb);
                if (strcmp(algo, "ba") == 0) {
                    N = bi_baseline_a(infile, k, &km, &ka, &cc, &pp, &ccb, &ppb);
                } else {
                    N = extract(infile, k, &km, &ka, u_flg);
                    if (strcmp(algo, "eu") == 0) {
                        printf("This option assumes an input be precomputed BiEulertigs. Proceed [y/n]: ");
                        fflush(stdout);
                        int c = getchar();
                        if (c != '\n' && c != EOF) {
                            int tmp;
                            while((tmp = getchar()) != '\n' && tmp != EOF);
                        }
                        if (c == 'y' || c == 'Y') {
                            N = bi_baseline_a(infile, k, &km, &ka, &cc, &pp, &ccb, &ppb);
                        } else {
                            return 0;
                        }
                    } else if (strcmp(algo, "gc") == 0) {
                        bi_greedy_cover(km, ka, &cc, &pp, &ccb, &ppb, k);
                    }
                }
                if (p_flg) {
                    r = bi_necklace_cover(km, ka, &cc, &pp, &ccb, &ppb, k);
                } else {
                    r = bi_pccover_to_cspss(ka, &cc, &pp, &ccb, &ppb, k);
                }
                free_vv(&cc); free_vv(&pp); free_vvb(&ccb); free_vvb(&ppb);
            }
        }
        free(ka); free_hm(&km);

        wrt(outfile, r);
        free(r);
    }

    free((char*)outfile);
    return 0;
}