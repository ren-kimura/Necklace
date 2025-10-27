#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "ds.h"
#include "utils.h"
#include "cov.h"
#include "out.h"
#include "write.h"
#include "veri.h"

static void usage(const char *s) {
    fprintf(stderr,
            "Usage:\n"
	        "\tgenerate: %s -i [in.fa] -k [k] -d [d] -c [c] -o [o]\n"
            "\tverify:   %s -m 1 -i [in.fa] -k [k] -d [d] -o [0] -f [target.str]\n\n"
            "Modes:\n"
            "\t-m MODE\t0:generate(default) 1:verify\n\n"
            "Options:\n"
	        "\t-i FILE\t\tinput FASTA file\n"
	        "\t-k INT\t\tk-mer length (>=2 && <=31)\n"
            "\t-d GRAPH TYPE\t0:unidirected 1:bidirected\n"
            "\t-c COVER TYPE\t0:matching(under d=0) 1:linearscan 2:greedycover 3:greedydfs(under o=2)\n"
	        "\t-o OPTION\t0:flat 1:pointer 2:bp\n"
            "\t-f TARGET FILE\t target file of verification (will replace .str with .arr when o ==1)\n",
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

int main(int argc, char *argv[]) {
    const char *infile = NULL;
    const char *outfile = NULL;
    int k, di, cov, out, opt;
    k = di = cov = out = -1;
    int m = 0;

	while ((opt = getopt(argc, argv, "i:k:d:c:o:m:f:")) != -1) {
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
                if (di != -1) usage(argv[0]);
                di = parse_int(optarg);
                if (di != 0 && di != 1) usage(argv[0]);
                break;
            case 'c':
                if (cov != -1) usage(argv[0]);
                cov = parse_int(optarg);
                if (cov != 0 && cov != 1 && cov != 2 && cov != 3) usage(argv[0]);
                break;
            case 'o':
                if (out != -1) usage(argv[0]);
                out = parse_int(optarg);
                if (out != 0 && out != 1 && out != 2)
                    usage(argv[0]);
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
        if (!infile || !outfile ||  k == -1 || di == -1 || out == -1) {
            fprintf(stderr, "Error: provide -i -k -d -o -f\n");
            usage(argv[0]);
        }
        return veri(infile, outfile, k, di, out);
    } else {
        if (!infile || k == -1 || di == -1 || cov == -1 || out == -1) {
            usage(argv[0]);
        }

        printf("infile = %s\n", infile);
        printf("k = %d\n", k);
        printf("di = %s\n", (di == 0) ? "uni" : "bi");
        printf("cov = %s\n", (cov == 0) ? "matching" : (cov == 1) ? "linearscan" : (cov == 2) ? "greedycover" : "greedydfs");
        printf("out = %s\n", (out == 0) ? "flat" : (out == 1) ? "pointer" : "bp");

        Hm *km = NULL; u64 *ka = NULL; u64 N;
        
        if (di == 0) {
            VV cc, pp; init_vv(&cc); init_vv(&pp);
            if (cov == 0) { // maximum matching
                N = extract(infile, k, &km, &ka, di);        
                u64 *mu = (u64*)malloc(N * sizeof(u64));
                u64 *mv = (u64*)malloc(N * sizeof(u64));
                if (mu == NULL || mv == NULL) {
                    fprintf(stderr, "Error: malloc failed for array mu or mv\n");
                    free(ka); free(mu); free(mv);
                    free_hm(&km);
                    return -1;
                }
                u64 M = mbm(km, ka, mu, mv, k, N);
                if (M == INF) {
                    fprintf(stderr, "Warning: too large matching\n");
                    free(ka); free(mu); free(mv);
                    free_hm(&km);
                    exit(EXIT_FAILURE);
                }
                decompose(mu, mv, &cc, &pp, N);
                free(mu); free(mv);
            } else if (cov == 1) { // directly find cover from infile
                N = dextract(infile, k, &km, &ka, &cc, &pp);
            } else if (cov == 2) { // greedy dfs from unvisited vertices
                N = extract(infile, k, &km, &ka, di);
                gcov(km, ka, &cc, &pp, k);
            } else {
                fprintf(stderr, "Error: invalid cover type\n");
                exit(EXIT_FAILURE);
            }
            
            Rep r; init_rep(&r);
            if (out == 0) { // flat
                r = flat(ka, &cc, &pp, k);
            } else if (out == 1) { // pointer
                r = ptr(km, ka, &cc, &pp, k);
            } else if (out == 2) { // bp
                r = bp(km, ka, &cc, &pp, k);
            } else {
                fprintf(stderr, "Error: invalid out arg\n");
                exit(EXIT_FAILURE);
            }
            u64 np = pp.size;
            free_vv(&cc); free_vv(&pp); 
            char* b = rm_ext(infile);
            wrt(b, &r, k, di, cov, out, np);
            free_rep(&r); free(b);
        } else {
            W w; init_w(&w);
            VV cc, pp; init_vv(&cc); init_vv(&pp);
            VVb ccb, ppb; init_vvb(&ccb); init_vvb(&ppb);
            if (cov == 0) { // greedy dfs
                fprintf(stderr, "under construction\n");
                exit(EXIT_FAILURE);
            } else if (cov == 1) { // directly find cover from infile
                N = bdextract(infile, k, &km, &ka, &w);
                // disp_hm(km, k);
            } else if (cov == 2) { // greedy covering from unvisited vertices
                N = extract(infile, k, &km, &ka, di);
                bgcov(km, ka, &cc, &pp, &ccb, &ppb, k);
            } else if (cov == 3) {
                if (out != 2) {
                    fprintf(stderr, "Error: cov=3 (greedydfs) is only compatible with out=2(bp)\n");
                    fprintf(stderr, "Overwritten out=%d -> 2\n", out);
                    out = 2;
                }
                N = extract(infile, k, &km, &ka, di);
                // bgdfs(km, ka, &cc, &pp, &ccb, &ppb, k);
            } else {
                fprintf(stderr, "Error: invalid cover type\n");
                exit(EXIT_FAILURE);
            }
            Rep r; init_rep(&r);
            if (out == 0) {
                r = bflat(ka, &cc, &pp, &ccb, &ppb, k);
            } else if (out == 1) {
                fprintf(stderr, "under construction\n");
                exit(EXIT_FAILURE);
            } else if (out == 2) {
                if (cov == 3) {
                    r = flat_w(&w);
                } else {
                    fprintf(stderr, "under construction\n");
                    exit(EXIT_FAILURE);
                }
            } else {
                fprintf(stderr, "Error: invalid out arg\n");
                exit(EXIT_FAILURE);
            }
            u64 np = 0;
            if (w.pp) for (np = 0; w.pp[np] != NULL; np++);
            char* b = rm_ext(infile);
            wrt(b, &r, k, di, cov, out, np);
            free_rep(&r);
            free_w(&w);
        }

        free(ka); free_hm(&km);   
        return 0;
    }
}