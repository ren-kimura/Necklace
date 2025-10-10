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

static void usage(const char *s) {
    fprintf(stderr,
	        "Usage: %s -i <file> -k <2-31> -d <0|1> -c <0|1> -o <0|1|2>\n"
	        "\t-i FILE\t\tinput FASTA file\n"
	        "\t-k INT\t\tk-mer length (>=2 && <=31)\n"
            "\t-d GRAPH TYPE\t0:unidirected 1:bidirected\n"
            "\t-c COVER TYPE\t0:matching 1:linearscan 2:greedydfs"
	        "\t-o OPTION\t0:flat 1:pointer 2:bp\n",
	        s);
    exit(EXIT_FAILURE);
}

static int parse_int(const char *s) {
    char *end;
    errno = 0;
    int v = strtol(s, &end, 10);
    if (errno || *end) usage("./necklace");
    return v;
}

int main(int argc, char *argv[]) {
    const char *infile = NULL;
    int k, di, cov, out, opt;
    k = di = cov = out = -1;

	while ((opt = getopt(argc, argv, "i:k:d:c:o:")) != -1) {
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
                if (cov != 0 && cov != 1) usage(argv[0]);
                break;
            case 'o':
                if (out != -1) usage(argv[0]);
                out = parse_int(optarg);
                if (out != 0 && out != 1 && out != 2)
                    usage(argv[0]);
                break;
            default:
                usage(argv[0]);
		}
	}

    if (!infile || k == -1 || di == -1 || cov == -1 || out == -1) {
        usage(argv[0]);
    }

    printf("infile = %s\n", infile);
    printf("k = %d\n", k);
    printf("di = %s\n", (di == 0) ? "uni" : "bi");
    printf("cov = %s\n", (cov == 0) ? "matching" : ((cov == 1) ? "linearscan" : "greedydfs"));
    printf("out = %s\n", (out == 0) ? "flat" : ((out == 1) ? "pointer" : "bp"));

    Hm *km = NULL; u64 *ka = NULL; u64 N;
    VV cc, pp; init_vv(&cc); init_vv(&pp);
    
    if (di == 0) {
        if          (cov == 0) { // maximum matching
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
                return -1;
            }
            decompose(mu, mv, &cc, &pp, N);
            free(mu); free(mv);
        } else if   (cov == 1) { // directly find cover from infile
            fprintf(stderr, "under construction\n");
            exit(EXIT_FAILURE);
        } else if   (cov == 2) { // greedy dfs cover
            fprintf(stderr, "under construction\n");
            exit(EXIT_FAILURE);
        } else {
            fprintf(stderr, "Error: invalid cover type\n");
            exit(EXIT_FAILURE);
        }
        
        if          (out == 0) { // flat
            Rep r = flat(ka, &cc, &pp, k);
            printf("%s\n", r.str);
        } else if   (out == 1) { // pointer
            Rep r = ptr(km, ka, &cc, &pp, k);
            printf("%s\n", r.str);
        } else if   (out == 2) { // bp
            fprintf(stderr, "Error: under construction\n");
            // Rep r = bp(km, ka, &cc, &pp, k);
        } else {
            fprintf(stderr, "Error: invalid out arg\n");
            exit(EXIT_FAILURE);
        }
    } else {
        fprintf(stderr, "under construction\n");
        exit(EXIT_FAILURE);
    }

    free(ka); free_hm(&km);
    free_vv(&cc); free_vv(&pp);    
	return 0;
}