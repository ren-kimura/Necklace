#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include "ds.h"
#include "utils.h"
#include "mbm.h"
#include "rep_utils.h"

void progress(size_t now, size_t total, const char* task) {
    if (total == 0) return; // avoid division by zero
    if (total > 100) {
        size_t interval = total / 100;
        if (now != 0 && now != total && (now % interval) != 0) {
            return;
        }
    }
    int percentage = (int)((now * 100) / total);
    printf("\r%s: %d%%", task, percentage);
    fflush(stdout);
}

void finished(const char* task) {
    printf("\r%s: 100%%\n", task);
    fflush(stdout);
}

static void usage(const char *prog) {
    fprintf(stderr,
	        "Usage: %s -i <file> -k <int> -o <0|1|2|10>\n"
	        "\t-i FILE\t\tinput FASTA file\n"
	        "\t-k INT\t\tk-mer length (>=1 && <=31)\n"
	        "\t-o OPTION\t0:naive 1:pointer 2:pseudoforest\n"
	        "\t\t\t10:eulertigs(uni-directed)\n",
	        prog);
    exit(EXIT_FAILURE);
}

static int parse_int(const char *s) {
    char *end;
    errno = 0;
    int v = strtol(s, &end, 10);
    if (errno || *end) usage("./rep");
    return v;
}

int main(int argc, char *argv[]) {
    const char *infile = NULL;
    int k, di, cover, mode, opt;
    k = di = cover = mode = -1;

	while ((opt = getopt(argc, argv, "i:k:d:c:o:")) != -1) {
	    switch (opt) {
		    case 'i':
			    if (infile) usage(argv[0]);
                infile = optarg;
                break;
            case 'k':
                if (k != -1) usage(argv[0]);
                k = parse_int(optarg);
                if (k < 1 || k > 31) usage(argv[0]);
                break;
            case 'd':
                if (di != -1) usage(argv[0]);
                di = parse_int(optarg);
                if (di != 0 && di != 1) usage(argv[0]);
                break;
            case 'c':
                if (cover != -1) usage(argv[0]);
                cover = parse_int(optarg);
                if (cover != 0 && cover != 1) usage(argv[0]);
                break;
            case 'o':
                if (mode != -1) usage(argv[0]);
                mode = parse_int(optarg);
                if (mode != 0 && mode != 1 && mode != 2)
                    usage(argv[0]);
                break;
            default:
                usage(argv[0]);
		}
	}

    if (!infile || k == -1 || di == -1 || cover == -1 || mode == -1) {
        usage(argv[0]);
    }

    printf("input = %s\n", infile);
    printf("k = %d\n", k);
    printf("di = %d\n", di);
    printf("cover = %d\n", cover);
    printf("mode = %d\n", mode);

    Hm *km = NULL;
    u64 *ka = NULL;
    u64 N = extract(infile, k, &km, &ka);
    
    u64 *mu = malloc(N * sizeof(u64));
    u64 *mv = malloc(N * sizeof(u64));
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

    VV cc, pp;
    init_vv(&cc); init_vv(&pp);
    decompose(mu, mv, &cc, &pp, N);
    
    if (mode == 0) {
        Rep r = flat(ka, &cc, &pp, k);
        printf("%s\n", r.str);
    } else if (mode == 1) {
        Rep r = ptr(km, ka, &cc, &pp, k);
        printf("Generated text:\n%s\n", r.str);
    } else if (mode == 2) {
        // Rep r = bp(km, ka, &cc, &pp, k);
    } else {
        fprintf(stderr, "Error: invalid mode for output format\n");
        exit(EXIT_FAILURE);
    }

    free(ka); free(mu); free(mv);
    free_hm(&km);
    free_vv(&cc); free_vv(&pp);
    
	return 0;
}