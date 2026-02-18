#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "ds.h"
#include "utils.h"
#include "veri.h"

static void usage(const char *s) {
    fprintf(stderr, "Usage: %s -i [original.fa] -v [target.fa] -k [k] (-u)\n", s);
    fprintf(stderr, "\t-u: distinguish a k-mer and its reverse complement\n");
    exit(EXIT_FAILURE);
}

int main (int argc, char *argv[]) {
    const char *infile = NULL;
    const char *tgfile = NULL;
    int k = -1;
    bool u_flg = false;
    int opt;

    while ((opt = getopt(argc, argv, "i:v:k:u")) != -1) {
        switch(opt) {
            case 'i': infile = optarg; break;
            case 'v': tgfile = optarg; break;
            case 'k': k = atoi(optarg); break;
            case 'u': u_flg = true; break;
            default: usage(argv[0]);
        }
    }

    if (!infile || !tgfile || k <= 0) usage(argv[0]);

    return veri_fa(infile, tgfile, k, u_flg);
}