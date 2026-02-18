#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "write.h"

char* rm_ext(const char* f) {
    if (f == NULL) {
        return NULL;
    }
    const char *dot = strrchr(f, '.');
    if (dot == NULL || dot == f) {
        char *b = (char*)malloc(strlen(f) + 1);
        if (b == NULL) return NULL;
        strcpy(b, f);
        return b;
    }
    size_t l = dot - f;
    char *b = (char*)malloc(l + 1);
    if (b == NULL) return NULL;
    memcpy(b, f, l);
    b[l] = '\0';
    return b;
}

void wrt(const char* outfile, const char* r) {
    if (r == NULL) return;
    FILE *fp = fopen(outfile, "w");
    if (fp == NULL) { fprintf(stderr, "Error: cannot open %s for writing\n", outfile); return; }

    fprintf(fp, "%s\n", r);
    fclose(fp);
    printf("%s written\n", outfile);
}