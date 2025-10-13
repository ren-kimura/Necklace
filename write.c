#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "write.h"

static void wrt_vle(FILE *fp, u64 val) {
    unsigned char buf[10];
    int i = 0;
    if (val == 0) { // handle the case of zero
        fputc(0, fp);
        return;
    }
    while (val > 0) {
        buf[i] = val & 0x7F; // get lower 7 bits
        val >>= 7;
        if (val > 0) {
            buf[i] |= 0x80; // set continuation bi
        }
        i++;
    }
    fwrite(buf, 1, i, fp);
}

u64 read_vle(FILE *fp) {
    u64 val = 0, s = 0;
    unsigned char byte;
    while (fread(&byte, 1, 1, fp) == 1) {
        val |= (u64)(byte & 0x7F) << s;
        if ((byte & 0x80) == 0) {
            return val;
        }
        s += 7;
        if (s >= 64) {
            return INF;
        }
    }
    return INF;
}

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

void wrt(const char* f, const Rep* r, int k, int di, int cov, int out, u64 np) {
    char sf[FILENAME_MAX];
    snprintf(sf, sizeof(sf), "%s-%d-%d-%d-%d.str", f, k, di, cov, out);
    FILE *sfp = fopen(sf, "w");
    if (sfp == NULL) {
        fprintf(stderr, "Error: cannot open string file for writing\n");
        return;
    }
    fprintf(sfp, "%s\n", r->str);
    fclose(sfp);
    printf("%s written\n", sf);

    if (r->arr != NULL && np > 0) {
        char af[FILENAME_MAX];
        snprintf(af, sizeof(af), "%s-%d-%d-%d-%d.arr", f, k, di, cov, out);
        FILE *afp = fopen(af, "w");
        if (afp == NULL) {
            fprintf(stderr, "Error: cannot open array file for writing\n");
            return;
        }
        wrt_vle(afp, np);
        for (u64 i = 0; i < np; i++) {
            wrt_vle(afp, r->arr[i]);
        }
        fclose(afp);
        printf("%s written\n", af);
    }
}

void vread(const char* f) {
    FILE *fp = fopen(f, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open array file\n");
        return;
    }
    u64 c = read_vle(fp);
    if (c == INF) {
        fprintf(stderr, "Error or file is empty\n");
        fclose(fp);
        return;
    }
    printf("%ld nums in arr\n", c);
    for (u64 i = 0; i < c; i++) {
        u64 v = read_vle(fp);
        if (v == INF) {
            fprintf(stderr, "Error: reached end of file prematurely\n");
            break;
        }
        printf("%ld ", v);
    }
    printf("\n");
    fclose(fp);
}