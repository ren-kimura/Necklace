#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "veri.h"
#include "utils.h"

static bool proc_rm(const char* s, int k, int di, Hs** ks);
static bool fveri(const char* cont, int k, int di, Hs** ks);

static bool proc_rm(const char* s, int k, int di, Hs** ks) {
    u64 h = enc(s, k);
    if (di) h = can(h, k);

    if (find_hs(*ks, h)) {
        del_hs(ks, h);
        return true;
    } else {
        char t[k + 1];
        dec(h, k, t);
        fprintf(stderr, "Error: target file has an extra k-mer: %s\n", t);
        return false;
    }
}

static bool fveri(const char* ss, int k, int di, Hs** ks) {
    char* tt = strdup(ss);
    if (!tt) return false;

    char* ct = NULL;
    char* pt = NULL;
    
    char* sp = strstr(tt, ",,");
    if (sp != NULL) {
        *sp = '\0';
        ct = tt;
        pt = sp + 2;
    } else {
        pt = (tt[0] == ',') ? tt + 1 : NULL;
        ct = (tt[0] != ',') ? tt : NULL;
    }

    char buf[k + 1];
    buf[k] = '\0';

    // process circular strings
    if (ct && *ct != '\0') {
        char* to = strtok(ct, ",");
        while (to) {
            size_t l = strlen(to);
            if (l > 0) {
                for (size_t i = 0; i < l; ++i) {
                    for (int j = 0; j < k; ++j) buf[j] = to[(i + j) % l];
                    if (!proc_rm(buf, k, di, ks)) {
                        free(tt);
                        return false;
                    }
                }
            }
            to = strtok(NULL, ",");
        }
    }

    // process non-circular strings
    if (pt && *pt != '\0') {
        char* to = strtok(pt, ",");
        while (to) {
            size_t l = strlen(to);
            if (l >= (size_t)k) {
                for (size_t i = 0; i <= l - k; ++i) {
                    strncpy(buf, to + i, k);
                    if (!proc_rm(buf, k, di, ks)) {
                        free(tt);
                        return false;
                    }
                }
            }
            to = strtok(NULL, ",");
        }
    }

    free(tt);
    return true;
}

static bool bveri(const char* ss, int k, int di, Hs** ks) {
    char* tt = strdup(ss);
    if (!tt) return false;

    char buf[k + 1];
    buf[k] = '\0';

    u64 mask = (1ULL << (2 * (k - 1))) - 1;

    char* to = strtok(tt, ",");
    while (to) {
        if (*to == '\0') {
            to = strtok(NULL, ",");
            continue;
        }

        size_t l = strlen(to);
        if (l < (size_t)k) {
            if (l > 0) {
                strncpy(buf, to, k);
                if (!proc_rm(buf, k, di, ks)) {
                    free(tt);
                    return false;
                }
            }
        } else {
            strncpy(buf, to, k);
            if (!proc_rm(buf, k, di, ks)) {
                free(tt);
                return false;
            }
            St s;
            init_st(&s);
            u64 pre = enc(buf, k) & mask;
            for (size_t i = 1; i <= l - k; ++i) {
                char cur = to[i + k - 1];
                if (isalpha(cur)) {
                    u64 h = (pre << 2) | (cur == 'A' ? 0 : cur == 'C' ? 1 : cur == 'G' ? 2 : 3);
                    dec(h, k, buf);
                    if (!proc_rm(buf, k, di, ks)) {
                        free(tt);
                        while(!is_empty_st(&s)) pop(&s);
                        return false;
                    }
                    pre = h & mask;
                } else if (cur == '(') {
                    push(&s, pre);
                } else if (cur == ')') {
                    if (!is_empty_st(&s)) {
                        pre = pop(&s);
                    }
                }
            }
            // スタックをクリーンアップ
            while(!is_empty_st(&s)) pop(&s);
        }
        to = strtok(NULL, ",");
    }

    free(tt);
    return true;
}

int veri(const char* of, const char* tf, int k, int di, int out) {
    fprintf(stdout, "verification mode\n");
    fprintf(stdout, "original file: %s\n", of);
    fprintf(stdout, "target file: %s\n", tf);
    fprintf(stdout, "\n[Step 1/3] extracting k-mers from original file\n");
    Hm *km = NULL;
    u64 *ka = NULL;
    u64 no = extract(of, k, &km, &ka, di);
    Hs *ks = NULL;
    Hm *s, *tmp;
    HASH_ITER(hh, km, s, tmp) { add_hs(&ks, s->key); }
    free(ka); free_hm(&km);
    fprintf(stdout, "%ld k-mers in %s\n", no, of);

    fprintf(stdout, "\n[Step 2/3] reconstructing k-mers from target file\n");
    FILE* fp = fopen(tf, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open target file\n");
        free_hs(&ks);
        return -1;
    }
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *sc = malloc(fsize + 1);
    size_t ir = fread(sc, 1, fsize, fp);
    if (ir != (size_t)fsize) {
        fprintf(stderr, "Error: cannot read whole file '%s', expected %ld bytes but got %ld\n", tf, fsize, ir);
        fclose(fp);
        free(sc);
        return -1;
    }
    fclose(fp);
    sc[fsize] = 0;
    if (fsize > 0 && sc[fsize - 1] == '\n') {
        sc[fsize - 1] = '\0';
    }
    bool z = false;
    switch (out) {
        case 0: z = fveri(sc, k, di, &ks); break;
        case 1:
        case 2: z = bveri(sc, k, di, &ks); break;
        default: fprintf(stderr, "veri for out=%d is not supported\n", out); break;
    }
    free(sc);
    if (!z) { free_hs(&ks); return -1; }

    fprintf(stdout, "\n[Step 3/3] final check\n");
    bool fz = false;
    if (HASH_COUNT(ks) == 0) {
        fprintf(stdout, "All original k-mers were found\n");
        fz = true;
    } else {
        fprintf(stderr, "Error: %d k-mer(s) were missing\n", HASH_COUNT(ks));
        Hs* i;
        int c = 0;
        for (i = ks; i != NULL && c < 5; i = i->hh.next) {
            char s[k + 1];
            dec(i->key, k, s);
            printf("  - missing: %s\n", s);
            c++;
        }
    }
    free_hs(&ks);
    return fz ? 0 : -1;
}

int veri_fa(const char *of, const char *tf, int k) {
    fprintf(stdout, "verification mode\n");
    fprintf(stdout, "original file: %s\n", of);
    fprintf(stdout, "target file: %s\n", tf);

    fprintf(stdout, "\n[Step 1/3] extracting k-mers from original file\n");
    Hm *km = NULL;
    u64 *ka = NULL;
    u64 no1 = extract(of, k, &km, &ka, 0);
    Hs *ks = NULL;
    Hm *s, *tmp;
    HASH_ITER(hh, km, s, tmp) { add_hs(&ks, s->key); }
    free(ka); free_hm(&km);
    fprintf(stdout, "%ld k-mers in %s\n", no1, of); 
    
    fprintf(stdout, "\n[Step 2/3] reconstructing k-mers from target file\n");
    km = NULL;
    ka = NULL;
    u64 no2 = extract(tf, k, &km, &ka, 0);
    fprintf(stdout, "%ld k-mers in %s\n", no2, tf);

    fprintf(stdout, "\n[Step 3/3] final check\n");
    bool fz = true;

    HASH_ITER(hh, km, s, tmp) {
        if (find_hs(ks, s->key)) {
            del_hs(&ks, s->key);
        } else {
            char t[k + 1];
            dec(s->key, k, t);
            fprintf(stderr, "Error: found a k-mer not in original: %s\n", t);
            fz = false;
            break;
        }
    }
    free(ka); free_hm(&km);

    if (fz && HASH_COUNT(ks) > 0) {
        fprintf(stderr, "Error: target file is missing %d k-mers from original\n", HASH_COUNT(ks));
        Hs* it;
        int c = 0;
        for (it = ks; it != NULL; it = it->hh.next) {
            char t[k + 1];
            dec(it->key, k, t);
            fprintf(stderr, "  - missing: %s\n", t);
            c++;
        }
        fz = false;
    }

    free_hs(&ks);

    if (fz) {
        fprintf(stdout, "All original k-mers were found\n");
        return 0;
    } else {
        fprintf(stderr, "spectrum not identical\n");
        return -1;
    }
}