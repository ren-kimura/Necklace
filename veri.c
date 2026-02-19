#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "veri.h"
#include "utils.h"

static bool proc_rm(const char* s, int k, bool u_flg, Hs** ks) {
    u64 h = enc(s, k);
    if (!u_flg) h = can(h, k);

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

static bool proc_bp_to(const char* to, int k, bool u_flg, Hs** ks, bool closed) {
    size_t l = strlen(to);
    if (closed) {
        l += (k - 1);
    }
    char* toc = (char*)malloc(l + 1);
    if (!toc) { fprintf(stderr, "Error: malloc failed for toc\n"); return false; }
    size_t i = 0;
    if (closed) {
        int depth = 0;
        Strbld ss; init_strbld(&ss);
        for (size_t x = 0; x < strlen(to); x++) {
            char c = (char)toupper(to[x]);
            if (c == '(') depth++;
            else if (c == ')') depth--;
            else if (isalpha(c) && depth == 0) {
                char tmp[2] = {c, '\0'};
                apnd_strbld(&ss, tmp);
            }
        }
        while (i < (size_t)(k - 1)) {
            int64_t j = ss.len - (k - 1) + i;
            while (j < 0) j += ss.len;
            toc[i] = ss.str[j % ss.len];
            i++;
        }
        free(ss.str);
    }
    memcpy(toc + i, to, strlen(to));
    toc[l] = '\0';

    char buf[k + 1]; // buffer for k-mer
    strncpy(buf, toc, k); buf[k] = '\0';
    u64 h = enc(buf, k);

    if (!proc_rm(buf, k, u_flg, ks)) { // handle the very first k-mer
        free(toc); return false;
    }
    u64 m = (1ULL << (2 * (k - 1))) - 1; // mask

    St s; init_st(&s);
    int depth = 0;

    for (size_t i = k; i < l; i++) {
        char c = (char)toupper(toc[i]);
        if (c == '(') {
            depth++;
            push(&s, h);
        } else if (c == ')') {
            depth--;
            h = pop(&s);
        } else if (isalpha(c)) {
            u64 g = (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : (c == 'T') ? 3 : INF;
            if (g == INF) {
                fprintf(stderr, "Error: non-ACGT appeared in token: %s", toc);
                free(toc); return false;
            }
            h = (h & m) << 2 | g;
            dec(h, k, buf);
            if (!proc_rm(buf, k, u_flg, ks)) {
                free(toc); return false;
            }
        }
    }
    free(toc);
    if (!is_empty_st(&s)) { fprintf(stderr, "Error: stack not empty at the end of verification of bp rep\n"); return false; }
    return true;
}

static bool bveri(const char* ss, int k, bool u_flg, Hs** ks) {
    char* tt = strdup(ss); // Duplicate string for strtok
    if (!tt) {
        perror("strdup failed in bveri");
        return false;
    }

    char* op = tt;
    char* cp = NULL;

    // Find the boundary ",,"
    char* bd = strstr(tt, ",,");
    if (bd != NULL) {
        *bd = '\0'; // Null-terminate the open part
        cp = bd + 2;
    } else {
        // No ",,", check if the whole string starts with ','
        if (tt[0] == ',') {
            cp = tt + 1; // Skip the first comma
            op = NULL; // No open part
        } else {
            // Assumes the whole string is an open path if no boundary and doesn't start with ','
            // open_part remains tt, closed_part remains NULL
        }
    }

    bool vf = true;

    // Process open necklaces (linear root)
    if (op != NULL) {
        char* to = strtok(op, ",");
        while (to != NULL) {
            if (*to != '\0') {
                 if (!proc_bp_to(to, k, u_flg, ks, false)) { // closed = false
                    vf = false;
                    goto clean;
                 }
            }
            to = strtok(NULL, ",");
        }
    }

    // Process closed necklaces (circular root)
    if (cp != NULL) {
        char* to = strtok(cp, ",");
        while (to != NULL) {
             if (*to != '\0') {
                 if (!proc_bp_to(to, k, u_flg, ks, true)) { // closed = true
                    vf = false;
                    goto clean;
                 }
            }
            to = strtok(NULL, ",");
        }
    }

clean:
    free(tt);
    return vf;
}

int veri(const char* of, const char* tf, int k, bool u_flg) {
    fprintf(stdout, "verification mode\n");
    fprintf(stdout, "original file: %s\n", of);
    fprintf(stdout, "target file: %s\n", tf);
    fprintf(stdout, "\n[Step 1/3] extracting k-mers from original file\n");
    Hm *km = NULL;
    u64 *ka = NULL;
    u64 no = extract(of, k, &km, &ka, u_flg);
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

    if (!bveri(sc, k, u_flg, &ks)) { free_hs(&ks); free(sc); return -1; }
    free(sc);

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

int veri_fa(const char *of, const char *tf, int k, bool u_flg) {
    fprintf(stdout, "verification mode\n");
    fprintf(stdout, "original file: %s\n", of);
    fprintf(stdout, "target file: %s\n", tf);

    fprintf(stdout, "\n[Step 1/3] extracting k-mers from original file\n");
    Hm *km = NULL;
    u64 *ka = NULL;
    u64 no1 = extract(of, k, &km, &ka, u_flg);
    Hs *ks = NULL;
    Hm *s, *tmp;
    HASH_ITER(hh, km, s, tmp) { add_hs(&ks, s->key); }
    free(ka); free_hm(&km);
    fprintf(stdout, "%ld k-mers in %s\n", no1, of); 
    
    fprintf(stdout, "\n[Step 2/3] reconstructing k-mers from target file\n");
    km = NULL;
    ka = NULL;
    u64 no2 = extract(tf, k, &km, &ka, u_flg);
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