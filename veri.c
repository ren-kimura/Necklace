#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include "veri.h"
#include "utils.h"
#include "stat.h"

static bool proc_rm_from_u64(u64 h, int k, bool u_flg, Hs **ks);

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

static bool proc_rm_rolling(const char *sq, size_t sq_len, int k, bool u_flg, Hs **ks) {
    u64 m = (1ULL << (2 * (k - 1))) - 1; //
    u64 j = next_pos(sq, 0, k); //
    
    while (j != INF && j <= sq_len - k) {
        char s[k + 1];
        strncpy(s, sq + j, k); s[k] = '\0';
        u64 h = enc(s, k); //

        if (!proc_rm_from_u64(h, k, u_flg, ks)) return false;

        while (++j <= sq_len - k) {
            char c = toupper(sq[j + k - 1]);
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
                h = ((h & m) << 2) | (c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3); //
                if (!proc_rm_from_u64(h, k, u_flg, ks)) return false;
            } else {
                j = next_pos(sq, j, k);
                break;
            }
        }
    }
    return true;
}

static bool proc_rm_from_u64(u64 h, int k, bool u_flg, Hs **ks) {
    if (!u_flg) h = can(h, k); //
    if (find_hs(*ks, h)) { //
        del_hs(ks, h); //
        return true;
    } else {
        char t[k + 1]; dec(h, k, t); //
        fprintf(stderr, "Error: target file has extra/duplicate k-mer: %s\n", t);
        return false;
    }
}

int veri_fa(const char *of, const char *tf, int k, bool u_flg) {
    fprintf(stdout, "verification mode (FASTA vs FASTA)\n");

    Hm *km = NULL; u64 *ka = NULL;
    extract(of, k, &km, &ka, u_flg);
    
    Hs *ks = NULL;
    Hm *s, *tmp_hm;
    HASH_ITER(hh, km, s, tmp_hm) { add_hs(&ks, s->key); }
    free(ka); free_hm(&km);

    fprintf(stdout, "\n[Step 2/3] verifying target k-mers with rolling hash\n");
    FILE *fp = fopen(tf, "rb");
    if (!fp) { free_hs(&ks); return -1; }

    struct stat st; stat(tf, &st);
    size_t fs = st.st_size;

    char *ln = NULL; size_t len = 0;
    char *buff = NULL; size_t buff_len = 0; size_t buff_cap = 0;
    bool success = true;

    while ((getline(&ln, &len, fp)) != -1) {
        prog(ftell(fp), fs, "verifying"); //
        if (ln[0] == '>') {
            if (buff_len >= (size_t)k) {
                if (!proc_rm_rolling(buff, buff_len, k, u_flg, &ks)) {
                    success = false; break;
                }
            }
            buff_len = 0;
        } else {
            ln[strcspn(ln, "\r\n")] = 0;
            size_t ln_len = strlen(ln);
            if (ln_len == 0) continue;
            
            if (buff_len + ln_len + 1 > buff_cap) {
                size_t ncap = (buff_cap == 0) ? 1024 : buff_cap * 2;
                while (ncap < buff_len + ln_len + 1) ncap *= 2;
                buff = realloc(buff, ncap);
                buff_cap = ncap;
            }
            memcpy(buff + buff_len, ln, ln_len);
            buff_len += ln_len;
            buff[buff_len] = '\0';
        }
    }
    // handle the last buffer
    if (success && buff_len >= (size_t)k) {
        success = proc_rm_rolling(buff, buff_len, k, u_flg, &ks);
    }

    free(ln); free(buff); fclose(fp);
    if (!success) { free_hs(&ks); return -1; }

    fprintf(stdout, "\n[Step 3/3] final check\n");
    if (HASH_COUNT(ks) == 0) {
        fprintf(stdout, "All original k-mers found exactly once.\n");
        free_hs(&ks); return 0;
    } else {
        fprintf(stderr, "Error: %d k-mers missing.\n", HASH_COUNT(ks));
        free_hs(&ks); return -1;
    }
}