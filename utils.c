#include "utils.h"
#include <string.h>
#include <ctype.h>

u64 next_pos(const char *read, u64 start, int k) {
    u64 len = strlen(read);
    int cnt = 0;
    for (u64 i = start; i < len; i++) {
        char c = toupper(read[i]);
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            cnt++;
            if (cnt == k) return i - k + 1;
        } else {
            cnt = 0; // reset cnt if non-ACGT found
        }
    }
    return INF; // no valid k-mer found after start
}

u64 enc(const char *s, int k) {
    u64 h = 0;
    for (int i = 0; i < k; i++) {
        char c = toupper(s[i]);
        if      (c == 'A') h = (h << 2) | 0; // A : 00
        else if (c == 'C') h = (h << 2) | 1; // C : 01
        else if (c == 'G') h = (h << 2) | 2; // G : 10
        else               h = (h << 2) | 3; // T : 11
    }
    return h;
}

void dec(u64 h, int k, char *s) {
    for (int i = 0; i < k; i++) {
        s[k - i - 1] = B[h & 3];
        h >>= 2;
        s[k] = '\0';
    }
}

char* dec_base(u64 h) {
    if      (h == 0) return "A";
    else if (h == 1) return "C";
    else if (h == 2) return "G";
    else             return "T";
    
}

void proc_sq(const char *sq, int k, Hm **km, u64 *id) {
    u64 sq_len = (u64)strlen(sq);
    if (sq_len < (u64)k) return; // too short
    u64 h;
    char s[k + 1]; // buffer for k-mer string
    
    u64 j = next_pos(sq, 0, k);
    if (j != INF) {
        strncpy(s, sq + j, k);
        s[k] = '\0';
        h = enc(s, k);
        if (add_hm(km, h, *id)) (*id)++;
        // rolling hash
        u64 m = (1ULL << (2 * (k - 1))) - 1; // clear two MSBs
        for (++j; j <= sq_len - k; ++j) {
            char c = toupper(sq[j + k - 1]);
            if      (c == 'A')  h = ((h & m) << 2) | 0;
            else if (c == 'C')  h = ((h & m) << 2) | 1;
            else if (c == 'G')  h = ((h & m) << 2) | 2; 
            else if (c == 'T')  h = ((h & m) << 2) | 3;
            else {
                j = next_pos(sq, j, k);
                if (j == INF) break; // no more k-mers in this sequence
                strncpy(s, sq + j, k);
                s[k] = '\0';
                h = enc(s, k);
            }         
            if (add_hm(km, h, *id)) (*id)++;
        }
    }
}

u64 extract(const char* infile, int k, Hm **km, u64 **ka) {
    FILE *fp = fopen(infile, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open file %s\n", infile);
        exit(EXIT_FAILURE);
    }
    printf("file %s opened\n", infile);

    u64 id = 0;
    char *ln = NULL;
    size_t len = 0;

    char *buff = NULL;
    size_t buff_len = 0; // current len of the string in the buffer
    size_t buff_cap = 0; // current allocated capacity for the buffer

    while ((getline(&ln, &len, fp)) != -1) {
        if (ln[0] == '>') {
            if (buff_len > 0) {
                // process the sequence accumulated so far
                proc_sq(buff, k, km, &id);
            }
            buff_len = 0;
        } else {
            ln[strcspn(ln, "\r\n")] = 0; // remove carriage returns and newlines
            size_t ln_len = strlen(ln);
            if (ln_len == 0) continue;
            // check if the buffer has enough capacity
            if (buff_len + ln_len + 1 > buff_cap) {
                size_t ncap = (buff_cap == 0) ? 256 : buff_cap * 2;
                while (ncap < buff_len + ln_len + 1) {
                    ncap *= 2;
                }
                char *nbuff = realloc(buff, ncap);
                if (nbuff == NULL) {
                    fprintf(stderr, "Error: realloc failed for sq_buff\n");
                    exit(EXIT_FAILURE);
                }
                buff = nbuff;
                buff_cap = ncap;
            }
            // append the line to the buffer
            memcpy(buff + buff_len, ln, ln_len);
            buff_len += ln_len;
            buff[buff_len] = '\0';
        }
    }

    // process the very last sequence in the file
    if (buff_len > 0) {
        proc_sq(buff, k, km, &id);
    }

    const u64 N = (u64)HASH_COUNT(*km);
    printf("total unique k-mers = %ld\n", N);
    
    // /* --- display km ---
    printf("key\t\tdec(key)\tval\n");
    printf("-----------------------------\n");
    Hm *s, *tmp;
    HASH_ITER(hh, *km, s, tmp) {
        char t[k + 1];
        dec(s->key, k, t);
        printf("%lu\t\t%s\t\t%lu\n", s->key, t, s->val);
    }
    printf("-----------------------------\n");
    // */

    *ka = malloc(N * sizeof(u64));
    if (ka == NULL) {
        fprintf(stderr, "Error: malloc failed for ka\n");
        exit(EXIT_FAILURE);
    }
    HASH_ITER(hh, *km, s, tmp) {
        (*ka)[s->val] = s->key;
    }

    free(ln);
    free(buff);
    fclose(fp);
    printf("file %s closed\n", infile);

    return N; // number of k-mers
}

u64 step(Hm *km, const u64 *ka, const int k, u64 id, int c, int is_fwd) {
    u64 h = ka[id];
    if (h == INF) return INF;
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return INF;
    if (is_fwd) {
        u64 m = (1ULL << (2 * (k - 1))) - 1;
        h = (h & m) << 2;
        if      (c == 'C')  h |= 1;
        else if (c == 'G')  h |= 2;
        else if (c == 'T')  h |= 3;
    } else {
        h >>= 2;
        if      (c == 'C')  h |= (1ULL << (2 * (k - 1)));
        else if (c == 'G')  h |= (2ULL << (2 * (k - 1)));
        else if (c == 'T')  h |= (3ULL << (2 * (k - 1)));
    }
    u64 nh = find_hm(km, h);
    if (nh == INF) return INF; // no branch to c    
    return nh;
}