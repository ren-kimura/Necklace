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
    }
    s[k] = '\0';
}

char* dec_base(u64 h) {
    if      (h == 0) return "A";
    else if (h == 1) return "C";
    else if (h == 2) return "G";
    else             return "T";
    
}

u64 rc(u64 h, int k)
{
    /* constant time calc of reverse complement */
    /* ref-> https://www.biostars.org/p/113640/#424278 */
    u64 nh = ~h; // take NOT
    // 2-bit-wise left-right flip
    nh = ((nh >> 2 & 0x3333333333333333) | (nh & 0x3333333333333333) << 2);
    nh = ((nh >> 4 & 0x0F0F0F0F0F0F0F0F) | (nh & 0x0F0F0F0F0F0F0F0F) << 4);
    nh = ((nh >> 8 & 0x00FF00FF00FF00FF) | (nh & 0x00FF00FF00FF00FF) << 8);
    nh = ((nh >> 16 & 0x0000FFFF0000FFFF) | (nh & 0x0000FFFF0000FFFF) << 16);
    nh = ((nh >> 32 & 0x00000000FFFFFFFF) | (nh & 0x00000000FFFFFFFF) << 32);
    
    return (nh >> (2 * (32 - k))); // extract MSBs
}

u64 can(u64 h, int k) {
    u64 rh = rc(h, k);
    return (h <= rh) ? h : rh;
}

void proc_sq(const char *sq, int k, Hm **km, u64 *id, int di) {
    u64 sq_len = (u64)strlen(sq);
    if (sq_len < (u64)k) return; // too short
    
    char s[k + 1]; // buffer for kmer
    u64 m = (1ULL << (2 * (k - 1))) - 1; // mask that clears two MSBs
    
    u64 j = next_pos(sq, 0, k); // look for the first kmer
    while (j != INF) {
        strncpy(s, sq + j, k);
        s[k] = '\0';
        u64 h = enc(s, k);
        if (di) h = can(h, k);
        if (add_hm(km, h, *id)) (*id)++;
        
        // rolling hash
        while (++j <= sq_len - k) {
            char c = toupper(sq[j + k - 1]);
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
                h = ((h & m) << 2) | (c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3);
                if (di) h = can(h, k);
                if (add_hm(km, h, *id)) (*id)++;
            } else {
                // if interrupted by a non-ACGT
                j = next_pos(sq, j, k);
                goto nxt;
            }
        }
        j = INF; // reached the end of sq
    nxt:;
    }
}

u64 extract(const char* infile, int k, Hm **km, u64 **ka, int di) {
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
                proc_sq(buff, k, km, &id, di);
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
        proc_sq(buff, k, km, &id, di);
    }

    const u64 N = (u64)HASH_COUNT(*km);
    printf("total unique k-mers = %ld\n", N);
    
    /* --- display km ---
    printf("key\t\tdec(key)\tval\n");
    printf("-----------------------------\n");
    Hm *s, *tmp;
    HASH_ITER(hh, *km, s, tmp) {
        char t[k + 1];
        dec(s->key, k, t);
        printf("%lu\t\t%s\t\t%lu\n", s->key, t, s->val);
    }
    printf("-----------------------------\n");
    */

    *ka = malloc(N * sizeof(u64));
    if (ka == NULL) {
        fprintf(stderr, "Error: malloc failed for ka\n");
        exit(EXIT_FAILURE);
    }
    Hm *s, *tmp; // comment out if [display km] activated
    HASH_ITER(hh, *km, s, tmp) {
        (*ka)[s->val] = s->key;
    }

    free(ln);
    free(buff);
    fclose(fp);
    printf("file %s closed\n", infile);

    return N; // number of k-mers
}

u64 step(Hm *km, const u64 *ka, const int k, u64 id, int c, bool is_fwd) {
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
    u64 nid = find_hm(km, h);
    if (nid == INF) return INF; // no branch to c    
    return nid;
}

u64 bstep(Hm *km, const u64 *ka, const int k, u64 id, int c, bool is_fwd, bool fromc, bool toc) {
    u64 h = ka[id]; // h is already canonical
    if (!fromc) h = rc(h, k); // take rc of h if starts from non-canonical

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
    u64 ch = can(h, k); // take can of h
    u64 nid = find_hm(km, ch);
    if (nid == INF) return INF; // no branch to c from the selected side
    if (ch != rc(ch, k)) {
        if ((ch == h) == toc) return nid; // positive if not self-complement and arrived at the selected side of ch
    } else {
        return nid; // positive regardless of toc if self-complement 
    }
    return INF;
}