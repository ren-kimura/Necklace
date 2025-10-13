#include "cov.h"
#include "utils.h"
#include "stat.h"
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>

int bfs(Hm *km, u64 *ka, const u64 *mu, const u64 *mv, u64 *dt, u64 k, u64 N) {
    Q q;
    init_q(&q);
    int found = 0;
    for (u64 u = 0; u < N; u++) {
        if (mu[u] == INF) {
            dt[u] = 0;
            enq(&q, u);
        } else dt[u] = INF;
    }
    while (!is_empty_q(&q)) {
        u64 u = deq(&q);
        for (int i = 0; i < 4; i++) {
            char c = B[i];
            u64 v = step(km, ka, k, u, c, 1);
            if (v == INF) continue;
            if (mv[v] == INF) {
                found = 1;
            } else if (dt[mv[v]] == INF) {
                dt[mv[v]] = dt[u] + 1;
                enq(&q, mv[v]);
            }
        }
    }
    return found;
}

int dfs(Hm *km, u64 *ka, u64 *mu, u64 *mv, u64 *dt, u64 u, u64 k) {
    for (int i = 0; i < 4; i++){
        char c = B[i];
        u64 v = step(km, ka, k, u, c, 1);
        if (v == INF) continue;
        if (mv[v] == INF || (dt[mv[v]] == dt[u] + 1 && 
            dfs(km, ka, mu, mv, dt, mv[v], k))) {
            mu[u] = v;
            mv[v] = u;
            return 1;
        }
    }
    dt[u] = INF; // found no augpath
    return 0;
}

u64 mbm(Hm *km, u64 *ka, u64 *mu, u64 *mv, u64 k, u64 N) {
    u64 *dt = malloc(N * sizeof(u64));
    if (dt == NULL) {
        printf("Error: malloc failed for array dt\n");
        free(dt);
        return INF;
    }
    u64 M = 0;
    for (u64 i = 0; i < N; i++) {
        mu[i] = INF;
        mv[i] = INF;
    }
    while (bfs(km, ka, mu, mv, dt, k, N)) {
        for (u64 u = 0; u < N; u++) {
            if (mu[u] == INF && dfs(km, ka, mu, mv, dt, u, k)) {
                M++; // found an augpath
            }
        }
        prog(M, N, "maximum matching");
    }
    fin("maximum matching");
    printf("matching size %ld\n", M);
    free(dt);
    return M;
}

void decompose(u64 *mu, u64 *mv, VV *cc, VV *pp, u64 N){
    u64 *su = (u64*)malloc(N * sizeof(u64));
    u64 *sv = (u64*)malloc(N * sizeof(u64));
    if (su == NULL || sv == NULL) {
        fprintf(stderr, "Error: malloc failed for su or sv\n");
        free(su); free(sv);
        return;
    }

    for (u64 i = 0; i < N; i++) su[i] = sv[i] = 0;
    for (u64 u = 0; u < N; u++) {
        if (su[u]) continue;
        V fv; init_v(&fv);
        u64 cur = u;
        while (1) {
            if (su[cur]) break;
            su[cur] = 1;
            push_back(&fv, cur);
            u64 next = mu[cur];
            if (next == INF) break;
            sv[next] = 1;
            cur = next;
        }
        if (fv.size > 0 && mu[fv.data[fv.size - 1]] == fv.data[0]) {
            push_backv(cc, fv);
        } else {
            V bv; init_v(&bv);
            cur = u;
            while (1) {
                if (sv[cur]) break;
                sv[cur] = 1;
                push_back(&bv, cur);
                u64 prev = mv[cur];
                if (prev == INF) break;
                su[prev] = 1;
                cur = prev;
            }
            V v; init_v(&v);
            if (bv.size > 0) {
                for (u64 i = bv.size - 1; i > 0; i--) {
                    push_back(&v, bv.data[i]);
                }
            }
            for (u64 i = 0; i < fv.size; i++) {
                push_back(&v, fv.data[i]);
            }
            push_backv(pp, v);
            free_v(&bv);
            free_v(&v);
        }
        prog(u, N, "decomposing");
        free_v(&fv);
    }
    fin("decomposing");
    free(su); free(sv);
}

void dproc_sq(const char *sq, int k, Hm **km, u64 *id, int di, VV *cc, VV *pp) {
    u64 sq_len = (u64)strlen(sq);
    if (sq_len < (u64)k) return; // too short

    u64 m = (1ULL << (2 * (k - 1))) - 1; // mask that clears two MSBs
    V v; init_v(&v);
    
    u64 j = 0;
    while ((j = next_pos(sq, j, k)) != INF) {
        char s[k + 1];
        memcpy(s, sq + j, k);
        s[k] = '\0';
        u64 h = enc(s, k);
        if (di) h = can(h, k);
        if (add_hm(km, h, *id)) {
            push_back(&v, *id);
            (*id)++;
        } else {
            if (v.size > 0) {
                u64 y = find_hm(*km, h);
                if (y == v.data[0]) {
                    push_backv(cc, v);
                } else {
                    push_backv(pp, v);
                }
            }
            init_v(&v);
        }
        
        // rolling hash
        j++;
        while (j <= sq_len - k) {
            char c = toupper(sq[j + k - 1]);
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
                h = ((h & m) << 2) | (c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3);
                if (di) h = can(h, k);
                if (add_hm(km, h, *id)) {
                    push_back(&v, *id);
                    (*id)++;
                } else {
                    if (v.size > 0) {
                        u64 y = find_hm(*km, h);
                        if (y == v.data[0]) {
                            push_backv(cc, v);
                        } else {
                            push_backv(pp, v);
                        }
                    }
                    init_v(&v);
                }
            } else {
                if (v.size > 0) {
                    push_backv(pp, v);
                }
                init_v(&v);
                break;
            }
            j++;
        }
    }
    if (v.size > 0) {
        push_backv(pp, v);
    }
    free_v(&v);
}

u64 dextract(const char* infile, int k, Hm **km, u64 **ka, int di, VV *cc, VV *pp) {
    FILE *fp = fopen(infile, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open file %s\n", infile);
        exit(EXIT_FAILURE);
    }
    printf("file %s opened\n", infile);

    struct stat st;
    stat(infile, &st);
    size_t fs = st.st_size;

    u64 id = 0;
    char *ln = NULL;
    size_t len = 0;

    char *buff = NULL;
    size_t buff_len = 0; // current len of the string in the buffer
    size_t buff_cap = 0; // current allocated capacity for the buffer

    while ((getline(&ln, &len, fp)) != -1) {
        prog(ftell(fp), fs, "extracting k-mers into cycles and paths");
        if (ln[0] == '>') {
            if (buff_len > 0) {
                // process the sequence accumulated so far
                dproc_sq(buff, k, km, &id, di, cc, pp);
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
        dproc_sq(buff, k, km, &id, di, cc, pp);
    }

    fin("extracting k-mers into cycles and paths");
    const u64 N = (u64)HASH_COUNT(*km);
    printf("total unique k-mers = %ld\n", N);
    
    // disp_hm(*km, k);

    *ka = malloc(N * sizeof(u64));
    if (ka == NULL) {
        fprintf(stderr, "Error: malloc failed for ka\n");
        exit(EXIT_FAILURE);
    }
    Hm *s, *tmp;
    HASH_ITER(hh, *km, s, tmp) {
        (*ka)[s->val] = s->key;
    }

    free(ln);
    free(buff);
    fclose(fp);
    printf("file %s closed\n", infile);

    u64 Z = 0;
    for (size_t i = 0; i < cc->size; i++) Z += cc->vs[i].size;
    for (size_t i = 0; i < pp->size; i++) Z += pp->vs[i].size;
    if (N != Z) {
        fprintf(stderr, "Error: total nodes in cycles and paths is not equal to k-mers count\n");
        exit(EXIT_FAILURE);
    }

    return N; // number of k-mers
}

void disp_cp(u64 *ka, VV *cc, VV *pp, int k) {
    printf("---------------------------------\n");
    printf("---cycles---\n");
    for (size_t i = 0; i < cc->size; i++) {
        V* c = &cc->vs[i];
        for (size_t j = 0; j < c->size; j++) {
            char s[k + 1];
            memset(s, 0, k + 1);
            dec(ka[c->data[j]], k, s);
            printf("%s ", s);
        }
        printf("\n");
    }
    printf("---paths---\n");
    for (size_t i = 0; i < pp->size; i++) {
        V* p = &pp->vs[i];
        for (size_t j = 0; j < p->size; j++) {
            char s[k + 1];
            memset(s, 0, k + 1);
            dec(ka[p->data[j]], k, s);
            printf("%s ", s);
        }
        printf("\n");
    }
    printf("---------------------------------\n");
}