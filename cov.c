#include "cov.h"
#include "utils.h"
#include "stat.h"
#include <string.h>

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