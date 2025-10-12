#include "out.h"
#include <string.h>

//---FSt---
void init_fst(FSt *s) {
    s->top = NULL;
}

int is_empty_fst(FSt *s) {
    return (s->top == NULL);
}

void push_fst(FSt *s, u64 cur, u64 idx, int b_idx) {
    F *nf = (F*)malloc(sizeof(F));
    if (nf == NULL) {
        fprintf(stderr, "Error: malloc failed for new frame\n");
        exit(EXIT_FAILURE);
    }
    nf->cur = cur;
    nf->idx = idx;
    nf->b_idx = b_idx;
    nf->next = s->top;
    s->top = nf;
}

void pop_fst(FSt *s) {
    if (is_empty_fst(s)) {
        fprintf(stderr, "Error: pop_fst called on an empty stack\n");
        return;
    }
    F *tmp = s->top;
    s->top = s->top->next;
    free(tmp);
}

//---FbSt---
void init_fbst(FbSt *s) {
    s->top = NULL;
}

int is_empty_fbst(FbSt *s) {
    return (s->top == NULL);
}

void push_fbst(FbSt *s, u64 cur, u64 idx, int b_idx, bool wrn) {
    Fb *nfb = (Fb*)malloc(sizeof(Fb));
    if (nfb == NULL) {
        fprintf(stderr, "Error: malloc failed for new frame\n");
        exit(EXIT_FAILURE);
    }
    nfb->cur = cur;
    nfb->idx = idx;
    nfb->b_idx = b_idx;
    nfb->wrn = wrn;
    nfb->next = s->top;
    s->top = nfb;
}

void pop_fbst(FbSt *s) {
    if (is_empty_fbst(s)) {
        fprintf(stderr, "Error: pop_fbst called on an empty stack\n");
        return;
    }
    Fb *tmp = s->top;
    s->top = s->top->next;
    free(tmp);
}

//---Strbld---
void init_strbld(Strbld *s) {
    s->cap = 64;
    s->len = 0;
    s->str = (char*)malloc(s->cap);
    if (s->str == NULL) {
        fprintf(stderr, "Error: malloc failed for s.str\n");
        exit(EXIT_FAILURE);
    }
    s->str[0] = '\0';
}

void apnd_strbld(Strbld *s, const char* t) {
    size_t lt = strlen(t);
    if (s->len + lt + 1 > s->cap) {
        size_t ncap = s->cap;
        while (ncap < s->len + lt + 1) ncap *= 2;
        char* nstr = (char*)realloc(s->str, ncap);
        if (nstr == NULL) {
            fprintf(stderr, "Error: realloc failed in apnd_strbld\n");
            exit(EXIT_FAILURE);
        }
        s->str = nstr;
        s->cap = ncap;
    }
    memcpy(s->str + s->len, t, lt);
    s->len += lt;
    s->str[s->len] = '\0';
}

//---Rep---
void init_rep(Rep *r) {
    r->str = NULL;
    r->arr = NULL;
}

void free_rep(Rep *r) {
    free(r->str);
    free(r->arr);
}

//---Flat---
Rep flat(u64 *ka, VV *cc, VV *pp, int k) {
    Rep r; init_rep(&r);
    size_t ll = 0;
    for (size_t i = 0; i < cc->size; i++) ll += cc->vs[i].size;
    for (size_t i = 0; i < pp->size; i++) ll += k - 1 + pp->vs[i].size;
    ll += cc->size + pp->size + 2; 
    r.str = malloc(ll);
    if (r.str == NULL) {
        fprintf(stderr, "Error: malloc failed for r.str\n");
        exit(EXIT_FAILURE);
    }
    char *ptr = r.str;
    ptr[0] = '\0';

    for (size_t i = 0; i < cc->size; i++) {
        V *cur = &cc->vs[i];
        for (size_t j = 0; j < cur->size; j++) {
            ptr += sprintf(ptr, "%s", dec_base(ka[cur->data[j]] % 4));
        }
        ptr += sprintf(ptr, ",");
    }
    ptr += sprintf(ptr, ",");

    for (size_t i = 0; i < pp->size; i++) {
        V *cur = &pp->vs[i];
        char s[k + 1];
        dec(ka[cur->data[0]], k, s);
        ptr += sprintf(ptr, "%s", s);
        for (size_t j = 1; j < cur->size; j++) {
            ptr += sprintf(ptr, "%s", dec_base(ka[cur->data[j]] % 4));
        }
        ptr += sprintf(ptr, ",");
    }

    size_t fl = strlen(r.str);
    if (fl > 0 && r.str[fl - 1] == ',') {
        r.str[fl - 1] = '\0';
    }
     if (fl > 1 && r.str[fl - 2] == ',' && r.str[fl - 1] == '\0') {
        r.str[fl - 2] = '\0';
    }
    return r;
}

V* findc(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *ino, bool *vis, u64 from) {
    V* nc = (V*)malloc(sizeof(V));
    if (nc == NULL) return NULL;
    init_v(nc);

    FSt s;
    init_fst(&s);

    vis[from] = true;
    push_fst(&s, from, 0, 0);

    while (!is_empty_fst(&s)) {
        if (s.top->idx >= (u64)pp->vs[s.top->cur].size) {
            pop_fst(&s);
            if (nc->size > 0) pop_back(nc);
            continue;
        }

        if (s.top->idx == 0 && s.top->b_idx == 0) push_back(nc, s.top->cur);
        
        if (s.top->b_idx < 4) {
            char c = B[s.top->b_idx];
            s.top->b_idx++;
            
            u64 cur = pp->vs[s.top->cur].data[s.top->idx];
            u64 nxt = step(km, ka, k, cur, c, 1);
            if (nxt == INF) continue;

            u64 ni = find_hm(hd, nxt);

            if (ni == INF || ino[ni] || vis[ni]) continue;

            if (ni == from) return nc; // cycle found
            
            vis[ni] = true;
            push_fst(&s, ni, 0, 0);
        } else {
            s.top->idx++;
            s.top->b_idx = 0;
        }
    }
    free_v(nc);
    free(nc);
    return NULL; // no cycle found
}

u64 nf(bool *v, u64 l) {
    u64 n = 0;
    for (u64 i = 0; i < l; i++) {
        if (v[i]) continue; // skip if true
        n++;
    }
    return n;
}

Rep ptr(Hm *km, u64 *ka, VV *cc, VV *pp, int k) {
    Hs *rp = NULL;
    Hm *hd = NULL; 
    for (u64 i = 0; i < pp->size; i++) {
        bool is_r = true;
        u64 h = pp->vs[i].data[0];
        for (int j = 0; j < 4; j++) {
            if (step(km, ka, k, h, B[j], 0) != INF) {
                is_r = false;
                break;
            }
        }
        if (is_r) add_hs(&rp, i);
        else      add_hm(&hd, h, i);
    }
    Rep r; init_rep(&r);
    r.arr = (u64*)malloc(pp->size * sizeof(u64));
    if (r.arr == NULL) {
        fprintf(stderr, "Error: malloc failed for r.arr\n");
        exit(EXIT_FAILURE);
    }
    for (u64 i = 0; i < pp->size; i++) { r.arr[i] = -1; }
    printf("# of (root paths, non-root paths) = (%d, %d)\n", HASH_COUNT(rp), HASH_COUNT(hd));

    Q q; init_q(&q);
    Vb o; init_vb(&o); // i-th item is false if i-th one is a cycle, true if it is a path
    V oi; init_v(&oi); // cycle/path id of i-th item
    bool ino[pp->size]; // false if the path is not in o yet, true otherwise
    for (size_t i = 0; i < pp->size; i++) { ino[i] = false; }
    
    u64 d = 0; // current cumulative distance
    u64 f = 0; // r.arr is initialized up to this idx

    Hs *s, *tmp;
    HASH_ITER(hh, rp, s, tmp) { // bfs from root paths
        enq(&q, s->key);
        push_backb(&o, true); // added a path
        push_back(&oi, s->key);
        ino[s->key] = true;
        r.arr[s->key] = d;
        d += k; // ',' and k-1 chars in the front

        while (!is_empty_q(&q)) {
            u64 i = deq(&q);
            for (size_t j = 0; j < pp->vs[i].size; j++) {
                u64 cur = pp->vs[i].data[j];
                for (int c = 0; c < 4; c++) {
                    u64 nxt = step(km, ka, k, cur, B[c], 1);
                    if (nxt == INF) continue;
                    u64 ni = find_hm(hd, nxt);
                    if (ni == INF) continue;
                    if (ino[ni]) continue;
                    enq(&q, ni);
                    push_backb(&o, true); // added a path
                    push_back(&oi, ni);
                    ino[ni] = true;
                    r.arr[f++] = d;
                } d++;
            } d++;
        }
    }

    for (size_t i = 0; i < cc->size; i++) { // bfs from cycles
        push_backb(&o, false); // added a cycle
        push_back(&oi, i);
        d++; // need flag '*' for a cycle to tell apart from a path
        
        for (size_t j = 0; j < cc->vs[i].size; j++) {
            u64 cur = cc->vs[i].data[j];
            for (int c = 0; c < 4; c++) {
                u64 nxt = step(km, ka, k, cur, B[c], 1);
                if (nxt == INF) continue;
                u64 ni = find_hm(hd, nxt);
                if (ni == INF) continue;
                if (ino[ni]) continue;
                enq(&q, ni);
                push_backb(&o, true); // added a path
                push_back(&oi, ni);
                ino[ni] = true;
                r.arr[f++] = d;
            } d++;
        } d++;
        
        while (!is_empty_q(&q)) {
            u64 l = deq(&q);
            for (size_t j = 0; j < pp->vs[l].size; j++) {
                u64 cur = pp->vs[l].data[j];
                for (int c = 0; c < 4; c++) {
                    u64 nxt = step(km, ka, k, cur, B[c], 1);
                    if (nxt == INF) continue;
                    u64 ni = find_hm(hd, nxt);
                    if (ni == INF) continue;
                    if (ino[ni]) continue;
                    enq(&q, ni);
                    push_backb(&o, true); // added a path
                    push_back(&oi, ni);
                    ino[ni] = true;
                    r.arr[f++] = d;
                } d++;
            } d++;
        }
    }

    u64 l;
    while ((l = nf(ino, pp->size))) {
        printf("%ld paths remaining\n", l);
        for (size_t i = 0; i < pp->size; i++) {
            if (ino[i]) continue;
            bool vis[pp->size];
            for (size_t i = 0; i < pp->size; i++) { vis[i] = false; }
            V *pis = findc(km, ka, hd, k, pp, ino, vis, i);
            if (pis != NULL && pis->size != 0) {
                printf("\rfound a cycle of %ld paths\n", pis->size);
                V nc; init_v(&nc);
                for (size_t j = 0; j < pis->size; j++) {
                    u64 pi = pis->data[j];
                    u64 fpi = pis->data[(j + 1) % pis->size];
                    V *p = &pp->vs[pi];
                    for (size_t x = 0; x < p->size; x++) {
                        u64 cur = p->data[x];
                        push_back(&nc, cur);
                        for (int c = 0; c < 4; c++) {
                            u64 nxt = step(km, ka, k, cur, B[c], 1);
                            if (nxt == INF) continue;
                            if (nxt == pp->vs[fpi].data[0]) {
                                // trim cycle from path, and renew hd
                                del_hm(&hd, pp->vs[pi].data[0]);
                                for (size_t y = 0; y <= x; y++) {
                                    pp->vs[pi].data[y] = pp->vs[pi].data[y + x + 1];
                                }
                                pp->vs[pi].size -= x + 1;
                                add_hm(&hd, pp->vs[pi].data[0], pi);
                                goto ej;
                            }
                        }
                    }
                    ej:;
                }
                push_backv(cc, nc);
                push_backb(&o, false);
                push_back(&oi, cc->size - 1);
                d++; // need '*' for cycle
                for (size_t j = 0; j < nc.size; j++) {
                    u64 cur = nc.data[j];
                    for (int c = 0; c < 4; c++) {
                        u64 nxt = step(km, ka, k, cur, B[c], 1);
                        if (nxt == INF) continue;
                        u64 ni = find_hm(hd, nxt);
                        if (ni == INF) continue;
                        if (ino[ni]) continue;
                        enq(&q, ni);
                        push_backb(&o, true);
                        push_back(&oi, ni);
                        ino[ni] = true;
                        r.arr[f++] = d;
                    } d++;
                } d++;
                while (!is_empty_q(&q)) {
                    u64 j = deq(&q);
                    for (size_t x = 0; x < pp->vs[j].size; x++) {
                        u64 cur = pp->vs[j].data[x];
                        for (int c = 0; c < 4; c++) {
                            u64 nxt = step(km, ka, k, cur, B[c], 1);
                            if (nxt == INF) continue;
                            u64 ni = find_hm(hd, nxt);
                            if (ni == INF) continue;
                            if (ino[ni]) continue;
                            enq(&q, ni);
                            push_backb(&o, true);
                            push_back(&oi, ni);
                            ino[ni] = true;
                            r.arr[f++] = d;
                        } d++;
                    } d++;
                }
                free_v(&nc);
                free_v(pis);
                free(pis);
                break;
            } else {
                free_v(pis);
                free(pis);
            }
        }
    }

    l = nf(ino, pp->size);
    if (l == 0) printf("All paths are pointed\n");
    else        printf("%ld path(s) are not pointed\n", l);

    // take difference of r.arr
    for (u64 i = pp->size; i > 0; --i) r.arr[i] -= r.arr[i - 1];

    size_t ll = 0;
    for (size_t i = 0; i < cc->size; i++) ll += cc->vs[i].size + 1; // # of nodes plus a flag '*'
    for (size_t i = 0; i < pp->size; i++) ll += pp->vs[i].size; // # of nodes
    for (size_t i = 0; i < HASH_COUNT(rp); i++) ll += k; // # of root paths: ',' and k - 1 chars in the front
    ll += cc->size + pp->size + 2; 
    r.str = malloc(ll);
    if (r.str == NULL) {
        fprintf(stderr, "Error: malloc failed for r.str\n");
        exit(EXIT_FAILURE);
    }
    char *ptr = r.str;
    ptr[0] = '\0';

    for (size_t i = 0; i < o.size; i++) {
        if (o.data[i]) { // path
            u64 id = oi.data[i];
            int isrp = find_hs(rp, id);
            V *cur = &pp->vs[id];
            if (isrp) { // root path
                char s[k];
                dec((ka[cur->data[0]] >> 2), k - 1, s);
                ptr += sprintf(ptr, ",%s", s);
            }
            for (size_t j = 0; j < cur->size; j++) {
                ptr += sprintf(ptr, "%s", dec_base(ka[cur->data[j]] % 4));
            }
        } else { // cycle
            V *cur = &cc->vs[oi.data[i]];
            ptr += sprintf(ptr, "*");
            for (size_t j = 0; j < cur->size; j++) {
                ptr += sprintf(ptr, "%s", dec_base(ka[cur->data[j]] % 4));
            }
        }
        ptr += sprintf(ptr, ","); // delim
    }

    size_t fl = strlen(r.str);
    if (fl > 0 && r.str[fl - 1] == ',') {
        r.str[fl - 1] = '\0';
    }
     if (fl > 1 && r.str[fl - 2] == ',' && r.str[fl - 1] == '\0') {
        r.str[fl - 2] = '\0';
    }
    return r;
}

char* subt(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *vis, u64 rpi) {
    Strbld sb; init_strbld(&sb);
    FbSt s; init_fbst(&s);

    vis[rpi] = true;
    push_fbst(&s, rpi, 0, 0, false);
    while (!is_empty_fbst(&s)) {
        Fb* f = s.top;
        V* p = &pp->vs[f->cur];
        if (f->idx >= p->size) {
            pop_fbst(&s);
            if (!is_empty_fbst(&s)) {
                apnd_strbld(&sb, ")");
            }
            continue;
        }
        u64 cur = p->data[f->idx];
        if (!f->wrn) {
            apnd_strbld(&sb, dec_base(ka[cur] % 4));
            f->wrn = true;
            f->b_idx = 0;
        }
        bool d = false;
        for (; f->b_idx < 4; f->b_idx++) {
            u64 nxt = step(km, ka, k, cur, B[f->b_idx], 1);
            if (nxt == INF) continue;
            u64 cpi = find_hm(hd, nxt);
            if (cpi == INF) continue;
            if (vis[cpi]) continue;
            vis[cpi] = true;
            f->b_idx++;
            apnd_strbld(&sb, "(");
            push_fbst(&s, cpi, 0, 0, false);
            d = true;
            break;
        }
        if (!d) {
            f->idx++;
            f->wrn = false;
        }
    }
    return sb.str;
}

Rep bp(Hm *km, u64 *ka, VV *cc, VV *pp, int k) {
    Hs *rp = NULL;
    Hm *hd = NULL; 
    for (u64 i = 0; i < pp->size; i++) {
        bool is_r = true;
        u64 h = pp->vs[i].data[0];
        for (int j = 0; j < 4; j++) {
            if (step(km, ka, k, h, B[j], 0) != INF) {
                is_r = false;
                break;
            }
        }
        if (is_r) add_hs(&rp, i);
        else      add_hm(&hd, h, i);
    }
    Rep r; init_rep(&r);
    printf("# of (root paths, non-root paths) = (%d, %d)\n", HASH_COUNT(rp), HASH_COUNT(hd));
    
    Strbld sb;
    init_strbld(&sb);

    bool vis[pp->size]; // 1: path already embedded
    for (size_t i = 0; i < pp->size; i++) vis[i] = false;

    Hs *s, *tmp;
    HASH_ITER(hh, rp, s, tmp) { // dfs from root paths
        char t[k]; // buffer for dec of head
        dec(ka[pp->vs[s->key].data[0]], k - 1, t);
        char* ss = subt(km, ka, hd, k, pp, vis, s->key);
        apnd_strbld(&sb, t);
        apnd_strbld(&sb, ss);
        apnd_strbld(&sb, ",");
        free(ss);
    }
    apnd_strbld(&sb, ","); // extra delim between open and closed necklaces
    for (size_t i = 0; i < cc->size; i++) {
        for (size_t j = 0; j < cc->vs[i].size; j++) {
            u64 cur = cc->vs[i].data[j];
            apnd_strbld(&sb, dec_base(ka[cur] % 4));
            for (int c = 0; c < 4; c++) {
                u64 nxt = step(km, ka, k, cur, B[c], 1);
                if (nxt == INF) continue;
                u64 ni = find_hm(hd, nxt);
                if (ni == INF) continue;
                if (vis[ni]) continue;
                char* ss = subt(km, ka, hd, k, pp, vis, ni);
                apnd_strbld(&sb, "(");
                apnd_strbld(&sb, ss);
                apnd_strbld(&sb, ")");
                free(ss);
            }
        }
        apnd_strbld(&sb, ",");
    }

    u64 l;
    while ((l = nf(vis, pp->size))) {
        printf("\r%ld paths remaining", l);
        bool fndc = false;
        for (u64 i = 0; i < pp->size; i++) {
            if (vis[i]) continue;
            bool tvis[pp->size];
            for (size_t j = 0; j < pp->size; j++) tvis[j] = false;
            V* pis = findc(km, ka, hd, k, pp, vis, tvis, i);
            if (pis != NULL && pis->size != 0) {
                fndc = true;
                V nc; init_v(&nc);
                for (size_t j = 0; j < pis->size; j++) {
                    u64 pi = pis->data[j];
                    u64 fpi = pis->data[(j + 1) % pis->size];
                    V *p = &pp->vs[pi];
                    for (size_t x = 0; x < p->size; x++) {
                        u64 cur = p->data[x];
                        push_back(&nc, cur);
                        u64 nxt_onp = (x + 1 < p->size) ? p->data[x + 1] : INF;
                        for (int c = 0; c < 4; c++) {
                            u64 nxt_br = step(km, ka, k, cur, B[c], 1);
                            if (nxt_br == INF || nxt_br == nxt_onp) continue;
                            if (nxt_br == pp->vs[fpi].data[0]) {
                                // trim cycle from path, and renew hd
                                del_hm(&hd, pp->vs[pi].data[0]);
                                for (size_t y = 0; y <= x; y++) {
                                    pp->vs[pi].data[y] = pp->vs[pi].data[y + x + 1];
                                }
                                pp->vs[pi].size -= x + 1;
                                add_hm(&hd, pp->vs[pi].data[0], pi);
                                goto ej;
                            }
                        }
                    }
                    ej:;
                }
                push_backv(cc, nc);
                for (size_t j = 0; j < nc.size; j++) {
                    u64 cur = nc.data[j];
                    for (int c = 0; c < 4; c++) {
                        u64 nxt = step(km, ka, k, cur, B[c], 1);
                        if (nxt == INF) continue;
                        u64 ni = find_hm(hd, nxt);
                        if (ni == INF) continue;
                        if (vis[ni]) continue;
                        char* ss = subt(km, ka, hd, k, pp, vis, ni);
                        apnd_strbld(&sb, "(");
                        apnd_strbld(&sb, ss);
                        apnd_strbld(&sb, ")");
                        free(ss);
                    }
                }
                free_v(&nc);
                free_v(pis);
                free(pis);
                break;
            } else {
                free_v(pis);
                free(pis);
            }
        }
        if (nf(vis, pp->size) == l) {
            fprintf(stderr, "Error: no new cycle found in remaining paths\n");
            exit(EXIT_FAILURE);
        }
    }
    if (sb.len > 0 && sb.str[sb.len - 1] == ',') {
        sb.len--;
        sb.str[sb.len] = '\0';
    }
    r.str = sb.str;
    return r;
}