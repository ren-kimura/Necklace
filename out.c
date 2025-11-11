#include "out.h"
#include "stat.h"
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

//---FcSt---
void init_fcst(FcSt *s) {
    s->top = NULL;
}

int is_empty_fcst(FcSt *s) {
    return (s->top == NULL);
}

void push_fcst(FcSt *s, u64 nid, int side, bool wrn, int b_idx, int nc) {
    Fc *nfc = (Fc*)malloc(sizeof(Fc));
    if (nfc == NULL) {
        fprintf(stderr, "Error: malloc failed for new frame\n");
        exit(EXIT_FAILURE);
    }
    nfc->nid = nid;
    nfc->side = side;
    nfc->wrn = wrn;
    nfc->b_idx = b_idx;
    nfc->nc = nc;
    nfc->next = s->top;
    s->top = nfc;
}

void pop_fcst(FcSt *s) {
    if (is_empty_fcst(s)) {
        fprintf(stderr, "Error: pop_fcst called on an empty stack\n");
        return;
    }
    Fc *tmp = s->top;
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
    Strbld sb; init_strbld(&sb);

    printf("# of (cycles, paths) = (%ld, %ld)\n", cc->size, pp->size);

    for (size_t i = 0; i < cc->size; i++) {
        V *cur = &cc->vs[i];
        for (size_t j = 0; j < cur->size; j++) {
            apnd_strbld(&sb, dec_base(ka[cur->data[j]] % 4));
        }
        apnd_strbld(&sb, ",");
    }
    apnd_strbld(&sb, ",");

    for (size_t i = 0; i < pp->size; i++) {
        V *cur = &pp->vs[i];
        char s[k + 1];
        dec(ka[cur->data[0]], k, s);
        apnd_strbld(&sb, s);
        for (size_t j = 1; j < cur->size; j++) {
            apnd_strbld(&sb, dec_base(ka[cur->data[j]] % 4));
        }
        apnd_strbld(&sb, ",");
    }

    r.str = sb.str;
    size_t fl = strlen(r.str);
    if (fl > 0 && r.str[fl - 1] == ',') {
        r.str[fl - 1] = '\0';
    }
     if (fl > 1 && r.str[fl - 2] == ',' && r.str[fl - 1] == '\0') {
        r.str[fl - 2] = '\0';
    }
    return r;
}

Rep flat_w(W *w) {
    Rep r; init_rep(&r);
    Strbld sb; init_strbld(&sb);

    size_t ncc = 0;
    if (w->cc) for (ncc = 0; w->cc[ncc] != NULL; ncc++);
    size_t npp = 0;
    if (w->pp) for (npp = 0; w->pp[npp] != NULL; npp++);

    printf("# of (cycles, paths) = (%ld, %ld)\n", ncc, npp);

    if (w->cc) {
        for (size_t i = 0; w->cc[i] != NULL; i++) {
            apnd_strbld(&sb, w->cc[i]);
            apnd_strbld(&sb, ",");
        }
    }
    apnd_strbld(&sb, ",");
    if (w->pp) {
        for (size_t i = 0; w->pp[i] != NULL; i++) {
            apnd_strbld(&sb, w->pp[i]);
            apnd_strbld(&sb, ",");
        }
    }
    r.str = sb.str;
    size_t fl = sb.len;

    if (fl > 0 && r.str[fl - 1] == ',') {
        r.str[fl - 1] = '\0';
    }
    if (fl > 1 && r.str[fl - 2] == ',' && r.str[fl - 1] == '\0') {
        r.str[fl - 2] = '\0';
    }
    return r;
}

Rep bflat(u64 *ka, VV *cc, VV *pp, VVb *ccb, VVb *ppb, int k) {
    Rep r; init_rep(&r);
    Strbld sb; init_strbld(&sb);

    printf("# of (cycles, paths) = (%ld, %ld)\n", cc->size, pp->size);

    for (size_t i = 0; i < cc->size; i++) {
        V *cur = &cc->vs[i]; Vb *curb = &ccb->vs[i];
        for (size_t j = 0; j < cur->size; j++) {
            u64 h = ka[cur->data[j]];
            if (curb->data[j] == false) h = rc(h, k);
            apnd_strbld(&sb, dec_base(h % 4));
        }
        apnd_strbld(&sb, ",");
    }
    apnd_strbld(&sb, ",");

    for (size_t i = 0; i < pp->size; i++) {
        V *cur = &pp->vs[i]; Vb *curb = &ppb->vs[i];
        char s[k + 1];
        u64 h = ka[cur->data[0]];
        if (curb->data[0] == false) h = rc(h, k);
        dec(h, k, s);
        apnd_strbld(&sb, s);
        for (size_t j = 1; j < cur->size; j++) {
            u64 h = ka[cur->data[j]];
            if (curb->data[j] == false) h = rc(h, k);
            apnd_strbld(&sb, dec_base(h % 4));
        }
        apnd_strbld(&sb, ",");
    }

    r.str = sb.str;
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
            if (ni == INF || ino[ni]) continue;
            if (ni == from) return nc; // cycle found
            if (vis[ni]) continue;            
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

//---Pointer---
Rep ptr(Hm *km, u64 *ka, VV *cc, VV *pp, int k) {
    size_t D = HASH_COUNT(km) + cc->size + pp->size;
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

    u64 *hi = (u64*)calloc(pp->size, sizeof(u64));
    if (hi == NULL) {
        fprintf(stderr, "Error: calloc failed for hi\n");
        exit(EXIT_FAILURE);
    }
    
    u64 d = 0; // current cumulative distance
    u64 f = 0; // r.arr is initialized up to this idx

    Hs *s, *tmp;
    HASH_ITER(hh, rp, s, tmp) { // bfs from root paths
        enq(&q, s->key);
        push_backb(&o, true); // added a path
        push_back(&oi, s->key);
        ino[s->key] = true;
        r.arr[f++] = d;
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
                } d++; prog(d, D, "pointing");
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
            } d++; prog(d, D, "pointing");
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
                } d++; prog(d, D, "pointing");
            } d++;
        }
    }

    u64 l;
    while ((l = nf(ino, pp->size))) {
        printf("\n%ld paths remaining", l);
        bool fndc = false;
        for (size_t i = 0; i < pp->size; i++) {
            if (ino[i]) continue;
            bool vis[pp->size];
            for (size_t i = 0; i < pp->size; i++) { vis[i] = false; }
            V *pis = findc(km, ka, hd, k, pp, ino, vis, i);
            if (pis != NULL && pis->size != 0) {
                fndc = true;
                printf("\rfound a cycle of %ld paths\n", pis->size);
                V nc; init_v(&nc);
                for (size_t j = 0; j < pis->size; j++) {
                    u64 pi = pis->data[j];
                    u64 fpi = pis->data[(j + 1) % pis->size];
                    V *p = &pp->vs[pi];
                    for (size_t x = hi[pi]; x < p->size; x++) {
                        u64 cur = p->data[x];
                        push_back(&nc, cur);
                        u64 nxt_onp = (x + 1 < p->size) ? p->data[x + 1] : INF;
                        for (int c = 0; c < 4; c++) {
                            u64 nxt_br = step(km, ka, k, cur, B[c], 1);
                            if (nxt_br != INF && nxt_br != nxt_onp && nxt_br == pp->vs[fpi].data[hi[fpi]]) {
                                hi[pi] = x + 1;
                                del_hm(&hd, pp->vs[pi].data[0]);
                                add_hm(&hd, pp->vs[pi].data[hi[pi]], pi);
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
                    } d++; prog(d, D, "pointing");
                } d++;
                while (!is_empty_q(&q)) {
                    u64 j = deq(&q);
                    for (size_t x = hi[j]; x < pp->vs[j].size; x++) {
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
                        } d++; prog(d, D, "pointing");
                    } d++;
                }
                free_v(&nc);
                free_v(pis);
                free(pis);
                break;
            }
            if (pis) { free_v(pis); free(pis); }
        }
        if (!fndc) break;
    }
    fin("pointing");
    l = nf(ino, pp->size);
    if (l == 0) printf("All paths are pointed\n");
    else        printf("%ld path(s) are not pointed\n", l);

    printf("# of (cycles, root paths, non-root paths) = (%ld, %d, %d)\n", cc->size, HASH_COUNT(rp), HASH_COUNT(hd));

    // take difference of r.arr
    if (pp->size > 1) {
        for (u64 i = pp->size - 1; i > 0; --i) {
            r.arr[i] -= r.arr[i - 1];
        }
    } else if (pp->size == 0) {
        free(r.arr);
        r.arr = NULL;
    }

    Strbld sb; init_strbld(&sb);
    for (size_t i = 0; i < o.size; i++) {
        if (o.data[i]) { // path
            u64 id = oi.data[i];
            V *cur = &pp->vs[id];
            u64 si = hi[id];
            if (find_hs(rp, id)) { // root path
                char s[k];
                dec(ka[cur->data[si]] >> 2, k - 1, s);
                apnd_strbld(&sb, ",");
                apnd_strbld(&sb, s);
            }
            for (size_t j = si; j < cur->size; j++) {
                apnd_strbld(&sb, dec_base(ka[cur->data[j]] % 4));
            }
        } else { // cycle
            V *cur = &cc->vs[oi.data[i]];
            apnd_strbld(&sb, "*");
            for (size_t j = 0; j < cur->size; j++) {
                apnd_strbld(&sb, dec_base(ka[cur->data[j]] % 4));
            }
        }
        apnd_strbld(&sb, ","); // delim
    }

    r.str = sb.str;
    size_t fl = strlen(r.str);
    if (fl > 0 && r.str[fl - 1] == ',') {
        r.str[fl - 1] = '\0';
    }
     if (fl > 1 && r.str[fl - 2] == ',' && r.str[fl - 1] == '\0') {
        r.str[fl - 2] = '\0';
    }
    free(hi); free_hs(&rp); free_hm(&hd); free_vb(&o); free_v(&oi);
    return r;
}

char* subt(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *vis, const u64 *hi, u64 rpi) {
    Strbld sb; init_strbld(&sb);
    FbSt s; init_fbst(&s);

    vis[rpi] = true;
    push_fbst(&s, rpi, hi[rpi], 0, false);
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
            push_fbst(&s, cpi, hi[cpi], 0, false);
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

//---Balanced parentheses---
Rep bp(Hm *km, u64 *ka, VV *cc, VV *pp, int k) {
    size_t Z = cc->size + pp->size, z = 0;
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
    u64 *hi = (u64*)calloc(pp->size, sizeof(u64));
    if (hi == NULL) {
        fprintf(stderr, "Error: calloc failed for hi\n");
        exit(EXIT_FAILURE);
    }

    Hs *s, *tmp;
    HASH_ITER(hh, rp, s, tmp) { // dfs from root paths
        char t[k]; // buffer for dec of head
        dec(ka[pp->vs[s->key].data[0]] >> 2, k - 1, t);
        char* ss = subt(km, ka, hd, k, pp, vis, hi, s->key);
        apnd_strbld(&sb, t);
        apnd_strbld(&sb, ss);
        apnd_strbld(&sb, ",");
        free(ss);
        z++; prog(z, Z, "embedding");
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
                char* ss = subt(km, ka, hd, k, pp, vis, hi, ni);
                apnd_strbld(&sb, "(");
                apnd_strbld(&sb, ss);
                apnd_strbld(&sb, ")");
                free(ss);
            }
        }
        apnd_strbld(&sb, ",");
        z++; prog(z, Z, "embedding");
    }

    u64 l;
    while ((l = nf(vis, pp->size))) {
        printf("\r%ld paths remaining\n", l);
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
                                hi[pi] = x + 1; // path pi's head is at index x + 1
                                del_hm(&hd, pp->vs[pi].data[0]);
                                add_hm(&hd, pp->vs[pi].data[x + 1], pi);
                                goto ej;
                            }
                        }
                    }
                    ej:;
                }
                push_backv(cc, nc);
                for (size_t j = 0; j < nc.size; j++) {
                    u64 cur = nc.data[j];
                    apnd_strbld(&sb, dec_base(ka[cur] % 4));
                    for (int c = 0; c < 4; c++) {
                        u64 nxt = step(km, ka, k, cur, B[c], 1);
                        if (nxt == INF) continue;
                        u64 ni = find_hm(hd, nxt);
                        if (ni == INF) continue;
                        if (vis[ni]) continue;
                        char* ss = subt(km, ka, hd, k, pp, vis, hi, ni);
                        apnd_strbld(&sb, "(");
                        apnd_strbld(&sb, ss);
                        apnd_strbld(&sb, ")");
                        free(ss);
                    }
                }
                apnd_strbld(&sb, ",");
                free_v(&nc);
                free_v(pis);
                free(pis);
                z++; prog(z, Z, "embedding");
                break;
            }
            if (pis) { free_v(pis); free(pis); }
        }
        if (!fndc) break;
    }
    fin("embedding");
    l = nf(vis, pp->size);
    if (l == 0) printf("All paths are embedded\n");
    else        printf("%ld path(s) are not embedded\n", l);

    printf("# of (cycles, root paths, non-root paths) = (%ld, %d, %d)\n", cc->size, HASH_COUNT(rp), HASH_COUNT(hd));

    if (sb.len > 0 && sb.str[sb.len - 1] == ',') {
        sb.len--;
        sb.str[sb.len] = '\0';
    }
    if (sb.len > 0 && sb.str[sb.len - 1] == ',') {
        sb.len--;
        sb.str[sb.len] = '\0';
    }
    free(hi); free_hs(&rp); free_hm(&hd);
    r.str = sb.str;
    return r;
}

static void subt_bgdfs(u64 tu, int c, Strbld *s, Hm *km, u64 *ka, int k, bool* vis) {
    FcSt t; init_fcst(&t);
    if (vis[tu]) return;
    
    vis[tu] = true;
    push_fcst(&t, tu, c, false, 0, 0);

    while(!is_empty_fcst(&t)) {
        Fc* f = t.top; // PEEK

        u64 h = ka[f->nid];
        if (f->side == 0) h = rc(h, k);

        if (!f->wrn) {
            if (f->next == NULL) { // root vertex
                char ts[k + 1];
                dec(h, k, ts);
                apnd_strbld(s, ts);
            } else {
                if (f->nc > 0) apnd_strbld(s, "(");
                apnd_strbld(s, dec_base(h % 4));
            }
            f->wrn = true;
        }

        if (f->b_idx > 0) {
            pop_fcst(&t);
            if (f->next != NULL && f->nc > 0) {
                 apnd_strbld(s, ")");
            }
            continue;
        }
        u64 chldn[4];
        int chldn_side[4];
        int nc = 0;

        for (uint8_t i = 0; i < 4; i++) { 
            for (int j = 1; j >= 0; j--) {
                u64 nxt = bstep(km, ka, k, f->nid, B[i], 1, f->side, j);
                if (nxt == INF) continue;
                if (vis[nxt]) continue;

                chldn[nc] = nxt;
                chldn_side[nc] = j;
                nc++;
            }
        }

        f->b_idx = 4;

        for (int i = 0; i < nc; i++) {
            vis[chldn[i]] = true;
            push_fcst(&t, chldn[i], chldn_side[i], false, 0, i);
        }
    }
}

Rep bgdfs(Hm *km, u64 *ka, int k) {
    Rep r; init_rep(&r);
    Strbld sb; init_strbld(&sb);
    
    u64 N = HASH_COUNT(km);
    bool* vis = (bool*)malloc(sizeof(bool) * N);
    if (vis == NULL) { fprintf(stderr, "Error: malloc failed for vis\n"); exit(EXIT_FAILURE); }
    for (u64 u = 0; u < N; u++) vis[u] = false;

    for (u64 u = 0; u < N; u++) {
        prog(u, N, "greedy BP dfs");
        if (vis[u]) continue;
        subt_bgdfs(u, 1, &sb, km, ka, k, vis);        
        apnd_strbld(&sb, ",");
    }
    
    fin("greedy BP dfs");
    free(vis);

    if (sb.len > 0 && sb.str[sb.len - 1] == ',') {
        sb.len--;
        sb.str[sb.len] = '\0';
    }
    if (sb.len > 0 && sb.str[sb.len - 1] == ',') {
        sb.len--;
        sb.str[sb.len] = '\0';
    }
    
    r.str = sb.str;
    return r;
}

static char* subt_bbp(u64 pid, bool forward, int offset, Hm *km, u64 *ka, Hm *hd, Hm *tl, int k, VV *pp, VVb *ppb, bool *vis, u64 *nvis) {
    prog(*nvis, pp->size, "embedding paths");
    Strbld sb; init_strbld(&sb);
    V *p = &pp->vs[pid];
    Vb *pb = &ppb->vs[pid];

    if (forward) {
        for (size_t j = offset; j < p->size; j++) {
            u64 cur = p->data[j];
            bool curside = pb->data[j];
            u64 h = ka[cur];
            if (curside == false) h = rc(h, k);
            apnd_strbld(&sb, dec_base(h % 4));
            u64 nxt_on_path = (j + 1 < p->size) ? p->data[j + 1] : INF;

            for (uint8_t i = 0; i < 4; i++) {
                for (int s = 1; s >= 0; s--) {
                    u64 nxt = bstep(km, ka, k, cur, B[i], 1, curside, s);
                    if (nxt == INF || nxt == nxt_on_path) continue;
                    u64 h = ka[nxt];
                    if (s == 0) h = rc(h, k);
                    u64 ni = find_hm(hd, h);
                    bool is_head = true;
                    if (ni == INF) {
                        ni = find_hm(tl, h);
                        is_head = false;
                    }
                    if (ni == INF) continue;
                    if (vis[ni]) continue;
                    vis[ni] = true; (*nvis)++;

                    char* ss = subt_bbp(ni, is_head, 0, km, ka, hd, tl, k, pp, ppb, vis, nvis);
                    apnd_strbld(&sb, "(");
                    apnd_strbld(&sb, ss);
                    apnd_strbld(&sb, ")");
                    free(ss);
                }
            }
        }
    } else {
        for (int64_t j = p->size - 1; j >= (int64_t)offset; j--) {
            u64 cur = p->data[j];
            bool curside = !pb->data[j];
            u64 h = ka[cur];
            if (curside == false) h = rc(h, k);
            apnd_strbld(&sb, dec_base(h % 4));
            u64 nxt_on_path = (j - 1 >= 0) ? p->data[j - 1] : INF;
            
            for (uint8_t i = 0; i < 4; i++) {
                for (int s = 1; s >= 0; s--) {
                    u64 nxt = bstep(km, ka, k, cur, B[i], 1, curside, s);
                    if (nxt == INF || nxt == nxt_on_path) continue;
                    u64 h = ka[nxt];
                    if (s == 0) h = rc(h, k);
                    u64 ni = find_hm(hd, h);
                    bool is_head = true;
                    if (ni == INF) {
                        ni = find_hm(tl, h);
                        is_head = false;
                    }
                    if (ni == INF) continue;
                    if (vis[ni]) continue;
                    vis[ni] = true; (*nvis)++;

                    char* ss = subt_bbp(ni, is_head, 0, km, ka, hd, tl, k, pp, ppb, vis, nvis);
                    apnd_strbld(&sb, "(");
                    apnd_strbld(&sb, ss);
                    apnd_strbld(&sb, ")");
                    free(ss);
                }
            }
        }
    }
    return sb.str;
}

Rep bbp(Hm *km, u64 *ka, VV *cc, VV *pp, VVb *ccb, VVb *ppb, int k) {
    Rep r; init_rep(&r);

    Strbld sbp; init_strbld(&sbp);
    Strbld sbc; init_strbld(&sbc);

    u64 Np = pp->size, nvis = 0; 
    bool *vis = (bool*)malloc(sizeof(bool) * Np);
    if (vis == NULL) { fprintf(stderr, "Error: malloc failed for vis\n"); exit(EXIT_FAILURE); }
    for (u64 i = 0; i < Np; i++) vis[i] = false;

    Hm *hd = NULL;
    Hm *tl = NULL;

    for (u64 i = 0; i < Np; i++) {
        V *p = &pp->vs[i];
        Vb *pb = &ppb->vs[i];
        if (p->size == 0) continue;

        u64 h = ka[p->data[0]];
        if (pb->data[0] == false) h = rc(h, k);
        add_hm(&hd, h, i); // raw hash of head considering the strand
        h = ka[p->data[p->size - 1]];
        if (pb->data[p->size - 1] == true) h = rc(h, k);
        add_hm(&tl, h, i); // raw hash of tail considering the strand
    }

    for (size_t i = 0; i < cc->size; i++) {
        V *c = &cc->vs[i];
        Vb *cb = &ccb->vs[i];
        if (c->size == 0) continue;

        for (size_t j = 0; j < c->size; j++) {
            u64 cur = c->data[j];
            bool curside = cb->data[j];
            u64 h = ka[cur];
            if (curside == false) h = rc(h, k);
            apnd_strbld(&sbc, dec_base(h % 4));
            u64 nxt_on_cyc = c->data[(j + 1) % c->size];

            for (uint8_t b = 0; b < 4; b++) {
                for (int toc = 1; toc >= 0; toc--) {
                    u64 nxt = bstep(km, ka, k, cur, B[b], 1, curside, toc);
                    if (nxt == INF || nxt == nxt_on_cyc) continue;
                    u64 h = ka[nxt];
                    if (toc == 0) h = rc(h, k);
                    u64 ni = find_hm(hd, h);
                    bool is_head = true;
                    if (ni == INF) {
                        ni = find_hm(tl, h);
                        is_head = false;
                    }
                    if (ni == INF) continue;
                    if (vis[ni]) continue;
                    vis[ni] = true; nvis++;

                    char* ss = subt_bbp(ni, is_head, 0, km, ka, hd, tl, k, pp, ppb, vis, &nvis);
                    apnd_strbld(&sbc, "(");
                    apnd_strbld(&sbc, ss);
                    apnd_strbld(&sbc, ")");
                    free(ss);
                }
            }
        }
        apnd_strbld(&sbc, ",");
    }

    for (u64 i = 0; i < Np; i++) {
        if (vis[i]) continue;
        vis[i] = true; nvis++;

        V *p = &pp->vs[i];
        Vb *pb = &ppb->vs[i];
        if (p->size == 0) continue;
        u64 head = p->data[0];
        bool headside = pb->data[0];
        u64 h = ka[head];
        if (headside == false) h = rc(h, k);
        char ts[k + 1];
        dec(h, k, ts);
        apnd_strbld(&sbp, ts);

        char* ss = subt_bbp(i, true, 1, km, ka, hd, tl, k, pp, ppb, vis, &nvis);
        apnd_strbld(&sbp, ss);
        free(ss);
        apnd_strbld(&sbp, ",");
    }
    fin("embedding paths");

    apnd_strbld(&sbp, ",");
    apnd_strbld(&sbp, sbc.str);
    if (sbp.str[sbp.len - 1] == ',') {
        sbp.len--;
        sbp.str[sbp.len] = '\0';
    }
    if (sbp.str[sbp.len - 1] == ',') {
        sbp.len--;
        sbp.str[sbp.len] = '\0';
    }
    free(sbc.str);

    free(vis);
    free_hm(&hd);
    free_hm(&tl);
    r.str = sbp.str;
    r.arr = NULL;
    return r;
}
