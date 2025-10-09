#include "rep_utils.h"
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

V* findc(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, u64 *in_pord, u64 *is_leaf, u64 *seen, u64 start) {
    V* nc = malloc(sizeof(V));
    if (nc == NULL) return NULL;
    init_v(nc);

    FSt s;
    init_fst(&s);

    seen[start] = 1;
    push_fst(&s, start, 0, 0);

    while (!is_empty_fst(&s)) {
        if (s.top->idx >= (u64)pp->vs[s.top->cur].size) {
            is_leaf[s.top->cur] = 1;
            pop_fst(&s);
            if (nc->size > 0) pop_back(nc);
            continue;
        }

        if (s.top->idx == 0 && s.top->b_idx == 0) push_back(nc, s.top->cur);
        
        if (s.top->b_idx < 4) {
            char c = B[s.top->b_idx];
            s.top->b_idx++;
            
            u64 cur_node = pp->vs[s.top->cur].data[s.top->idx];
            u64 next_node = step(km, ka, k, cur_node, c, 1);
            if (next_node == INF) continue;

            u64 next_pid = find_hm(hd, next_node);

            if (next_pid == INF || in_pord[next_pid] || seen[next_pid] || is_leaf[next_pid]) {
                continue;
            }

            if (next_pid == start) { // cycle found
                return nc;
            }
            
            seen[next_pid] = 1;
            push_fst(&s, next_pid, 0, 0);
        } else {
            s.top->idx++;
            s.top->b_idx = 0;
        }
    }
    free_v(nc);
    free(nc);
    return NULL; // no cycle found
}

u64 nt(bool *v, u64 l) {
    u64 n = 0;
    for (u64 i = 0; i < l; i++) {
        if (!v[i]) continue;
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

    while (nt(ino, pp->size)) {
        // from here
    }
    return r;
}

Rep bp(Hm *km, u64 *ka, VV *cc, VV *pp, int k) {
    Rep r; init_rep(&r);
    return r;
}