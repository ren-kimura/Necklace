#include "out.h"
#include "stat.h"
#include <string.h>

//---frame stack---
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

//---frame stack---
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

//---frame stack---
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

//---bool stack---
void init_blst(BlSt *s) {
    s->top = NULL;
}

int is_empty_blst(BlSt *s) {
    return (s->top == NULL);
}

void push_blst(BlSt *s, bool is_main) {
    Bl *nbl = (Bl*)malloc(sizeof(Bl));
    if (nbl == NULL) {
        fprintf(stderr, "Error: malloc failed for new Bool element\n");
        exit(EXIT_FAILURE);
    }
    nbl->is_main = is_main;
    nbl->next = s->top;
    s->top = nbl;
}

void pop_blst(BlSt *s) {
    if (is_empty_blst(s)) {
        fprintf(stderr, "Error: pop_blst called on an empty stack\n");
        return;
    }
    Bl *tmp = s->top;
    s->top = s->top->next;
    free(tmp);
}

// --- incoming and outgoing edges of nodes in rbp ---
typedef struct Link {
    u64 u_id;      // Node ID (Parent)
    V p_idxs;      // List of Path IDs starting from this node
    UT_hash_handle hh;
} Link;

void add_link(Link **links, u64 u_id, u64 p_idx) {
    Link *l;
    HASH_FIND(hh, *links, &u_id, sizeof(u64), l);
    if (!l) {
        l = (Link*)malloc(sizeof(Link));
        l->u_id = u_id;
        init_v(&l->p_idxs);
        HASH_ADD(hh, *links, u_id, sizeof(u64), l);
    }
    push_back(&l->p_idxs, p_idx);
}

void free_links(Link **links) {
    Link *curr, *tmp;
    HASH_ITER(hh, *links, curr, tmp) {
        free_v(&curr->p_idxs);
        HASH_DEL(*links, curr);
        free(curr);
    }
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

// spelling functions for FirstVar
char* pccover_to_cspss(u64 *ka, VV *cc, VV *pp, int k) {
    char* r;
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

    r = sb.str;
    size_t fl = strlen(r);
    if (fl > 0 && r[fl - 1] == ',') {
        r[fl - 1] = '\0';
    }
     if (fl > 1 && r[fl - 2] == ',' && r[fl - 1] == '\0') {
        r[fl - 2] = '\0';
    }
    return r;
}

char* bi_pccover_to_cspss(u64 *ka, VV *cc, VV *pp, VVb *ccb, VVb *ppb, int k) {
    char* r;
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

    r = sb.str;
    size_t fl = strlen(r);
    if (fl > 0 && r[fl - 1] == ',') {
        r[fl - 1] = '\0';
    }
     if (fl > 1 && r[fl - 2] == ',' && r[fl - 1] == '\0') {
        r[fl - 2] = '\0';
    }
    return r;
}

// spelling functions for SecondVar
V* find_new_cycle(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *ino, bool *vis, u64 from) {
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

u64 num_falses(bool *v, u64 l) {
    u64 n = 0;
    for (u64 i = 0; i < l; i++) {
        if (v[i]) continue; // skip if true
        n++;
    }
    return n;
}

char* find_subtree(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *vis, const u64 *hi, u64 rpi) {
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

char* necklace_cover(Hm *km, u64 *ka, VV *cc, VV *pp, int k) {
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
    char* r = NULL;
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
        char* ss = find_subtree(km, ka, hd, k, pp, vis, hi, s->key);
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
                char* ss = find_subtree(km, ka, hd, k, pp, vis, hi, ni);
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
    while ((l = num_falses(vis, pp->size))) {
        printf("\r%ld paths remaining\n", l);
        bool fndc = false;
        for (u64 i = 0; i < pp->size; i++) {
            if (vis[i]) continue;
            bool tvis[pp->size];
            for (size_t j = 0; j < pp->size; j++) tvis[j] = false;
            V* pis = find_new_cycle(km, ka, hd, k, pp, vis, tvis, i);
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
                        char* ss = find_subtree(km, ka, hd, k, pp, vis, hi, ni);
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
    l = num_falses(vis, pp->size);
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
    r = sb.str;
    return r;
}

// variation of necklace_cover: attach non-primitive paths to its arbitrary parent
// may yield more closed necklaces than necklace_cover
// is_head_full: if true, print the full k-mer of the first node. Otherwise print only the last char.
void subr_necklace_cover2(u64 pid, VV *pp, Link **links, bool *vis, Strbld *sb, u64 *ka, int k, bool is_head_full, u64 *parent, size_t *z, size_t Z) {
    if (vis[pid]) return;
    vis[pid] = true;
    if (parent[pid] != INF) {
        (*z)++; // attached a non-root open path
        prog(*z, Z, "embedding");
    }

    V *p = &pp->vs[pid];
    for (size_t i = 0; i < p->size; i++) {
        u64 u = p->data[i];

        // Output the node string
        if (i == 0 && is_head_full) {
            char s[k + 1];
            dec(ka[u], k, s);
            apnd_strbld(sb, s);
        } else {
            apnd_strbld(sb, dec_base(ka[u] % 4));
        }

        // Check if this node has attached child paths (branches)
        Link *l;
        HASH_FIND(hh, *links, &u, sizeof(u64), l);
        if (l) {
            for (size_t j = 0; j < l->p_idxs.size; j++) {
                u64 child_pid = l->p_idxs.data[j];
                if (!vis[child_pid]) { // Avoid infinite recursion if cycle exists
                    apnd_strbld(sb, "(");
                    // Children are attached, so they share k-1 overlap. 
                    // Their head is already represented by u's suffix. 
                    // But in standard notation, if we branch, we usually continue the sequence.
                    // Here, child starts with successor of u. We print it incrementally.
                    subr_necklace_cover2(child_pid, pp, links, vis, sb, ka, k, false, parent, z, Z);
                    apnd_strbld(sb, ")");
                }
            }
        }
    }
}

char* necklace_cover2(Hm *km, u64 *ka, VV *cc, VV *pp, int k) {
    char* r;
    Strbld sb; init_strbld(&sb);
    u64 Np = pp->size;
    u64 Nc = cc->size;
    
    // Setup Data Structures
    Link *links = NULL; // Map: Node ID -> List of Path IDs (children)
    u64 *parent = (u64*)malloc(Np * sizeof(u64));
    bool *vis = (bool*)calloc(Np, sizeof(bool));
    
    if (!parent || !vis) {
        fprintf(stderr, "Error: malloc failed in rbp\n");
        exit(EXIT_FAILURE);
    }
    for(size_t i = 0; i < Np; i++) { parent[i] = INF; }

    Hm *node2pid = NULL; // add Np to cycle id to distinguish paths and cycles
    for (u64 i = 0; i < Np; i++) {
        for (u64 j = 0; j < pp->vs[i].size; j++) {
            add_hm(&node2pid, pp->vs[i].data[j], i);
        }
    }
    for (u64 i = 0; i < Nc; i++) {
        for (u64 j = 0; j < cc->vs[i].size; j++) {
            add_hm(&node2pid, cc->vs[i].data[j], Np + i);
        }
    }

    // For each path, try to find *one* incoming edge from the graph
    size_t Z = 0;
    for (size_t i = 0; i < Np; i++) {
        if (pp->vs[i].size == 0) continue;
        u64 head = pp->vs[i].data[0];

        // backward search for a predecessor
        for (int j = 0; j < 4; j++) {
            u64 pred = step(km, ka, k, head, B[j], 0); // 0 = backward
            if (pred != INF) {
                u64 p_pid = find_hm(node2pid, pred);
                if (p_pid != INF) {
                    parent[i] = p_pid;
                    add_link(&links, pred, i);
                    Z++;
                    break;
                }
            }
        }
    }

    free_hm(&node2pid);

    size_t z = 0;
    printf("# of non-root open paths: %ld\n", Z);

    // from root open paths
    for (size_t i = 0; i < Np; i++) {
        if (parent[i] == INF && !vis[i]) {
            subr_necklace_cover2(i, pp, &links, vis, &sb, ka, k, true, parent, &z, Z);
            apnd_strbld(&sb, ",");
        }
    }

    apnd_strbld(&sb, ","); // separator between open and closed necklaces

    // from closed paths
    for (size_t i = 0; i < Nc; i++) {
        V *c = &cc->vs[i];
        if (c->size == 0) continue;
        
        for (size_t j = 0; j < c->size; j++) {
            u64 u = c->data[j];                 
            apnd_strbld(&sb, dec_base(ka[u] % 4));

            // Check for attached paths
            Link *l;
            HASH_FIND(hh, links, &u, sizeof(u64), l);
            if (l) {
                for (size_t x = 0; x < l->p_idxs.size; x++) {
                    u64 cpid = l->p_idxs.data[x];
                    if (!vis[cpid]) {
                        apnd_strbld(&sb, "(");
                        subr_necklace_cover2(cpid, pp, &links, vis, &sb, ka, k, false, parent, &z, Z);
                        apnd_strbld(&sb, ")");
                    }
                }
            }
        }
        apnd_strbld(&sb, ",");
    }

    // from remaining non-root open paths
    u64 *stack = (u64*)malloc(Np * sizeof(u64));
    bool *in_stack = (bool*)calloc(Np, sizeof(bool));

    for (u64 i = 0; i < Np; i++) {
        if (vis[i]) continue;

        u64 cur = i;
        size_t top = 0;
        while (cur != INF && !vis[cur] && !in_stack[cur]) {
            in_stack[cur] = true;
            stack[top++] = cur;
            cur = parent[cur];
        }
        if (cur != INF && in_stack[cur]) {
            size_t start_idx = 0;
            while (stack[start_idx] != cur) start_idx++;
            // arrange in forward order (stack[start_idx...top-1] is in reverse)
            u64 c_len = top - start_idx;
            u64 *c_pids = (u64*)malloc(c_len * sizeof(u64));
            for (size_t j = 0; j < c_len; j++) c_pids[j] = stack[top - 1 - j];

            // build closed necklace
            for (size_t j = 0; j < c_len; j++) {
                u64 pid = c_pids[j];
                u64 next_pid = c_pids[(j + 1) % c_len];
                vis[pid] = true;
                if (parent[pid] != INF) {
                    z++; // attached a non-root open path
                    prog(z, Z, "embedding");
                }

                V *p = &pp->vs[pid];

                bool opened = false;
                for (size_t x = 0; x < p->size; x++) {
                    u64 u = p->data[x];
                    apnd_strbld(&sb, dec_base(ka[u] % 4));

                    Link *l;
                    HASH_FIND(hh, links, &u, sizeof(u64), l);
                    bool to_next = false;
                    if (l) {
                        for (size_t b = 0; b < l->p_idxs.size; b++) {
                            u64 child_pid = l->p_idxs.data[b];
                            if (child_pid == next_pid) {
                                to_next = true;
                            } else if (!vis[child_pid]) {
                                apnd_strbld(&sb, "(");
                                subr_necklace_cover2(child_pid, pp, &links, vis, &sb, ka, k, false, parent, &z, Z);
                                apnd_strbld(&sb, ")");
                            }
                        }
                    }
                    if (to_next) {
                        apnd_strbld(&sb, "(");
                        opened = true;
                    }
                }
                if (opened) {
                    apnd_strbld(&sb, ")");
                }
            }
            apnd_strbld(&sb, ",");
            free(c_pids);
        }
        for (size_t j = 0; j < top; j++) in_stack[stack[j]] = false;
    }
    fin("embedding");

    // Cleanup string
    if (sb.len > 0 && sb.str[sb.len - 1] == ',') {
        sb.len--;
        sb.str[sb.len] = '\0';
    }

    free(parent); free(vis); free_links(&links);
    free(stack); free(in_stack);

    r = sb.str;
    return r;
}

// greedy baseline: iterate DFS from any unvisited vertex to form necklaces
static void subr_greedy_baseline(u64 tu, Strbld *s, Hm *km, u64 *ka, int k, bool* vis) {
    FcSt t; init_fcst(&t);
    if (vis[tu]) return;

    vis[tu] = true;
    push_fcst(&t, tu, 0, false, 0, 0);

    while (!is_empty_fcst(&t)) {
        Fc* f = t.top; // PEEK

        u64 h = ka[f->nid];
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
        int nc = 0;

        for (uint8_t i = 0; i < 4; i++) {
            u64 nxt = step(km, ka, k, f->nid, B[i], 1);
            if (nxt == INF) continue;
            if (vis[nxt]) continue;

            chldn[nc] = nxt;
            nc++;
        }

        f->b_idx = 4;

        for (int i = 0; i < nc; i++) {
            vis[chldn[i]] = true;
            push_fcst(&t, chldn[i], 0, false, 0, i);
        }
    }
}

char* greedy_baseline(Hm *km, u64 *ka, int k) {
    char* r;
    Strbld sb; init_strbld(&sb);

    u64 N = HASH_COUNT(km);
    bool* vis = (bool*)calloc(N, sizeof(bool));
    if (!vis) { fprintf(stderr, "Error: calloc failed for vis\n"); exit(EXIT_FAILURE); }

    for (u64 u = 0; u < N; u++) {
        prog(u, N, "greedy BP dfs");
        if (vis[u]) continue;
        subr_greedy_baseline(u, &sb, km ,ka, k, vis);
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

    r = sb.str;
    return r;
}

static bool subr_greedy_baseline_close(u64 tu, Strbld *s, Hm *km, u64 *ka, int k, bool* vis) {
    vis[tu] = true;
    FcSt t; init_fcst(&t); BlSt b; init_blst(&b);
    push_fcst(&t, tu, 0, false, 0, 0);
    push_blst(&b, true);

    bool main_is_cycle = false;

    while (!is_empty_fcst(&t)) {
        Fc* f = t.top; 
        bool is_main = b.top->is_main;
        u64 h = ka[f->nid];

        if (!f->wrn) {
            if (f->next == NULL) {
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
            pop_blst(&b);
            if (f->next != NULL && f->nc > 0) apnd_strbld(s, ")");
            continue;
        }

        u64 chldn[4];
        int nc = 0;
        bool cycle_found_now = false;

        for (uint8_t i = 0; i < 4; i++) {
            u64 nxt = step(km, ka, k, f->nid, B[i], 1);
            if (nxt == INF) continue;

            if (is_main && nxt == tu) {
                main_is_cycle = true;
                cycle_found_now = true;
                break; 
            }

            if (vis[nxt]) continue;
            chldn[nc++] = nxt;
        }

        f->b_idx = 4;
        if (cycle_found_now) continue;

        for (int i = 0; i < nc; i++) {
            vis[chldn[i]] = true; 
            bool next_is_main = (is_main && i == 0);
            push_fcst(&t, chldn[i], 0, false, 0, i);
            push_blst(&b, next_is_main);
        }
    }

    return main_is_cycle;
}

char* greedy_baseline_close(Hm *km, u64 *ka, int k) {
    char* r;
    Strbld sb_p; init_strbld(&sb_p);
    Strbld sb_c; init_strbld(&sb_c);

    u64 N = HASH_COUNT(km);
    bool* vis = (bool*)calloc(N, sizeof(bool));
    if (!vis) { fprintf(stderr, "Error: calloc failed for vis\n"); exit(EXIT_FAILURE); }

    for (u64 u = 0; u < N; u++) {
        prog(u, N, "greedy BP dfs");
        if (vis[u]) continue;

        Strbld sb_tmp; init_strbld(&sb_tmp);
        bool is_cycle = subr_greedy_baseline_close(u, &sb_tmp, km, ka, k, vis);
        if (is_cycle) {
            apnd_strbld(&sb_c, sb_tmp.str + (k - 1));
        } else {
            apnd_strbld(&sb_p, sb_tmp.str);
        }
        apnd_strbld(&sb_p, ",");
        free(sb_tmp.str);
    }

    Strbld final_sb; init_strbld(&final_sb);
    apnd_strbld(&final_sb, sb_p.str);
    apnd_strbld(&final_sb, ",");
    apnd_strbld(&final_sb, sb_c.str);

    fin("greedy BP dfs");
    free(vis);
    free(sb_p.str); free(sb_c.str);

    if (final_sb.len > 0 && final_sb.str[final_sb.len - 1] == ',') {
        final_sb.len--;
        final_sb.str[final_sb.len] = '\0';
    }
    if (final_sb.len > 0 && final_sb.str[final_sb.len - 1] == ',') {
        final_sb.len--;
        final_sb.str[final_sb.len] = '\0';
    }

    r = final_sb.str;
    return r;
}

// full greedy: for every non-root vertex, choose arbitrary one of its parents to form necklaces
void subr_full_greedy(u64 u, Hm *child, bool *vis, Strbld *sb, Hm *km, u64 *ka, int k, bool is_head_full, u64 *parent, size_t *n, size_t N) {
    if (vis[u]) return;
    vis[u] = true;

    if (parent[u] != INF) {
        prog((*n)++, N, "finding necklace");
    }

    if (is_head_full) {
        char s[k + 1];
        dec(ka[u], k, s);
        apnd_strbld(sb, s);
    } else {
        apnd_strbld(sb, dec_base(ka[u] % 4));
    }

    u64 branches = find_hm(child, u);
    if (branches != INF) {
        for (int b = 0; b < 4; b++) {
            if (branches & (1 << b)) {
                u64 cid = step(km, ka, k, u, B[b], 1);
                if (cid != INF && !vis[cid]) {
                    apnd_strbld(sb, "(");
                    subr_full_greedy(cid, child, vis, sb, km, ka, k, false, parent, n, N);
                    apnd_strbld(sb, ")");
                }
            }
        }
    }
}

char* full_greedy(Hm *km, u64 *ka, int k) {
    size_t N = HASH_COUNT(km);
    char* r;
    Strbld sb; init_strbld(&sb);
    Hm *child = NULL; // map: node id -> integer 0~15 representing branches
    u64 *parent = (u64*)malloc(N * sizeof(u64));
    bool *vis = (bool*)calloc(N, sizeof(bool));

    if (!parent || !vis) exit(EXIT_FAILURE);
    for (size_t i = 0; i < N; i++) { parent[i] = INF; }

    for (size_t i = 0; i < N; i++) {
        for (int j = 0; j < 4; j++) {
            u64 pred = step(km, ka, k, i, B[j], 0); // 0: backward
            if (pred != INF) {
                parent[i] = pred;
                update_branch_hm(&child, pred, j);
                break;
            }
        }
    }

    size_t n = 0;

    // from root nodes
    for (size_t i = 0; i < N; i++) {
        if (parent[i] == INF && !vis[i]) {
            subr_full_greedy(i, child, vis, &sb, km, ka, k, true, parent, &n, N);
            apnd_strbld(&sb, ",");
        }
    }
    apnd_strbld(&sb, ",");

    // detect cycle and find necklace
    for (size_t i = 0; i < N; i++) {
        if (vis[i]) continue;
        St s; init_st(&s); // for backtracking
        Hs *in_s = NULL; // for backtracking
        u64 ti = i;

        while (ti != INF && !vis[ti] && !find_hs(in_s, ti)) {
            add_hs(&in_s, ti);
            push(&s, ti);
            ti = parent[ti];
        }

        if (ti != INF && find_hs(in_s, ti)) {
            V c; init_v(&c);
            u64 cur;

            V tmp; init_v(&tmp);
            while (!is_empty_st(&s)) {
                cur = pop(&s);
                push_back(&tmp, cur);
            }
            
            for (int j = tmp.size - 1; j >= 0; j--) {
                push_back(&c, tmp.data[j]);
            }
            free_v(&tmp);

            // build closed necklace
            for (size_t j = 0; j < c.size; j++) {
                u64 id = c.data[j];
                u64 nid = c.data[(j + 1) % c.size];
                vis[id] = true;
                apnd_strbld(&sb, dec_base(ka[id] % 4));
                prog(n++, N, "finding necklace");

                u64 branches = find_hm(child, id);
                if (branches != INF) {
                    for (int b = 0; b < 4; b++) {
                        if (branches & (1 << b)) {
                            u64 cid = step(km, ka, k, id, B[b], 1);
                            if (cid != nid && cid != INF && !vis[cid]) {
                                apnd_strbld(&sb, "(");
                                subr_full_greedy(cid, child, vis, &sb, km, ka, k, false, parent, &n, N);
                                apnd_strbld(&sb, ")");
                            }
                        }
                    }
                }
            }
            apnd_strbld(&sb, ",");
            free_v(&c);
        }
        free_hs(&in_s);
    }
    fin("finding necklace");

    if (sb.len > 0 && sb.str[sb.len - 1] == ',') {
        sb.len--;
        sb.str[sb.len] = '\0';
    }

    free(parent); free(vis); free_hm(&child);
    
    r = sb.str;
    return r;
}

// greedy baseline (bi-directed)
static void subr_bi_greedy_baseline(u64 tu, int c, Strbld *s, Hm *km, u64 *ka, int k, bool* vis) {
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

char* bi_greedy_baseline(Hm *km, u64 *ka, int k) {
    char* r;
    Strbld sb; init_strbld(&sb);
    
    u64 N = HASH_COUNT(km);
    bool* vis = (bool*)malloc(sizeof(bool) * N);
    if (vis == NULL) { fprintf(stderr, "Error: malloc failed for vis\n"); exit(EXIT_FAILURE); }
    for (u64 u = 0; u < N; u++) vis[u] = false;

    for (u64 u = 0; u < N; u++) {
        prog(u, N, "greedy BP dfs");
        if (vis[u]) continue;
        subr_bi_greedy_baseline(u, 1, &sb, km, ka, k, vis);        
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
    
    r = sb.str;
    return r;
}

// necklace_cover (bi-directed)
static void subr_bi_necklace_cover(Strbld *sb, u64 pid, bool forward, int offset, Hm *km, u64 *ka, Hm *hd, Hm *tl, int k, VV *pp, VVb *ppb, bool *vis, u64 *nvis) {
    prog(*nvis, pp->size, "embedding paths");
    V *p = &pp->vs[pid];
    Vb *pb = &ppb->vs[pid];

    if (forward) {
        for (size_t j = offset; j < p->size; j++) {
            u64 cur = p->data[j];
            bool curside = pb->data[j];
            u64 h = ka[cur];
            if (curside == false) h = rc(h, k);
            apnd_strbld(sb, dec_base(h % 4));
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

                    apnd_strbld(sb, "(");
                    subr_bi_necklace_cover(sb, ni, is_head, 0, km, ka, hd, tl, k, pp, ppb, vis, nvis);
                    apnd_strbld(sb, ")");
                }
            }
        }
    } else {
        for (int64_t j = p->size - 1; j >= (int64_t)offset; j--) {
            u64 cur = p->data[j];
            bool curside = !pb->data[j];
            u64 h = ka[cur];
            if (curside == false) h = rc(h, k);
            apnd_strbld(sb, dec_base(h % 4));
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

                    apnd_strbld(sb, "(");
                    subr_bi_necklace_cover(sb, ni, is_head, 0, km, ka, hd, tl, k, pp, ppb, vis, nvis);
                    apnd_strbld(sb, ")");
                }
            }
        }
    }
}

char* bi_necklace_cover(Hm *km, u64 *ka, VV *cc, VV *pp, VVb *ccb, VVb *ppb, int k) {
    char* r = NULL;

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

                    apnd_strbld(&sbc, "(");
                    subr_bi_necklace_cover(&sbc, ni, is_head, 0, km, ka, hd, tl, k, pp, ppb, vis, &nvis);
                    apnd_strbld(&sbc, ")");
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

        subr_bi_necklace_cover(&sbp, i, true, 1, km, ka, hd, tl, k, pp, ppb, vis, &nvis);
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
    r = sbp.str;
    return r;
}
