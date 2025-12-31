#include "eutils.h"
#include <string.h>

static void reverse_v(V *v) {
    if (v->size < 2) return;
    size_t s = 0, e = v->size - 1;
    while (s < e) {
        u64 tmp = v->data[s];
        v->data[s] = v->data[e];
        v->data[e] = tmp;
        s++;
        e--;
    }
}

Node* build(Hm *km, int k, u64 N) {
    Node *g = NULL;
    u64 m = (1ULL << (2 * (k - 1))) - 1;
    u64 n = 0; // for prog

    Hm *s, *tmp;
    HASH_ITER(hh, km, s, tmp) {
        prog(++n, N, "building graph");
        u64 h = s->key;
        u64 pre = h >> 2;
        u64 suf = h & m;
        uint8_t f = (h >> (2 * (k - 1))) & 3;
        uint8_t l = h & 3;

        Node *np;
        HASH_FIND(hh, g, &pre, sizeof(u64), np);
        if (!np) {
            np = malloc(sizeof(Node));
            np->key = pre;
            np->a = 0;
            np->o = 0;
            init_v(&np->dm);
            HASH_ADD(hh, g, key, sizeof(u64), np);
        }

        Node *ns;
        HASH_FIND(hh, g, &suf, sizeof(u64), ns);
        if (!ns) {
            ns = malloc(sizeof(Node));
            ns->key = suf;
            ns->a = 0;
            ns->o = 0;
            init_v(&ns->dm);
            HASH_ADD(hh, g, key, sizeof(u64), ns);
        }
        
        np->a |= (1 << (l + 4)); // record outneighbor
        ns->a |= (1 << f); // record inneighbor
        np->o++;
    }
    fin("building graph");

    printf("%d k-1-mers\n", HASH_COUNT(g));
    return g;
}

void balance(Node **g) {
    V src, snk;
    init_v(&src); init_v(&snk);
    size_t nn = 0, NN = HASH_COUNT(*g);

    Node *s, *tmp;
    HASH_ITER(hh, *g, s, tmp) {
        prog(++nn, NN, "counting imbalances");
        uint8_t i = 0, o = 0;
        for (uint8_t c = 0; c < 4; c++) {
            if (s->a & (1 << c)) i++;
            if (s->a & (1 << (c + 4))) o++;
        }

        if (o > i) {
            for (uint8_t d = 0; d < o - i; d++) push_back(&src, s->key);
        } else if (o < i) {
            for (uint8_t d = 0; d < i - o; d++) push_back(&snk, s->key);
        }
    }
    fin("counting imbalances");

    if (src.size != snk.size) {
        fprintf(stderr, "Error: different number of sources and sinks\n");
        exit(EXIT_FAILURE);
    }

    printf("Adding %ld dummy edges to balance the graph\n", src.size);
    for (size_t i = 0; i < src.size; i++) {
        prog(i + 1, src.size, "adding dummy edges");
        u64 ksnk = snk.data[i];
        u64 ksrc = src.data[i];
        
        Node *nsnk;
        HASH_FIND(hh, *g, &ksnk, sizeof(u64), nsnk);
        if (nsnk) {
            push_back(&nsnk->dm, ksrc);
            nsnk->o++;
        }
    }
    free_v(&src); free_v(&snk);
    fin("adding dummy edges");
}

void etigs(Node **g, VV *tt, int k) {
    size_t ne = 0;
    Node *s, *tmp;
    HASH_ITER(hh, *g, s, tmp) {
        ne += s->o;
    }

    size_t nue = 0;
    while (nue < ne) {
        Node *ns = NULL;
        HASH_ITER(hh, *g, s, tmp) {
            if (s->o > 0) {
                ns = s;
                break;
            }
        }
        if (!ns) break; // visited all edges

        // Hierholzer

        St p;
        init_st(&p);
        push(&p, ns->key);
        V t; init_v(&t);

        while (!is_empty_st(&p)) {
            u64 ku = p.top->data; // peek
            Node *nu;
            HASH_FIND(hh, *g, &ku, sizeof(u64), nu);
            if (nu && nu->o > 0) {
                nu->o--;
                nue++;
                prog(nue, ne, "computing Eulertigs");
                if (nu->dm.size > 0) {
                    u64 kv = pop_back(&nu->dm);
                    push(&p, kv);
                } else {
                    for (uint8_t b = 0; b < 4; b++) {
                        if (nu->a & (1 << (b + 4))) {
                            u64 m = (1ULL << (2 * (k - 2))) - 1; // take k-2-mer
                            u64 kv = ((ku & m) << 2) | b;
                            push(&p, kv);
                            nu->a &= ~(1 << (b + 4));
                            break;
                        }
                    }
                }
            } else {
                ku = pop(&p);
                push_back(&t, ku);
            }
        }

        reverse_v(&t);
        push_backv(tt, t);
        free_v(&t);
    }
    if (ne > 0) fin("computing Eulertigs");
    printf("%ld tours\n", tt->size);
}

void free_ss(char **ss, size_t ns) {
    if (!ss) return;
    for (size_t i = 0; i < ns; i++) {
        free(ss[i]);
    }
    free(ss);
}

char** spell(VV *tt, int k, size_t *ns) {
    *ns = 0;
    size_t scap = 8;
    char **ss = malloc(scap * sizeof(char*));
    if (!ss) return NULL;

    size_t bufcap = 256;
    char *buf = malloc(bufcap);
    if (!buf) { free(ss); return NULL; }

    for (size_t i = 0; i < tt->size; i++) {
        V *t = &tt->vs[i];
        if (tt->size == 0) continue;
        dec(t->data[0], k - 1, buf);
        size_t buflen = strlen(buf);

        for (size_t j = 0; j < t->size - 1; j++) {
            u64 ku = t->data[j];
            u64 kv = t->data[j + 1];
            u64 musuf = (1ULL << (2 * (k - 2))) - 1;
            u64 usuf = ku & musuf;
            u64 vpre = kv >> 2;

            if (usuf == vpre) {
                if (buflen + 1 + 1 > bufcap) {
                    bufcap *= 2;
                    buf = realloc(buf, bufcap);
                    if (buf == NULL) {
                        fprintf(stderr, "Error: realloc failed for buffer\n");
                        exit(EXIT_FAILURE);
                    }
                }
                buf[buflen++] = B[kv & 3];
                buf[buflen] = '\0';
            } else {
                if (*ns >= scap) {
                    scap *= 2;
                    ss = realloc(ss, scap * sizeof(char*));
                    if (ss == NULL) {
                        fprintf(stderr, "Error: realloc failed for contigs\n");
                        exit(EXIT_FAILURE);
                    }
                }
                ss[*ns] = strdup(buf);
                (*ns)++;

                dec(kv, k - 1, buf);
                buflen = strlen(buf);
            }
        }
        
        if (*ns >= scap) {
            scap *= 2;
            ss = realloc(ss, scap * sizeof(char*));
            if (ss == NULL) {
                fprintf(stderr, "Error: realloc failed for contigs\n");
                exit(EXIT_FAILURE);
            }
        }
        ss[*ns] = strdup(buf);
        (*ns)++;
    }
        
    free(buf);
    return ss;
}

void tt_to_cc_and_pp(VV *tt, Hm *km, VV *cc, VV *pp) {
    for (size_t i = 0; i < tt->size; i++) {
        V *t = &tt->vs[i];
        if (t->size < 2) continue; // there is no k-mer in t

        V p; init_v(&p);
        bool has_dummy = false;

        for (size_t j = 0; j < t->size - 1; j++) {
            u64 h = ((t->data[j] << 2) | (t->data[j + 1] & 3)); // edge hash
            u64 id = find_hm(km, h);

            if (id != INF) {
                push_back(&p, id);
            } else {
                if (p.size > 0) {
                    push_backv(pp, p);
                    free_v(&p);
                    init_v(&p);
                }
                has_dummy = true;
            }
        }

        if (p.size > 0) {
            if (!has_dummy && t->data[0] == t->data[t->size - 1]) {
                push_backv(cc, p);
            } else {
                push_backv(pp, p);
            }
            free_v(&p);
        }
    }
}

void free_g(Node **g) {
    Node *s, *tmp;
    HASH_ITER(hh, *g, s, tmp) {
        free_v(&s->dm);
        HASH_DEL(*g, s);
        free(s);
    }
}