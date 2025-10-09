#ifndef REP_UTILS_H
#define REP_UTILS_H

#include "ds.h"
#include "utils.h"

typedef struct F {
    u64 cur;
    u64 idx;
    int b_idx;
    struct F *next;
} F;

typedef struct {
    F *top;
} FSt;

void init_fst(FSt *s);
int is_empty_fst(FSt *s);
void push_fst(FSt *s, u64 cur, u64 idx, int b_idx);
void pop_fst(FSt *s);

typedef struct {
    char *str;
    u64 *arr;
} Rep;

void init_rep(Rep *r);
void free_rep(Rep *r);

Rep flat(u64 *ka, VV *cc, VV *pp, int k);

V* findc(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, u64 *in_pord, u64 *is_leaf, u64 *visited, u64 start);
u64 nt(bool *v, size_t l);

Rep ptr(Hm *km, u64 *ka, VV *cc, VV *pp, int k);
Rep bp(Hm *km, u64 *ka, VV *cc, VV *pp, int k);

#endif // REP_UTILS_H