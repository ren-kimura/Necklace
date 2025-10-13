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

typedef struct Fb {
    u64 cur;
    u64 idx;
    int b_idx;
    bool wrn;
    struct Fb *next;
} Fb;

typedef struct {
    Fb *top;
} FbSt;

void init_fbst(FbSt *s);
int is_empty_fbst(FbSt *s);
void push_fbst(FbSt *s, u64 cur, u64 idx, int b_idx, bool wrn);
void pop_fbst(FbSt *s);

typedef struct {
    char* str;
    size_t len;
    size_t cap;
} Strbld;

void init_strbld(Strbld *s);
void apnd_strbld(Strbld *s, const char* t);

typedef struct {
    char *str;
    u64 *arr;
} Rep;

void init_rep(Rep *r);
void free_rep(Rep *r);

Rep flat(u64 *ka, VV *cc, VV *pp, int k);

V* findc(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *ino, bool *vis, u64 from);
u64 nf(bool *v, size_t l);

Rep ptr(Hm *km, u64 *ka, VV *cc, VV *pp, int k);

char* subt(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *vis, const u64 *hi, u64 rpi);
Rep bp(Hm *km, u64 *ka, VV *cc, VV *pp, int k);

#endif // REP_UTILS_H