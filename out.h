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

typedef struct Fc {
    u64 nid;
    int side;
    bool wrn;
    int b_idx;
    int nc;
    struct Fc *next;
} Fc;

typedef struct {
    Fc *top;
} FcSt;

void init_fcst(FcSt *s);
int is_empty_fcst(FcSt *s);
void push_fcst(FcSt *s, u64 nid, int side, bool wrn, int b_idx, int nc);
void pop_fcst(FcSt *s);

typedef struct Bl {
    bool is_main;
    struct Bl *next;
} Bl;

typedef struct {
    Bl *top;
} BlSt;

void init_blst(BlSt *s);
int is_empty_blst(BlSt *s);
void push_blst(BlSt *s, bool is_main);
void pop_blst(BlSt *s);

typedef struct {
    char* str;
    size_t len;
    size_t cap;
} Strbld;

void init_strbld(Strbld *s);
void apnd_strbld(Strbld *s, const char* t);

char* pccover_to_cspss(u64 *ka, VV *cc, VV *pp, int k);
char* bi_pccover_to_cspss(u64 *ka, VV *cc, VV *pp, VVb *ccb, VVb *ppb, int k);

V* find_new_cycle(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *ino, bool *vis, u64 from);
u64 num_falses(bool *v, size_t l);

char* find_subtree(Hm *km, u64 *ka, Hm *hd, int k, VV *pp, bool *vis, const u64 *hi, u64 rpi);
char* necklace_cover(Hm *km, u64 *ka, VV *cc, VV *pp, int k);
char* necklace_cover2(Hm *km, u64 *ka, VV *cc, VV *pp, int k);

char* greedy_baseline(Hm *km, u64 *ka, int k);
char* greedy_baseline_close(Hm *km, u64 *ka, int k);

char* full_greedy(Hm *km, u64 *ka, int k);

char* bi_greedy_baseline(Hm *km, u64 *ka, int k);

char* bi_necklace_cover(Hm *km, u64 *ka, VV *cc, VV *pp, VVb *ccb, VVb *ppb, int k);

#endif // REP_UTILS_H
