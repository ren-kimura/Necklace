#ifndef MBM_H
#define MBM_H

#include "ds.h"

int bfs(Hm *km, u64 *ka, const u64 *mu, const u64 *mv, u64 *dt, int k, u64 N);
int dfs(Hm *km, u64 *ka, u64 *mu, u64 *mv, u64 *dt, u64 u, int k);
u64 maximum_matching(Hm *km, u64 *ka, u64 *mu, u64 *mv, int k, u64 N);
void decompose(u64 *mu, u64 *mv, VV *cc, VV *pp, u64 N);

u64 baseline_a(const char* infile, int k, Hm **km, u64 **ka, VV *cc, VV *pp);
u64 bi_baseline_a(const char* infile, int k, Hm **km, u64 **ka, VV *cc, VV *pp, VVb *ccb, VVb *ppb);

void greedy_cover(Hm *km, u64 *ka, VV *cc, VV *pp, int k);
void bi_greedy_cover(Hm *km, u64 *ka, VV *cc, VV *pp, VVb *ccb, VVb *ppb, int k);

#endif // MBM_H