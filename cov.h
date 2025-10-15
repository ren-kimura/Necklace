#ifndef MBM_H
#define MBM_H

#include "ds.h"

int bfs(Hm *km, u64 *ka, const u64 *mu, const u64 *mv, u64 *dt, int k, u64 N);
int dfs(Hm *km, u64 *ka, u64 *mu, u64 *mv, u64 *dt, u64 u, int k);
u64 mbm(Hm *km, u64 *ka, u64 *mu, u64 *mv, int k, u64 N);
void decompose(u64 *mu, u64 *mv, VV *cc, VV *pp, u64 N);

u64 dextract(const char* infile, int k, Hm **km, u64 **ka, int di, VV *cc, VV *pp);

void disp_cp(u64 *ka, VV *cc, VV *pp, int k);

#endif // MBM_H