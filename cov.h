#ifndef MBM_H
#define MBM_H

#include "ds.h"

int bfs(Hm *km, u64 *ka, const u64 *mu, const u64 *mv, u64 *dt, u64 k, u64 N);
int dfs(Hm *km, u64 *ka, u64 *mu, u64 *mv, u64 *dt, u64 u, u64 k);
u64 mbm(Hm *km, u64 *ka, u64 *mu, u64 *mv, u64 k, u64 N);
void decompose(u64 *mu, u64 *mv, VV *cc, VV *pp, u64 N);

#endif // MBM_H