#define _POSIX_C_SOURCE 200809L
#ifndef UTILS_H
#define UTILS_H

#include "ds.h"

u64 next_pos(const char *read, u64 start, int k);
u64 enc(const char *s, int k);
void dec(u64 h, int k, char *s);
char* dec_base(u64 h);
void proc_sq(const char *sq, int k, Hm **km, u64 *id);
u64 extract(const char* infile, int k, Hm **km, u64 **ka);
u64 step(Hm *km, const u64 *ka, const int k, u64 id, int c, int is_forward);

#endif // UTILS_H