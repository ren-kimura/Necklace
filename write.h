#ifndef WRITE_H
#define WRITE_H

#include "ds.h"
#include "out.h"

u64 read_vle(FILE *fp);
char* rm_ext(const char* f);

void wrt(const char* f, const Rep* r, int k, int di, int cov, int out, u64 np);
void vread(const char* f);

#endif // WRITE_H