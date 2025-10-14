#ifndef EUTILS_H
#define EUTILS_H

#include "ds.h"
#include "utils.h"
#include "uthash.h"
#include "stat.h"

typedef struct {
    u64 key; // k-1-mer
    uint8_t a; // bit mask representing adjacency (in/out)
    uint8_t o; // outdegree including dummy edges
    V dm; // destination of dummy edges from this k-1-mer
    UT_hash_handle hh;
} Node;

Node* build(Hm *km, int k, u64 N);
void balance(Node **g);
void etigs(Node **g, VV *tt, int k);
char** spell(VV *tt, int k, size_t *ns);

void free_ss(char **ss, size_t ns);
void free_g(Node **g);

#endif // EUTILS_H

