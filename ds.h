#ifndef DS_H
#define DS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "uthash.h"

#define INF UINT64_MAX
typedef uint64_t u64;
extern char B[4];

//--- vector of bool ---
typedef struct {
    bool *data;
    size_t size;
    size_t cap;
} Vb;

void init_vb(Vb *v);
void push_backb(Vb *v, bool b);
void free_vb(Vb *v);

//--- vector of uint64_t ---
typedef struct {
    u64 *data;
    size_t size;
    size_t cap;
} V;

void init_v(V *v);
void push_back(V *v, u64 val);
u64 pop_back(V *v);
void free_v(V *v);

//--- vector of V ---
typedef struct {
    V *vs;
    size_t size;
    size_t cap;
} VV;

void init_vv(VV *vv);
void push_backv(VV *vv, V v);
void free_vv(VV *vv);
void print_vv(const char *title, VV *vv);

//---Stack and Queue of uint64_t---
typedef struct El {
    u64 data;
    struct El *next;
} El;

typedef struct {
    El *top;
} St;

void init_st(St *s);
int is_empty_st(St *s);
void push(St *s, u64 val);
u64 pop(St *s);

typedef struct {
    El *front;
    El *rear;
} Q;

void init_q(Q *q);
int is_empty_q(Q *q);
void enq(Q *q, u64 val);
u64 deq(Q *q);

//---Hash set of u64---
typedef struct {
    u64 key;
    UT_hash_handle hh;
} Hs;

int add_hs(Hs **s, u64 k);
int find_hs(Hs *s, u64 k);
void del_hs(Hs **s, u64 k);
void free_hs(Hs **s);

//---Hash map from u64 to u64---
typedef struct {
    u64 key;
    u64 val;
    UT_hash_handle hh;
} Hm;

int add_hm(Hm **m, u64 k, u64 v);
u64 find_hm(Hm *m, u64 k);
void del_hm(Hm **m, u64 k);
void free_hm(Hm **m);

#endif // DS_H