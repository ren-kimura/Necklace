#include "ds.h"
#include <string.h>

char B[4] = {'A', 'C', 'G', 'T'};

//---Vb---
void init_vb(Vb *v) {
    v->data = NULL;
    v->size = 0;
    v->cap = 0;
}

void push_backb(Vb *v, bool b) {
    if (v->size >= v->cap) {
        size_t ncap = (v->cap == 0) ? 8 : v->cap * 2;
        bool *ndata = realloc(v->data, ncap * sizeof(bool));
        if (ndata == NULL) {
            fprintf(stderr, "Error: malloc failed for new bool\n");
            exit(EXIT_FAILURE);
        }
        v->data = ndata;
        v->cap = ncap;
    }
    v->data[v->size++] = b;
}

bool pop_backb(Vb *v) {
    if (v->size == 0) {
        fprintf(stderr, "Error: pop_backb called on an empty vector\n");
        return false;
    }
    v->size--;
    return v->data[v->size];
}

void free_vb(Vb *v) {
    free(v->data);
    init_vb(v);
}

//---V---
void init_v(V *v) {
    v->data = NULL;
    v->size = 0;
    v->cap = 0;
}

void push_back(V *v, u64 val) {
    if (v->size >= v->cap) {
        size_t ncap = (v->cap == 0) ? 8 : v->cap * 2;
        u64 *ndata = realloc(v->data, ncap * sizeof(u64));
        if (ndata == NULL) {
            fprintf(stderr, "Error: malloc failed for new u64\n");
            exit(EXIT_FAILURE);
        }
        v->data = ndata;
        v->cap = ncap;
    }
    v->data[v->size++] = val;
}

u64 pop_back(V *v) {
    if (v->size == 0) {
        fprintf(stderr, "Error: pop_back called on an empty vector");
        return INF;
    }
    v->size--;
    return v->data[v->size];
}

void free_v(V *v) {
    free(v->data);
    init_v(v);
}

//---VV---
void init_vv(VV *vv) {
    vv->vs = NULL;
    vv->size = 0;
    vv->cap = 0;
}

void push_backv(VV *vv, V v) {
    if (vv->size >= vv->cap) {
        size_t ncap = (vv->cap == 0) ? 8 : vv->cap * 2;
        V *nvs = realloc(vv->vs, ncap * sizeof(V));
        if (nvs == NULL) {
            fprintf(stderr, "Error: malloc failed for new vector\n");
            exit(EXIT_FAILURE);
        }
        vv->vs = nvs;
        vv->cap = ncap;
    }
    V *nv = &vv->vs[vv->size];
    nv->size = v.size;
    nv->cap = v.size;
    if (v.size > 0) {
        nv->data = malloc(nv->cap * sizeof(u64));
        memcpy(nv->data, v.data, nv->size * sizeof(u64));
    } else {
        nv->data = NULL;
    }
    vv->size++;
}

void free_vv(VV *vv) {
    for (size_t i = 0; i < vv->size; i++) {
        free_v(&vv->vs[i]);
    }
    free(vv->vs);
    init_vv(vv);
}

//---St---
void init_st(St *s) {
    s->top = NULL;
}

int is_empty_st(St *s) {
    return (s->top == NULL);
}

void push(St *s, u64 val) {
    El *nel = (El*)malloc(sizeof(El));
    if (nel == NULL) {
        fprintf(stderr, "Error: malloc failed for new el\n");
        exit(EXIT_FAILURE);
    }
    nel->data = val;
    nel->next = s->top;
    s->top = nel;
}

u64 pop(St *s) {
    if (is_empty_st(s)) {
        fprintf(stderr, "Error: pop called on an empty stack\n");
        return INF;
    }
    El *t = s->top;
    u64 val = t->data;
    s->top = s->top->next;
    free(t);
    return val;
}

//---Q---
void init_q(Q *q) {
    q->front = NULL;
    q->rear = NULL;
}

int is_empty_q(Q *q) {
    return (q->front == NULL);
}

void enq(Q *q, u64 val) {
    El *nel = (El*)malloc(sizeof(El));
    if (nel == NULL) {
        fprintf(stderr, "Error: malloc failed for new el\n");
        return;
    }
    nel->data = val;
    nel->next = NULL;

    if (is_empty_q(q)) {
        q->front = nel;
        q->rear = nel;
    } else {
        q->rear->next = nel;
        q->rear = nel;
    }
}

u64 deq(Q *q) {
    if (is_empty_q(q)) {
        fprintf(stderr, "Error: deq called on an empty queue\n");
        return INF;
    }
    El *t = q->front;
    u64 val = t->data;
    q->front = q->front->next;
    if (q->front == NULL) {
        q->rear = NULL;
    }
    free(t);
    return val;
}

//---Hs---
int add_hs(Hs **hs, u64 k) {
    Hs *e;
    HASH_FIND(hh, *hs, &k, sizeof(k), e);
    if (!e) {
        e = malloc(sizeof(Hs));
        e->key = k;
        HASH_ADD(hh, *hs, key, sizeof(k), e);
        return 1;
    }
    return 0;
}

int find_hs(Hs *hs, u64 k) {
    Hs *e;
    HASH_FIND(hh, hs, &k, sizeof(k), e);
    return (e != NULL);
}

void del_hs(Hs **hs, u64 k) {
    Hs *e;
    HASH_FIND(hh, *hs, &k, sizeof(k), e);
    if (e) { HASH_DEL(*hs, e); free(e); }
}

void free_hs(Hs **hs) {
    Hs *cur, *tmp;
    HASH_ITER(hh, *hs, cur, tmp) {
        HASH_DEL(*hs, cur);
        free(cur);
    }
}

//---Hm---
int add_hm(Hm **hm, u64 k, u64 v) {
    Hm *e;
    HASH_FIND(hh, *hm, &k, sizeof(k), e);
    if (!e) {
        e = malloc(sizeof *e);
        e->key = k;
        e->val = v;
        HASH_ADD(hh, *hm, key, sizeof(k), e);
        return 1; // newly added
    }
    return 0; // already exists
}

u64 find_hm(Hm *hm, u64 k) {
    Hm *e;
    HASH_FIND(hh, hm, &k, sizeof(k), e);
    if (!e) return INF;
    return e->val; 
}

void del_hm(Hm **hm, u64 k) {
    Hm *e;
    HASH_FIND(hh, *hm, &k, sizeof(k), e);
    if (e) { HASH_DEL(*hm, e); free(e); }
}

void free_hm(Hm **hm) {
    Hm *cur, *tmp;
    HASH_ITER(hh, *hm, cur, tmp) {
        HASH_DEL(*hm, cur);
        free(cur);
    }
}