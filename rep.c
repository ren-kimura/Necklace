#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <ctype.h>
#include "uthash.h"

#define INF UINT64_MAX
char B[4] = {'A', 'C', 'G', 'T'};

typedef struct {
    uint64_t *data;
    size_t size;
    size_t cap;
} Vec;

typedef struct {
    Vec *vecs;
    size_t size;
    size_t cap;
} Vvec;

void init_vec(Vec *v) {
    v->data = NULL;
    v->size = 0;
    v->cap = 0;
}

void init_vvec(Vvec *vv) {
    vv->vecs = NULL;
    vv->size = 0;
    vv->cap = 0;
}

void free_vec(Vec *v) {
    free(v->data);
    init_vec(v);
}

void free_vvec(Vvec *vv) {
    for (size_t i = 0; i < vv->size; i++) {
        free_vec(&vv->vecs[i]);
    }
    free(vv->vecs);
    init_vvec(vv);
}

void print_vvec(const char *title, Vvec *vv) {
    printf("\n--- Debug Output: %s ---\n", title);
    if (vv->size == 0) {
        printf("No items found.\n");
        printf("---------------------------------\n");
        return;
    }
    for (size_t i = 0; i < vv->size; i++) {
        Vec *curr_vec = &vv->vecs[i];
        printf("Item %zu (%zu elements): ", i, curr_vec->size);
        for (size_t j = 0; j < curr_vec->size; j++) {
            printf("%lu ", (uint64_t)curr_vec->data[j]);
        }
        printf("\n");
    }
    printf("--------------------------------\n");
}

void add_vec(Vvec *vv) {
    if (vv->size >= vv->cap) {
        size_t new_cap = (vv->cap == 0) ? 8 : vv->cap * 2;
        Vec *new_vecs = realloc(vv->vecs, new_cap * sizeof(Vec));
        if (new_vecs == NULL) {
            perror("Failed to reallocate memory for vectors\n");
            exit(EXIT_FAILURE);
        }
        vv->vecs = new_vecs;
        vv->cap = new_cap;
    }
    Vec *new_vec = &vv->vecs[vv->size];
    new_vec->data = NULL;
    new_vec->size = 0;
    new_vec->cap = 0;
    vv->size++;
}

void push_back(Vec *v, uint64_t val) {
    if (v->size >= v->cap) {
        size_t new_cap = (v->cap == 0) ? 8 : v->cap * 2;
        uint64_t *new_data = realloc(v->data, new_cap * sizeof(uint64_t));
        if (new_data == NULL) {
            perror("Failed to reallocate memory for data\n");
            exit(EXIT_FAILURE);
        }
        v->data = new_data;
        v->cap = new_cap;
    }
    v->data[v->size++] = val;
}

uint64_t pop_back(Vec *v) {
    if (v->size == 0) {
        fprintf(stderr, "Error: pop_back called on an empty vector");
        return INF;
    }
    v->size--;
    return v->data[v->size];
}

// for stack and queue
typedef struct Elem {
    uint64_t data;
    struct Elem *next;
} Elem;

// stack
typedef struct {
    Elem *top;
} Stack;

void init_stack(Stack *s) {
    s->top = NULL;
}

int is_empty_stack(Stack *s) {
    return (s->top == NULL);
}

void push(Stack *s, uint64_t val) {
    Elem *new_elem = (Elem *)malloc(sizeof(Elem));
    if (new_elem == NULL) {
        printf("memory allocation failed\n");
        return;
    }
    new_elem->data = val;
    new_elem->next = s->top;
    s->top = new_elem;
    printf("pushed %ld\n", val);
}

uint64_t pop(Stack *s) {
    if (is_empty_stack(s)) {
        printf("Error: stack is empty\n");
        return -1;
    }
    Elem *tmp = s->top;
    uint64_t val = tmp->data;
    s->top = s->top->next;
    free(tmp);
    return val;
}

// queue
typedef struct {
    Elem *front;
    Elem *rear;
} Queue;

void init_queue(Queue *q) {
    q->front = NULL;
    q->rear = NULL;
}

int is_empty_queue(Queue *q) {
    return (q->front == NULL);
}

void enqueue(Queue *q, uint64_t val) {
    Elem *new_elem = (Elem *)malloc(sizeof(Elem));
    if (new_elem == NULL) {
        printf("memory allocation failed\n");
        return;
    }
    new_elem->data = val;
    new_elem->next = NULL;

    if (is_empty_queue(q)) {
        q->front = new_elem;
        q->rear = new_elem;
    } else {
        q->rear->next = new_elem;
        q->rear = new_elem;
    }
    // printf("enqueued %ld\n", val);
}

uint64_t dequeue(Queue *q) {
    if (is_empty_queue(q)) {
        printf("Error: queue is empty\n");
        return -1;
    }
    Elem *tmp = q->front;
    uint64_t val = tmp->data;
    q->front = q->front->next;
    if (q->front == NULL) {
        q->rear = NULL;
    }
    free(tmp);
    return val;
}

// hashset
typedef struct {
    uint64_t key;
    UT_hash_handle hh;
} set_t;

int set_add(set_t **set, uint64_t k) {
    set_t *e;
    HASH_FIND(hh, *set, &k, sizeof(k), e);
    if (!e) {
        e = malloc(sizeof *e);
        e->key = k;
        HASH_ADD(hh, *set, key, sizeof(k), e);
        return 1;
    }
    return 0;
}

int set_get(set_t *set, uint64_t k) {
    set_t *e;
    HASH_FIND(hh, set, &k, sizeof(k), e);
    return (e != NULL);
}

void free_set(set_t **set) {
    set_t *cur, *tmp;
    HASH_ITER(hh, *set, cur, tmp) {
        HASH_DEL(*set, cur);
        free(cur);
    }
}

// hashmap
typedef struct {
    uint64_t key;
    uint64_t val;
    UT_hash_handle hh;
} map_t;

int map_add(map_t **map, uint64_t k, uint64_t v) {
    map_t *e;
    HASH_FIND(hh, *map, &k, sizeof(k), e);
    if (!e) {
        e = malloc(sizeof *e);
	    e->key = k;
        e->val = v;
	    HASH_ADD(hh, *map, key, sizeof(k), e);
        return 1; // newly added
    }
    return 0; // already exists
}

uint64_t map_get(map_t *map, uint64_t k) {
    map_t *e;
    HASH_FIND(hh, map, &k, sizeof(k), e);
    if (!e) return INF;
    return e->val;
}

void map_del(map_t **map, uint64_t k) {
    map_t *e;
    HASH_FIND(hh, *map, &k, sizeof(k), e);
    if (e) { HASH_DEL(*map, e); free(e); }
}

void free_map(map_t **map) {
    map_t *cur, *tmp;
    HASH_ITER(hh, *map, cur, tmp) {
        HASH_DEL(*map, cur);
	    free(cur);
    }
}

typedef struct {
    char* str;
    uint64_t* arr;
} Rep;

// usage
static void usage(const char *prog) {
    fprintf(stderr,
	        "Usage: %s -i <file> -k <int> -o <0|1|2|10>\n"
	        "\t-i FILE\t\tinput FASTA file\n"
	        "\t-k INT\t\tk-mer length (>=1 && <=31)\n"
	        "\t-o OPTION\t0:naive 1:pointer 2:pseudoforest\n"
	        "\t\t\t10:eulertigs(uni-directed)\n",
	        prog);
    exit(EXIT_FAILURE);
}

static int parse_int(const char *s) {
    char *end;
    errno = 0;
    int v = strtol(s, &end, 10);
    if (errno || *end) usage("./rep");
    return v;
}

uint64_t enc(const char *s, int k) {
    uint64_t h = 0;
    for (int i = 0; i < k; ++i) {
        char c = toupper(s[i]);
        if      (c == 'A')  h = (h << 2) | 0; // A -> 00
        else if (c == 'C')  h = (h << 2) | 1; // C -> 01
        else if (c == 'G')  h = (h << 2) | 2; // G -> 02
        else if (c == 'T')  h = (h << 2) | 3; // T -> 03
        else return UINT64_MAX;
    }
    return h;
}

void dec(uint64_t h, int k, char *s) {
    for (int i = 0; i < k; i++) {
        s[k - i - 1] = B[h & 3];
        h >>= 2;
        s[k] = '\0';
    }
}

char* dec_base(uint64_t h) {
    if (h == 0)         return "A";
    else if (h == 1)    return "C";
    else if (h == 2)    return "G";
    else                return "T";
}

long next_valid(const char *read, long start, long k) {
    long len = strlen(read);
    long count = 0;
    for (long j = start; j < len; j++) {
        char c = toupper(read[j]); // to uppercase
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            ++count;
            if (count == k) return j - k + 1;
        } else {
            count = 0; // reset count if non-ACGT char found
        }
    }
    return -1; // no valid k-mer found after `start`
}

uint64_t extract_kmers(const char* infile, int k, map_t **kmap, uint64_t **karr) {
    FILE *fp = fopen(infile, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open file %s\n", infile);
        exit(EXIT_FAILURE);
    }
    printf("file %s opened\n", infile);

        uint64_t id = 0;
    char *line = NULL;
    size_t len = 0;
    ssize_t read_len;

    while ((read_len = getline(&line, &len, fp)) != -1) {
        if (read_len > 0 && line[read_len - 1] == '\n') {
            line[read_len - 1] = '\0';
        }
        if (read_len > 1 && line[read_len - 2] == '\r') {
            line[read_len - 2] = '\0';
        }
        if (strlen(line) == 0 || line[0] == '>') {
            continue;
        }
        long line_len = strlen(line);
        if (line_len < k) {
            continue; // too short to contain a k-mer
        }
        uint64_t h;
        char s[k + 1]; // buffer for k-mer string
        long j = next_valid(line, 0, k);
        if (j != -1) {
            strncpy(s, line + j, k);
            s[k] = '\0';
            h = enc(s, k);
            if (map_add(&kmap, h, id)) id++;
            uint64_t m = (1ULL << (2 * (k - 1))) - 1; // clear two MSBs
            for (++j; j <= line_len - k; ++j) {
                char c = toupper(line[j + k - 1]);
                if (c == 'A')       h = ((h & m) << 2) | 0;
                else if (c == 'C')  h = ((h & m) << 2) | 1;
                else if (c == 'G')  h = ((h & m) << 2) | 2;
                else if (c == 'T')  h = ((h & m) << 2) | 3;
                else {
                    j = next_valid(line, j + 1, k);
                    if (j == -1) break; // no more k-mer on this line
                    strncpy(s, line + j, k);
                    s[k] = '\0';
                    h = enc(s, k);
                }
                if (map_add(&kmap, h, id)) id++;
            }
        }
    }
    const uint64_t N = (uint64_t)HASH_COUNT(*kmap); // number of k-mers

    printf("total unique k-mers = %ld\n", N);
    /* display kmap */
    // /*
    printf("key\t\tdec(key)\tvalue\n");
    printf("-------------------------------------------\n");
    for (map_t *s = kmap; s != NULL; s = (map_t*)(s->hh.next)) {
        char t[k + 1];
        dec(s->key, k, t);
        t[k] = '\0';
        printf("%lu\t\t%s\t\t%lu\n", s->key, t, s->val);
    }
    printf("-------------------------------------------\n");
    // */

    *karr = malloc(N * sizeof(uint64_t));
    if (karr == NULL) {
        fprintf(stderr, "Memory allocation failed for karr\n");
        exit(EXIT_FAILURE);
    }
    for (map_t *s = *kmap; s != NULL; s = (map_t*)(s->hh.next)) {
        (*karr)[s->val] = s->key;
    }

    free(line);
    fclose(fp);
    printf("file %s closed\n", infile);

    return N;
}

uint64_t step(map_t *kmap, const uint64_t *karr, const int k,
              uint64_t id, int c, int is_forward) {
    uint64_t h = karr[id];
    if (h == INF) return INF;
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return INF;
    if (is_forward) {
        uint64_t m = (1ULL << (2 * (k - 1))) - 1;
        h = (h & m) << 2;
        if (c == 'C')       h |= 1;
        else if (c == 'G')  h |= 2;
        else if (c == 'T')  h |= 3;
    } else {
        h >>= 2;
        if (c == 'C')       h |= (1ULL << (2 * (k - 1)));
        else if (c == 'G')  h |= (2ULL << (2 * (k - 1)));
        else if (c == 'T')  h |= (3ULL << (2 * (k - 1)));
    }

    uint64_t h_new = map_get(kmap, h);
    if (h_new == INF || h_new == h) return INF; // no branch to c
    
    return h_new;
}

int bfs(map_t *kmap, uint64_t *karr, const uint64_t *mu,
        const uint64_t *mv, uint64_t *dt, uint64_t k, uint64_t N) {
    Queue q;
    init_queue(&q);
    int found = 0;
    for (uint64_t u = 0; u < N; u++) {
        if (mu[u] == INF) {
            dt[u] = 0;
            enqueue(&q, u);
        } else dt[u] = INF;
    }
    while (!is_empty_queue(&q)) {
        uint64_t u = dequeue(&q);
        for (int i = 0; i < 4; i++) {
            char c = B[i];
            uint64_t v = step(kmap, karr, k, u, c, 1);
            if (v == INF) continue;
            if (mv[v] == INF) {
                found = 1;
            } else if (dt[mv[v]] == INF) {
                dt[mv[v]] = dt[u] + 1;
                enqueue(&q, mv[v]);
            }
        }
    }
    return found;
}

int dfs(map_t *kmap, uint64_t *karr, uint64_t *mu, uint64_t *mv,
        uint64_t *dt, uint64_t u, uint64_t k) {
    for (int i = 0; i < 4; i++){
        char c = B[i];
        uint64_t v = step(kmap, karr, k, u, c, 1);
        if (v == INF) continue;
        if (mv[v] == INF || (dt[mv[v]] == dt[u] + 1 && 
            dfs(kmap, karr, mu, mv, dt, mv[v], k))) {
            mu[u] = v;
            mv[v] = u;
            return 1;
        }
    }
    dt[u] = INF; // found no augpath
    return 0;
}

uint64_t hopcroft_karp(map_t *kmap, uint64_t *karr, uint64_t *mu,
                       uint64_t *mv, uint64_t k, uint64_t N) {
    uint64_t *dt = malloc(N * sizeof(uint64_t));
    if (dt == NULL) {
        printf("Memory allocation failed for array dt\n");
        free(dt);
        return INF;
    }

    uint64_t M = 0;
    for (uint64_t i = 0; i < N; i++) {
        mu[i] = INF;
        mv[i] = INF;
    }
    while (bfs(kmap, karr, mu, mv, dt, k, N)) {
        for (uint64_t u = 0; u < N; u++) {
            if (mu[u] == INF && dfs(kmap, karr, mu, mv, dt, u, k)) {
                M++; // found an augpath
            }
        }
    }
    printf("matching size %ld\n", M);
    free(dt);
    return M;
}

void decompose(uint64_t *mu, uint64_t *mv, Vvec *cc, Vvec *pp, uint64_t N){
    uint64_t seen_u[N]; uint64_t seen_v[N];
    for (uint64_t i = 0; i < N; i++) seen_u[i] = seen_v[i] = 0;
    for (uint64_t u = 0; u < N; u++) {
        if (seen_u[u]) continue;
        Vec fv; init_vec(&fv);
        uint64_t crnt = u;
        while (1) {
            if (seen_u[crnt]) break;
            seen_u[crnt] = 1;
            push_back(&fv, crnt);
            uint64_t next = mu[crnt];
            if (next == INF) break;
            seen_v[next] = 1;
            crnt = next;
        }
        if (fv.size > 1 && mu[fv.data[fv.size - 1]] == fv.data[0]) {
            add_vec(cc);
            Vec *c = &cc->vecs[cc->size - 1];
            for (uint64_t i = 0; i < fv.size; i++) {
                push_back(c, fv.data[i]);
            }
        } else {
            Vec bv; init_vec(&bv);
            crnt = u;
            while (1) {
                if (seen_v[crnt]) break;
                seen_v[crnt] = 1;
                push_back(&bv, crnt);
                uint64_t prev = mv[crnt];
                if (prev == INF) break;
                seen_u[prev] = 1;
                crnt = prev;
            }
            add_vec(pp);
            Vec *p = &pp->vecs[pp->size - 1];
            if (bv.size > 0) {
                for (uint64_t i = bv.size - 1; i > 0; i--) {
                    push_back(p, bv.data[i]);
                }
            }
            for (uint64_t i = 0; i < fv.size; i++) {
                push_back(p, fv.data[i]);
            }
            free(bv.data);
        }
        free(fv.data);
    }
    return;
}

Rep plain(uint64_t *karr, Vvec *cc, Vvec *pp, uint64_t k, uint64_t N) {
    Rep r;
    r.str = (char*)malloc(cc->size + pp->size + N + pp->size * (k - 1));
    if (r.str == NULL) {
        perror("Failed to allocate memory in plain()\n");
        exit(1);
    }
    r.str[0] = '\0';
    r.arr = NULL;
    printf("Generating plaintext representation...\n");
    for (size_t i = 0; i < cc->size; i++) {
       Vec *curr_vec = &cc->vecs[i];
       for (size_t j = 0; j < curr_vec->size; j++) {
           strcat(r.str, dec_base(karr[curr_vec->data[j]] % 4));
       }
       strcat(r.str, ",");
    }
    strcat(r.str, ",");
    for (size_t i = 0; i < pp->size; i++) {
        Vec *curr_vec = &pp->vecs[i];
        char s[k + 1]; s[k] = '\0';
        dec(karr[curr_vec->data[0]], k, s);
        strcat(r.str, s);
        for (size_t j = 1; j < curr_vec->size; j++) {
            strcat(r.str, dec_base(karr[curr_vec->data[j]] % 4));
        }
        strcat(r.str, ",");
    }
    if (strlen(r.str)) r.str[strlen(r.str) - 1] = '\0';
    if (pp->size == 0) r.str[strlen(r.str) - 1] = '\0';
    return r;
}

Rep unsorted(map_t *kmap, uint64_t *karr, Vvec *cc, Vvec *pp, uint64_t k, uint64_t N) {
    Rep r;
    r.str[0] = '\0';
    uint64_t pred[pp->size];
    for (uint64_t i = 0; i < pp->size; i++) { pred[i] = INF; }
    for (uint64_t i = 0; i < pp->size; i++) {
        for (int c_idx = 0; c_idx < 4; c_idx++) {
            if (step(kmap, karr, k, pp->vecs[i].data[0], B[c_idx], 0) != INF) {
                break;
            }
        }
    }
    map_t *hd = NULL;
    set_t *rp = NULL;
    for (uint64_t i = 0; i < pp->size; i++) {
        int is_rp = 1;
        uint64_t hid = pp->vecs[i].data[0];
        for (int c_idx = 0; c_idx < 4; c_idx++) {
            if (step(kmap, karr, k, hid, B[c_idx], 0) != INF) {
                is_rp = 0;
                break;
            }
        }
        if (is_rp) {
            set_add(&rp, i);
        } else {
            map_add(&hd, hid, i);
        }
    }
    // pointer vector
    r.arr = (uint64_t*)malloc(sizeof(uint64_t) * pp->size);
    for (uint64_t i = 0; i < pp->size; i++) { r.arr[i] = -1; }

    uint64_t pos = 0;
    for (uint64_t i = 0; i < cc->size; i++) {
        for (uint64_t j = 0; j < cc->vecs[i].size; j++) {
            for (int c_idx = 0; c_idx < 4; c_idx++) {
                uint64_t next = step(kmap, karr, k, cc->vecs[i].data[j], B[c_idx], 1);
                if (next == INF) continue;
                uint64_t npid = map_get(hd, next);
                if (npid == INF) continue;
                r.arr[npid] = pos;
                map_del(&hd, next);
                if(hd == NULL) goto end;
            } pos++;
        } pos++;
    }
    for (uint64_t i = 0; i < pp->size; i++) {
        for (uint64_t j = 0; j < pp->vecs[i].size; j++) {
            for (int c_idx = 0; c_idx < 4; c_idx++) {
                uint64_t next = step(kmap, karr, k, pp->vecs[i].data[j], B[c_idx], 1);
                if (next == INF) continue;
                uint64_t npid = map_get(hd, next);
                if (npid == INF) continue;
                r.arr[npid] = pos;
                map_del(&hd, next);
                if (hd == NULL) goto end;
            } pos++;
        } pos++;
    }
    end:
    r.str[0] = '\0';
    for (size_t i = 0; i < cc->size; i++) {
       Vec *curr_vec = &cc->vecs[i];
       for (size_t j = 0; j < curr_vec->size; j++) {
           strcat(r.str, dec_base(karr[curr_vec->data[j]] % 4));
       }
       strcat(r.str, ",");
    }
    strcat(r.str, ",");
    for (size_t i = 0; i < pp->size; i++) {
        Vec *curr_vec = &pp->vecs[i];
        char s[k + 1]; s[k] = '\0';
        dec(karr[curr_vec->data[0]], k, s);
        strcat(r.str, s);
        for (size_t j = 1; j < curr_vec->size; j++) {
            strcat(r.str, dec_base(karr[curr_vec->data[j]] % 4));
        }
        strcat(r.str, ",");
    }
    if (strlen(r.str)) r.str[strlen(r.str) - 1] = '\0';
    if (pp->size == 0) r.str[strlen(r.str) - 1] = '\0';
    return r; 
}

typedef struct Frame {
    uint64_t current;
    uint64_t index;
    int b_index;
    struct Frame *next;
} Frame;

typedef struct {
    Frame *top;
} Fstack;

void init_fstack(Fstack *s) {
    s->top = NULL;
}

int is_empty_fstack(Fstack *s) {
    return (s->top == NULL);
}

void push_frame(Fstack *s, uint64_t current, uint64_t index, int b_index) {
    Frame *new_frame = (Frame*)malloc(sizeof(Frame));
    if (new_frame == NULL) {
        fprintf(stderr, "Memory allocation failed for new_frame\n");
        exit(EXIT_FAILURE);
    }
    new_frame->current = current;
    new_frame->index = index;
    new_frame->b_index = b_index;
    new_frame->next = s->top;
    s->top = new_frame;
    printf("pushed a frame {%ld, %ld, %d}\n", current, index, b_index);
}

void pop_frame(Fstack *s) {
    if (is_empty_fstack(s)) {
        printf("Error: stack is empty\n");
        return;
    }
    Frame *tmp = s->top;
    s->top = s->top->next;
    free(tmp);
}

Vec* find_new_cycle(map_t *kmap, uint64_t *karr, map_t *hd, int k, Vvec *pp, uint64_t *in_pord, uint64_t *is_leaf, uint64_t *visited, uint64_t start) {
    Vec* cycle = malloc(sizeof(Vec));
    if (cycle == NULL) {
        return NULL;
    }
    cycle->data = NULL;
    cycle->size = cycle->cap = 0;

    /* implement the process here */
    Fstack s;
    init_fstack(&s);

    visited[start] = 1;
    push_frame(&s, start, 0, 0);

    while (!is_empty_fstack(&s)) {
        // rm from stack if explored all the nodes in the current path
        if (s.top->index >= (uint64_t)pp->vecs[s.top->current].size) {
            is_leaf[s.top->current] = 1; // memo as a leaf
            pop_frame(&s);
            if (cycle->size != 0) pop_back(&cycle);
            continue;
        }
        // when entering a new node (on start of this frame)
        if (s.top->index == 0 && s.top->b_index == 0) push_back(&cycle, s.top->current);
        // try next base
        if (s.top->b_index < 4) {
            char c = B[s.top->b_index];
            s.top->b_index++;
            uint64_t current_node = pp->vecs[s.top->current].data[s.top->index];
            uint64_t next_node = step(kmap, karr, k, current_node, c, 1);
            if (next_node == INF) continue;
            uint64_t next_pid = map_get(hd, next_node);
            if (next_pid == INF) continue;
            if (in_pord[next_pid]) continue;
            if (next_pid == start) return cycle;
            if (visited[next_pid]) continue; // skip visited paths
            if (is_leaf[next_pid]) continue; // skip leaves  
            visited[next_pid] = 1;
            push_frame(&s, next_pid, 0, 0);          
        } else {
            s.top->index++;
            s.top->b_index = 0;
        }
    }
    init_vec(cycle);
    return cycle;
}

int main(int argc, char *argv[]) {
    const char *infile = NULL;
    int k, di, cover, mode, opt;
    k = di = cover = mode = -1;

	while ((opt = getopt(argc, argv, "i:k:d:c:o:")) != -1) {
	    switch (opt) {
		    case 'i':
			    if (infile) usage(argv[0]);
                infile = optarg;
                break;
            case 'k':
                if (k != -1) usage(argv[0]);
                k = parse_int(optarg);
                if (k < 1 || k > 31) usage(argv[0]);
                break;
            case 'd':
                if (di != -1) usage(argv[0]);
                di = parse_int(optarg);
                if (di != 0 && di != 1) usage(argv[0]);
                break;
            case 'c':
                if (cover != -1) usage(argv[0]);
                cover = parse_int(optarg);
                if (cover != 0 && cover != 1) usage(argv[0]);
                break;
            case 'o':
                if (mode != -1) usage(argv[0]);
                mode = parse_int(optarg);
                if (!(mode == 0 || mode == 1 || mode == 2 || mode == 3 || mode == 10))
                    usage(argv[0]);
                break;
            default:
                usage(argv[0]);
		}
	}

    if (!infile || k == -1 || di == -1 || cover == -1 || mode == -1) {
        usage(argv[0]);
    }

    printf("input = %s\n", infile);
    printf("k = %d\n", k);
    printf("di = %d\n", di);
    printf("cover = %d\n", cover);
    printf("mode = %d\n", mode);

    map_t *kmap = NULL;
    map_t *karr = NULL;
    uint64_t N = extract_kmers(infile, k, &kmap, &karr);
    
    uint64_t *mu = malloc(N * sizeof(uint64_t));
    uint64_t *mv = malloc(N * sizeof(uint64_t));
    if (mu == NULL || mv == NULL) {
        printf("Memory allocation failed for array mu or mv\n");
        free(karr); free(mu); free(mv);
        free_map(&kmap);
        return -1;
    }
    uint64_t M = hopcroft_karp(kmap, karr, mu, mv, k, N);
    if (M == INF) {
        printf("Exiting due to bad memory allocation\n");
        free(karr); free(mu); free(mv);
        free_map(&kmap);
        return -1;
    }
    Vvec cc, pp;
    init_vvec(&cc); init_vvec(&pp);
    decompose(mu, mv, &cc, &pp, N);
    print_vvec("Cycles (cc)", &cc);
    print_vvec("Paths (pp)", &pp);
    
    if (mode == 0) {
        Rep r = plain(karr, &cc, &pp, k, N);
        printf("Generated plaintext:\n%s\n", r.str);
    } else if (mode == 1) {
        Rep r = unsorted(kmap, karr, &cc, &pp, k, N);
        printf("Generated text:\n%s\n", r.str);
    } else if (mode == 2) {
        printf("under construction\n");
    } else if (mode == 3) {
        printf("under construction\n");
    } else if (mode == 10) {
        printf("under construction\n");
    } else {
        printf("Error: invalid mode for output format\n");
        exit(1);
    }

    free(karr); free(mu); free(mv);
    free_map(&kmap);
    free_vvec(&cc); free_vvec(&pp);
    
	return 0;
}
