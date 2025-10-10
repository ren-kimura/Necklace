#include <stdio.h>

void prog(size_t n, size_t t, const char* s) {
    if (t == 0) return; // avoid division by zero
    if (t > 100) {
        size_t vl = t / 100;
        if (n != 0 && n != t && (n % vl) != 0) {
            return;
        }
    }
    int pc = (int)((n * 100) / t);
    printf("\r%s: %d%%", s, pc);
    fflush(stdout);
}

void fin(const char* s) {
    printf("\r%s: 100%%\n", s);
    fflush(stdout);
}
