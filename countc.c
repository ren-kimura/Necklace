#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFFER_SIZE (1024 * 1024) // 1MB

// loopup table
static int base_table[256];
static int target_table[256];
static int delim_table[256];

void init_tables() {
    for (int i = 0; i < 256; i++) {
        base_table[i] = -1;
        target_table[i] = 0;
        delim_table[i] = 0;
    }
    base_table['A'] = base_table['a'] = 0;
    base_table['C'] = base_table['c'] = 1;
    base_table['G'] = base_table['g'] = 2;
    base_table['T'] = base_table['t'] = 3;
    
    target_table['A'] = target_table['a'] = 1;
    target_table['C'] = target_table['c'] = 1;
    target_table['G'] = target_table['g'] = 1;
    target_table['T'] = target_table['t'] = 1;
    target_table['('] = target_table[')'] = 1;

    delim_table['\n'] = delim_table['\r'] = delim_table[','] = delim_table['$'] = 1;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s <path> <opt=0|1> <k>\n", argv[0]);
        return 1;
    }

    init_tables();
    const char *filename = argv[1];
    int opt = atoi(argv[2]);
    int k = (opt == 1 && argc == 4) ? atoi(argv[3]) : 0;

    FILE *fp = fopen(filename, "rb"); // binary mode
    if (!fp) { perror("Error opening file"); return 1; }

    unsigned long long count = 0, implicit_cycle_count = 0;
    unsigned char *buf = malloc(BUFFER_SIZE);
    size_t bytes_read;
    
    int is_header_line = 0;
    int is_first_char_of_line = 1;
    
    // vars for opt=1
    unsigned long long pre = 0, suf = 0;
    unsigned long long mask = (opt == 1) ? (1ULL << (2 * (k - 1))) - 1 : 0;
    unsigned long long token_char_count = 0;

    while ((bytes_read = fread(buf, 1, BUFFER_SIZE, fp)) > 0) {
        for (size_t i = 0; i < bytes_read; i++) {
            unsigned char c = buf[i];

            if (is_first_char_of_line) {
                is_header_line = (c == '>');
                is_first_char_of_line = 0;
            }

            if (is_header_line) {
                if (c == '\n') is_first_char_of_line = 1;
                continue;
            }

            if (opt == 1) {
                int val = base_table[c];
                if (delim_table[c] || val == -1) {
                    if (token_char_count > 0) {
                        count += token_char_count;
                        if (token_char_count >= (size_t)k - 1 && pre == suf) {
                            count -= (k - 1);
                            implicit_cycle_count++;
                        }
                    }
                    token_char_count = 0;
                    pre = suf = 0;
                    if (c == '\n') is_first_char_of_line = 1;
                } else {
                    token_char_count++;
                    if (token_char_count < (unsigned long long)k) {
                        pre = (pre << 2) | val;
                        suf = pre;
                    } else {
                        suf = ((suf << 2) | val) & mask;
                    }
                }
            } else { // opt == 0
                if (c == '\n') {
                    is_first_char_of_line = 1;
                } else if (target_table[c]) {
                    count++;
                }
            }
        }
    }

    // opt=1
    if (opt == 1 && token_char_count > 0) {
        count += token_char_count;
        if (token_char_count >= (size_t)k - 1 && pre == suf) {
            count -= (k - 1);
            implicit_cycle_count++;
        }
    }

    free(buf);
    fclose(fp);
    printf("Total base count: %llu\n", count);
    if (opt) printf("Implicit cycle count: %llu\n", implicit_cycle_count);

    return 0;
}
