# About
This project implements edge-centric and node-centric De Bruijn Graphs and generates a cycle/path cover from FASTA-format input.

# Building
Run the following command in the same directory as `max.cpp` to compile the code:

```bash
make
```

This will generate the executable `max.out`.

# Usage
To run the program, use the following command:

```bash
./max.out [path to <input.fa>] [k] [option]
```

- `[path to <input.fa>]` is a path to an input FASTA file.

- `[k]` is the length of the k-mer $(2\leq k\leq 32)$.

- `[option]` is the output option:
    - 0 : Outputs plain text representations without pointers.
    - 1 : Outputs representations with unsorted pointers.
    - 2 : Outputs representations with sorted pointers (expected to be the smallest in size).

# Example
```bash
./max.out sample.fa 11 2
```
This command processes the input file `sample.fa` with k-mer length 11 and generates a representation with sorted pointers.

# Requirements
- C++17 or later
- Standard C++ libraries
- POSIX compliant system for memory management functions