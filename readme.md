# About
This project implements edge-centric and node-centric De Bruijn Graphs and generates a cycle/path cover from FASTA-format input.

# Building
Run the following command in the same directory as `rep.cpp` to compile the code:

```bash
make
```

This will generate the executable `rep.out`.

# Usage
To run the program, use the following command:

```bash
./rep.out [flag] [path to <input.fa>] [k] [option]
```

- set `[flag]` with respect to the purpose of computation
    - 0 : Computes representations of node-centric path/cycle cover.
    - 1 : Computes the cumulative length of Eulertigs.

- `[path to <input.fa>]` is a path to an input FASTA file.

- `[k]` is the length of the k-mer $(2\leq k\leq 32)$.

- `[option]` is the output option (only needed when `[flag] = 0`):
    - 0 : Outputs plain text representation without pointers.
    - 1 : Outputs representation with unsorted pointers.
    - 2 : Outputs representation with sorted pointers (expected to be the smallest in size).
    - 3 : Outputs tree representation with BP.

# Example
```bash
./rep.out 0 input.fa 11 2
```
This command processes the input file `input.fa` with k-mer length 11 and generates a representation with sorted pointers.

```bash
./rep.out 1 input.fa 11
```
This command processes the input file `input.fa` with k-mer length 11 and computes the cumulative length of its Eulertigs (not bidirected).

# Requirements
- C++17 or later
- Standard C++ libraries
- POSIX compliant system for memory management functions