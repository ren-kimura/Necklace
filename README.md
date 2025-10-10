## Licenses

- This project uses `uthash`, which is licensed under the 1-clause BSD license. The license text is included in the `uthash.h` file.

## About
This project implements uni-directed/bi-directed node-centric de Bruijn Graphs and generates a cycle/path cover from FASTA-format input.

## Building
Run the following command in the same directory as `nkl.c` to compile the code:

```bash
make
```

This will generate the executable `nkl`.

## Usage
To run the program, use the following command after adding a PATH to the directory (or simply specify the full or relative path to the execution file):

```bash
nkl -i [infile] -k [2-31] -d [0|1] -c [0|1|2] -o [0|1|2]
```
- `-i`: set a path to an input FASTA file after this
- `-k`: choose $k$ s.t. $2\leq k\leq 31$
- `-d`: choose either 0: uni-directed or 1: bi-directed
- `-c`: choose how to cover the de Bruijn graph, 0: maximum matching (ONLY WHEN d == 0) 1: greedy linear scan of `infile` 2: greedy dfs from unvisited vertices
- `-o`: choose representation 0: flat 1: pointer 2: balanced parentheses 

## Example
```bash
nkl -i input.fa -k 11 -d 0 -c 0 -o 2
```
This command processes the input file `input.fa` with k-mer length 11, finds a path cover of uni-directed de Bruijn graph with maximum matching and generates a balanced parentheses representation.

```bash
nkl -i input.fa -k 31 -d 0 -c 1 -o 0
```
This command processes the input file `input.fa` with k-mer length 31, finds a path cover of uni-directed de Bruijn graph by greedy linear scan of `input.fa` and generates a flat representation.

```bash
nkl -i input.fa -k 11 -d 1 -c 2 -o 1
```
This command processes the input file `input.fa` with k-mer length 11, finds a path cover of bi-directed de Bruijn graph by greedy dfs and generates a pointer representation.

## Requirements
- A C compiler supporting the C99 standard or later (e.g., GCC, Clang)
- Standard C libraries
- A POSIX-compliant system (e.g., Linux, macOS) for functions like `getopt` (or WSL in Windows)