## Licenses

- This project uses `uthash`, which is licensed under the 1-clause BSD license. The license text is included in the `uthash.h` file.

## About
This project implements uni-directed/bi-directed node-centric de Bruijn Graphs and generates a cycle/path cover and a necklace cover from FASTA-format input.

## Building
Run the following command in the same directory as `nkl.c` to compile the code:

```bash
make
```

This will generate the executable(s) `nkl` (and `vfa`).

## Usage
To run the program, use the following command after adding a PATH to the directory (or simply specify the full or relative path to the execution file):
### generate a representation
```bash
nkl -i [in.fa] -k [k] -a [a] (-p) (-u)
```
- `-i`: set a path to an input FASTA file
- `-o`: specify output str filename (optional, default:`in_k.str`)
- `-k`: choose k-mer length $k$ s.t. $2\leq k\leq 31$
- `-a`: specify an algorithm to be run (eu:Eulertigs(default) fg:FullGreedy (needs -u -p) gb:GreedyBaseline (needs -p) ba:BaselineA gc:GreedyCover)
- `-p`: parenthesis representation (optional, computed with necklace_cover2)
- `-P`: parenthesis representation (optional, computed with necklace_cover(alg1))
- `-u`: distinguish a k-mer and its reverse complement (optional)

- `-h`: print help

### verify an output
```bash
nkl -i [in.fa] -k [k] (-u) -v [target.str]
```
- `-v`: verify the specified output

## Troubleshoot
Our implementation uses recursion in the source codes. It might cause stack overflow in some environment with some specific setting. In that case, the following is worth trying.

```bash
ulimit -s unlimited
```

## Requirements
- GCC supporting the C99 standard or later
- Standard C libraries
- A POSIX-compliant system (e.g., Linux, macOS) for functions like `getopt` (or WSL in Windows)