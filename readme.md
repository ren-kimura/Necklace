# About
File <max.cpp> exports edge-cetric/node-centric De Bruijn Graphs and one necklace cover with respect to FASTA-format input.

# Preliminary Setup
Run the following command in the same directory as <max.cpp>.

```make max.out```

# Usage
```./max.out [path to <input.fa>] [k] [option]```

- [path to <input.fa>] is a path to an input FASTA file.

- [k] is $k$ of $k$-mer.

- [option] is an option value, where:
    - 0 : outputs plain text representations without pointers
    - 1 : outputs representations with pointers not being sorted
    - 2 : outputs representations with sorted pointers (expected to be the smallest among these in size)
