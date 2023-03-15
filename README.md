# kcount

A kmer counting tool

**kcount** is a single fast kmer counting tool. It generates a kmer histogram and a database. It can merge multiple kmer databases. It is based on parallel hashing.

## Installation

Either download one of the releases or `git clone https://github.com/vgl-hub/kcount.git --recursive` and `make -j` in the `kcount` folder.

## Usage

```
kcount count -f input.[fasta|fastq][.gz] -k 21
```

It accepts multiple files as input, separated by space. To check out all options and flags use `kcount -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
kcount union -d testFiles/random2.kc testFiles/random2.kc -o random2.hist
```

## How to cite

This tool is part of the **gfastar** tool suite. If you use **kcount** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460

