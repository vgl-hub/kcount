# kcount

A kmer counting tool

**kcount** is a single fast kmer counting tool. It generates the kmer histogram.

## Installation

Either download one of the releases or `git clone https://github.com/vgl-hub/kcount.git --recursive` and `make -j` in `kcount` folder.

## Usage

`kcount -f input.[fasta|fastq|gfa][.gz] -k 21

To check out all options and flags use `kcount -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
kcount -f testFiles/random1.fasta -k 3
```

## How to cite

This tools is part of the **gfastar** tool suite. If you use **kcount** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460

