#!/usr/bin/env python3

import os
import sys

infile = sys.argv[1]

gene_names = []
cov = []
span = []

for line in open(infile, "r"):

    if line.startswith("chr"):

        line = line.strip().split()
       # print(line)
        gene_info = line[5]
        index1 = gene_info.find(";")
        gene = gene_info[8:index1]
       # print(gene)
        gene_names.append(gene)

        cov.append(int(line[-1]))

        span.append(int(line[2]) - int(line[1])+1)

# print(gene_names)

gene_names_unique = sorted(set(gene_names))

# print(gene_names_unique)
# print(cov)

for i in gene_names_unique:

    region_count = gene_names.count(i)

    index2 = gene_names.index(i)

    gene_cov_total = sum(cov[index2:index2+region_count])

    span_total = sum(span[index2:index2+region_count])

    if span_total != 0:
        gene_cov = round(gene_cov_total/span_total, 1)

        print(str(gene_cov))

    else:
        print("0")


# g = open("gene_symbols.txt", "w")

# g.write("\n".join(gene_names_unique))
