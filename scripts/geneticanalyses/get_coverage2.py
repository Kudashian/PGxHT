#!/usr/bin/env python3

import os
import sys

infile = sys.argv[1]

gene_names = []
cov = []
span = []

amp_cov = []

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

        amp_cov.append(str(round(int(line[-1])/(int(line[2])-int(line[1])),1)))
# print(gene_names)

print("\n".join(amp_cov))
