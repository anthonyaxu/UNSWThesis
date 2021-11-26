#!/usr/bin/python3

import subprocess
import os
import re

cell = [
    "HiCPrEC HiCLNCaP",
    "HiCPrEC HiCPC3",
    "HiCLNCaP HiCPC3"
]

resolution = [
    "5000",
    "10000",
    "25000",
    "40000",
    "500000"
]

for c in cell:
    for r in resolution:
        cs = c.split(" ")
        a = cs[0]
        b = cs[1]

        tad_path_a = "data/{}/{}/TADs/".format(a, r)
        tad_path_b = "data/{}/{}/TADs/".format(b, r)

        file_a = tad_path_a + "{}_{}.tad".format(a, r)
        file_b = tad_path_b + "{}_{}.tad".format(b, r)

        out_file = "output/overlaps/{}_{}_{}_overlap.bed".format(a, b, r)

        subprocess.run("bedtools intersect -a {} -b {} -bed > {}".format(file_a, file_b, out_file), shell = True)

cells = "HiCPrEC HiCLNCaP HiCPC3"
for r in resolution:
    cs = cells.split(" ")
    a = cs[0]
    b1 = cs[1]
    b2 = cs[2]

    tad_path_a = "data/{}/{}/TADs/".format(a, r)
    tad_path_b1 = "data/{}/{}/TADs/".format(b1, r)
    tad_path_b2 = "data/{}/{}/TADs/".format(b2, r)

    file_a = tad_path_a + "{}_{}.tad".format(a, r)
    file_b1 = tad_path_b1 + "{}_{}.tad".format(b1, r)
    file_b2 = tad_path_b2 + "{}_{}.tad".format(b2, r)

    out_file = "output/overlaps/{}_overlap.bed".format(r)

    subprocess.run("bedtools intersect -a {} -b {} {} -bed > {}".format(file_a, file_b1, file_b2, out_file), shell = True)

total_gene_length = 3000000000

overlap_dir = "output/overlaps/"
for r, d, f in os.walk(overlap_dir):
    for file in f:
        resolution = int(re.search(r"\d+\d+", file).group())

        read_file = open(overlap_dir + file, "r")
        contents = read_file.readlines()
        read_file.close()

        sum = 0
        for line in contents:
            l = line.split("\t")
            start = int(l[1])
            end = int(l[2])
            diff = (end - start)*resolution
            sum += diff
            
        percent = diff/total_gene_length*100
        print("{}: {}, {:.3f}%".format(file, diff, percent))