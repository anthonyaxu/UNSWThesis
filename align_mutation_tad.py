#!/usr/bin/python3

import sys, csv, re, glob

csv.field_size_limit(sys.maxsize)

# tad_file = "../data/ontad/output/ontad_healthy_40000.tad"
tad_file = "../data/ontad/output/ontad_cancer_40000.tad"

mutation_file = "../data/mutations/top_mutations_5.csv"

dict = {}
with open(tad_file) as file:
    f = csv.reader(file, delimiter = "\t")
    next(f)

    for row in f:
        chr = re.sub("chr", "", row[0])
        seq = [int(row[1]), int(row[2])]

        if not chr in dict:
            dict[chr] = [seq]
        else:
            dict[chr].append(seq)

f = open(mutation_file, "r")
rows = f.readlines()
content = rows[1:]
f.close()

# write_file = "../data/mutations/healthy_tad_mutations.csv"
write_file = "../data/mutations/cancer_tad_mutations.csv"
wfile = open(write_file, "w+")
wfile.writelines(rows[0])

for r in content:
    seq = r.split(",")[1]
    seq = seq.split(":")
    chr = seq[0]
    chr = re.sub("chr", "", chr)
    start = int(seq[1])
    end = int(seq[2])

    if chr in dict:
        seq_list = dict[chr]
        for s in seq_list:
            tad_start = s[0]
            tad_end = s[1]

            if start >= tad_start and end <= tad_end:
                wfile.writelines(r)

wfile.close()