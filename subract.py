#!/usr/bin/python3

import subprocess
import csv
import sys
import re

csv.field_size_limit(sys.maxsize)

# For 25000 resolution
RESOLUTION = 25000

def subtract_overlaps(a, b, outfile):
    subprocess.run("bedtools subtract -a {} -b {} > {}".format(a, b, outfile), shell = True)


def get_tads(file_path):
    t = {}
    with open(file_path) as file:
        f = csv.reader(file, delimiter = "\t")

        for row in f:
            chr = re.sub("chr", "", row[0])
            seq = [int(row[1])*RESOLUTION, int(row[2])*RESOLUTION]

            if not chr in t:
                t[chr] = [seq]
            else:
                t[chr].append(seq)
    return t


def get_mutations():
    m_dict = {}
    count_dict = {}
    with open("mutation_samples_count.csv") as csv_file:
            f_m = csv.reader(csv_file, delimiter = ",")
            # Skip csv header
            next(f_m)

            for row in f_m:
                id = row[0]
                seq = row[1].split(":")
                chr = re.sub("chr", "", seq[0])
                start = int(seq[1])
                end = int(seq[2])
                m_dict[id] = [chr, start, end]
                count = int(row[3])
                count_dict[id] = count

    count_dict = dict(sorted(count_dict.items(), key = lambda x:x[1], reverse = True))
    return m_dict, count_dict


def cancer_mutations(tads, mutations):
    c_mutations = []
    for m in mutations:
        seq = mutations[m]
        chr = seq[0]
        start = int(seq[1])
        end = int(seq[2])

        try:
            tds = tads[chr]

            for t in tds:
                tad_start = t[0]
                tad_end = t[1]
                if start >= tad_start and end <= tad_end:
                    c_mutations.append(m)
        except:
            continue

    return c_mutations


def top_cancer_mutations(count, mutations):
    cancer_muts = {}
    for c in count:
        if len(cancer_muts) > 50:
            break
        if c in mutations:
            cancer_muts[c] = count[c]
    return cancer_muts


def get_cnv_samples():
    dict = {}
    with open("copy_number_count.csv") as csv_file:
            f_cnvc = csv.reader(csv_file, delimiter = ",")
            # Skip csv header
            next(f_cnvc)

            for row in f_cnvc:
                sample_id = row[0]
                cnv = row[1].split(",")
                if not len(cnv) == int(row[2]):
                    print("ERROR")
                dict[sample_id] = cnv
    return dict


def get_cancer_cnvs(tads, cnvs):
    cnv_dict = {}
    for c in cnvs:
        for sequences in cnvs[c]:
            sequences = re.sub("\s+", "", sequences)
            seq = sequences.split(":")
            chr = re.sub("chr", "", seq[0])
            start = int(seq[1])
            end = int(seq[2])

            try:
                tds = tads[chr]

                for t in tds:
                    tad_start = t[0]
                    tad_end = t[1]
                    if start >= tad_start and end <= tad_end:
                        if not c in cnv_dict:
                            cnv_dict[c] = [sequences]
                        else:
                            cnv_dict[c].append(sequences)
            except:
                continue
    
    return cnv_dict


def write_file(data, name):
    f = open(name, "w+")
    for d in data:
        if isinstance(data[d], list):
            f.write("{}\t".format(d))
            i = 0
            for s in data[d]:
                f.write("{}".format(s))
                i += 1
                if i < len(data[d]):
                    f.write(",")
            f.write("\n")
        else:
            f.write("{}\t{}\n".format(d, data[d]))
    f.close()


def main():
    cancer = "output/overlaps/HiCLNCaP_HiCPC3_25000_overlap.bed"
    healthy = "data/HiCPrEC/25000/TADs/HiCPrEC_25000.tad"
    # outfile = "output/cancer_tads_25000.bed"
    outfile = "output/healthy_tads_25000.bed"
    subtract_overlaps(healthy, cancer, outfile)
    tads = get_tads(outfile)
    # mutations, mutation_count = get_mutations()
    # cancer_muts = cancer_mutations(tads, mutations)
    # cancer_muts_count = top_cancer_mutations(mutation_count, cancer_muts)
    # write_file(cancer_muts_count, "output/top_cancer_mutations.tsv")
    # cnvs = get_cnv_samples()
    # cancer_cnvs = get_cancer_cnvs(tads, cnvs)
    # write_file(cancer_cnvs, "output/cancer_cnvs.tsv")    

if __name__ == "__main__":
    main()