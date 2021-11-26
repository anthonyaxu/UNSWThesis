#!/usr/bin/python3

import csv
import sys
import re
from typing import DefaultDict

csv.field_size_limit(sys.maxsize)

RESOLUTION = 40000

data_path = "../data/genes/"

gene_file = data_path + "hg19/genes_hg19.txt"

cnv_tad_file = "output/cancer_tad_cnv_interest.csv"

def get_genes(file):
    dict = {}
    with open(file) as csv_file:
        f = csv.reader(csv_file, delimiter = "\t")
        next(f)

        for row in f:
            chr = row[1]
            start = int(row[2])
            end = int(row[3])
            name = row[7]
            symbol = re.sub("^'", "(", row[9])
            symbol = re.sub("'$", ")", symbol)
            id = name + " " + symbol
            seq = [start, end, id]

            if not chr in dict:
                dict[chr] = [seq]
            else:
                dict[chr].append(seq)
    return dict


def associate(tad_file, genes):
    f = open(tad_file)
    rows = f.readlines()[1:]
    f.close()

    writename = "output/cancer_tad_cnv_interest_genes.tsv"
    writefile = open(writename, "w+")
    writefile.writelines("chr\tstart\tend\tTAD length\t# healthy CNV overlapping\t# cancer CNV overlapping\tFraction for healthy\tFraction for cancer\tEnrichment\tAssociated Genes\tCancer associated genes\n")

    for r in rows:
        r = r.strip()
        r = r.split(",")
        chr = r[0]
        tad_start = int(r[1])
        tad_end = int(r[2])
        length = int(r[3])
        healthy_count = int(r[4])
        cancer_count = int(r[5])
        healthy_fraction = r[6]
        cancer_fraction = r[7]
        enrichment = r[8]

        gene_list = ""
        try:
            for g in genes[chr]:
                gene_start = g[0]
                gene_end = g[1]
                id = g[2]

                if gene_start >= tad_start and gene_end <= tad_end:
                    gene_list += id + ", "
        except:
            pass

        gene_list = gene_list[:-2]     
        writefile.writelines("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr, tad_start, tad_end, length, healthy_count, cancer_count, healthy_fraction, cancer_fraction, enrichment, gene_list, ""))
    writefile.close()


def main():
    genes = get_genes(gene_file)
    associate(cnv_tad_file, genes)


if __name__ == "__main__":
    main()