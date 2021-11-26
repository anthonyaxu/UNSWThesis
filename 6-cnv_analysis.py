#!/usr/bin/python3

import sys, csv, re, glob

csv.field_size_limit(sys.maxsize)

RESOLUTION = 40000
cnv_path = "../data/cnv/"


def tad_dictionary(file_path):
    tads = {}
    for fp in file_path:
        dict = {}
        with open(fp) as file:
            f = csv.reader(file, delimiter = "\t")
            next(f)

            for row in f:
                chr = re.sub("chr", "", row[0])
                seq = [int(row[1]), int(row[2])]

                if not chr in dict:
                    dict[chr] = [seq]
                else:
                    dict[chr].append(seq)
        tads[fp] = dict
    return tads


def healthy_cnv_locations(healthy_dataset):
    dict = {}
    count = 0
    with open(healthy_dataset) as csv_file:
        f = csv.reader(csv_file, delimiter = "\t")
        next(f)

        for row in f:
            count += 1
            chr = row[1]
            start = int(row[2])
            end = int(row[3])
            id = row[4]
            seq = [start, end, id]

            if not chr in dict:
                dict[chr] = [seq]
            else:
                dict[chr].append(seq)
    return dict, count


def cancer_cnv_locations(cancer_dataset):
    dict = {}
    count = 0
    with open(cancer_dataset) as csv_file:
        f = csv.reader(csv_file, delimiter = ",")
        next(f)

        for row in f:
            count += 1
            id = row[3]
            chr = row[11]
            start = int(row[12])
            end = int(row[13])
            seq = [start, end, id]

            if not chr in dict:
                dict[chr] = [seq]
            else:
                dict[chr].append(seq)
    return dict, count

def get_tads():
    path = "../data/ontad/output/ontad_healthy_40000.tad"
    # path = "../data/ontad/output/ontad_cancer_40000.tad"
    f = open(path)
    next(f)
    rows = f.readlines()
    f.close()

    data = []
    for r in rows:
        r = r.strip()
        data.append(r)
    return data


def calculate_overlaps(x, y):
    return range(max(x[0], y[0]), min(x[-1], y[-1]) + 1)


def count_overlaps(tads, healthy, cancer):
    # outfile = "../data/ontad/output/healthy_tad_cnvs.tsv"
    outfile = "../data/ontad/output/cancer_tad_cnvs.tsv"
    writefile = open(outfile, "w+")
    writefile.writelines("chr\tstart\tend\t# healthy CNV overlapping\t# cancer CNV overlapping\n")
    for t in tads:
        t = t.split("\t")
        tad_chr = re.sub("chr", "", t[0])
        tad_start = int(t[1])
        tad_end = int(t[2])
        tad_range = range(tad_start, tad_end)

        healthy_count = 0
        cancer_count = 0

        for s in healthy[tad_chr]:
            healthy_start = s[0]
            healthy_end = s[1]
            healthy_range = range(healthy_start, healthy_end)
            healthy_length = len(healthy_range)

            if healthy_length == 0:
                healthy_length = 1
                if healthy_start >= tad_start and healthy_end <= tad_end:
                    overlap = 1
                else:
                    overlap = 0
            else:
                overlap_range = calculate_overlaps(tad_range, healthy_range)
            overlap = len(overlap_range)

            # Get percentage of CNV located in TAD region
            cnv_coverage = overlap/healthy_length*100

            if cnv_coverage >= 50:
                healthy_count += 1
        
        for s in cancer[tad_chr]:
            cancer_start = s[0]
            cancer_end = s[1]
            cancer_range = range(cancer_start, cancer_end)
            cancer_length = len(cancer_range)

            if cancer_length == 0:
                cancer_length = 1
                if cancer_start >= tad_start and cancer_end <= tad_end:
                    overlap = 1
                else:
                    overlap = 0
            else:
                overlap_range = calculate_overlaps(tad_range, cancer_range)
                overlap = len(overlap_range)

            # Get percentage of CNV located in TAD region
            cnv_coverage = overlap/cancer_length*100

            if cnv_coverage >= 50:
                cancer_count += 1
           
        writefile.writelines("{}\t{}\t{}\t{}\t{}\n".format(tad_chr, tad_start, tad_end, healthy_count, cancer_count))
    writefile.close()


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

    writename = "../data/ontad/output/cancer_tad_cnvs_interest_genes.tsv"
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


def get_cancer_cnv_tads(tads, cancer):
    outfile = "../data/ontad/output/overlap_cancercnvs_healthytads.tsv"
    # outfile = "../data/ontad/output/overlap_healthycnvs_healthytads.tsv"
    writefile = open(outfile, "w+")
    writefile.writelines('track name="overlap cancer cnvs" description="cancer cnvs that overlap with cancer TADs" visibility=3 itemRgb="On"\n')

    counter = 1
    for t in tads:
        t = t.split("\t")
        tad_chr = re.sub("chr", "", t[0])
        tad_start = int(t[1])
        tad_end = int(t[2])
        tad_range = range(tad_start, tad_end)

        for s in cancer[tad_chr]:
            cancer_start = s[0]
            cancer_end = s[1]
            cancer_range = range(cancer_start, cancer_end)
            cancer_length = len(cancer_range)

            if cancer_length == 0:
                cancer_length = 1
                if cancer_start >= tad_start and cancer_end <= tad_end:
                    overlap = 1
                else:
                    overlap = 0
            else:
                overlap_range = calculate_overlaps(tad_range, cancer_range)
                overlap = len(overlap_range)

            # Get percentage of CNV located in TAD region
            cnv_coverage = overlap/cancer_length*100

            if cnv_coverage >= 50:
                writefile.writelines("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(tad_chr, cancer_start, cancer_end, counter, "0", ".", cancer_start, cancer_end, "15,206,162"))
                counter += 1
    writefile.close()


def main():
    healthy_cnv_dataset = cnv_path + "CNV_dataset_healthy.txt"
    cancer_cnv_dataset = cnv_path + "copy_number_somatic_mutations.csv"
    tads = get_tads()
    healthy_cnvs, total_healthy_cnvs = healthy_cnv_locations(healthy_cnv_dataset)  
    cancer_cnvs, total_cancer_cnvs = cancer_cnv_locations(cancer_cnv_dataset)
    # count_overlaps(tads, healthy_cnvs, cancer_cnvs)

    # gene_file = "../data/genes/hg19/genes_hg19.txt"
    # genes = get_genes(gene_file)

    # cnv_tad_file = "../data/ontad/output/cancer_tad_cnvs_interest.csv"
    # associate(cnv_tad_file, genes)

    get_cancer_cnv_tads(tads, cancer_cnvs)


if __name__ == "__main__":
    main()