#!/usr/bin/python3

import csv
import re
import sys

csv.field_size_limit(sys.maxsize)

# For 25000 resolution
RESOLUTION = 25000

def tad_dictionary(names, file_type):
    tads = {}
    for n in names:
        dict = {}
        if file_type == "overlap":
            file_path = "output/overlaps/{}_{}_overlap.bed".format(n, RESOLUTION)
        else:
            file_path = "data/{}/{}/TADs/{}_{}.tad".format(n, RESOLUTION, n, RESOLUTION)
        
        with open(file_path) as file:
            f = csv.reader(file, delimiter = "\t")

            for row in f:
                chr = re.sub("chr", "", row[0])
                seq = [int(row[1])*RESOLUTION, int(row[2])*RESOLUTION]

                if not chr in dict:
                    dict[chr] = [seq]
                else:
                    dict[chr].append(seq)
        tads[n] = dict
    return tads


def get_mutations():
    dict = {}
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
                dict[id] = [chr, start, end]
    return dict


def get_samples():
    dict = {}
    with open("sample_mutations_count.csv") as csv_file:
            f_s = csv.reader(csv_file, delimiter = ",")
            # Skip csv header
            next(f_s)

            for row in f_s:
                sample_id = row[0]
                mutations = row[1].split(",")
                if not len(mutations) == int(row[2]):
                    print("ERROR")
                dict[sample_id] = mutations
    return dict

def total_mutation_count(names, tad_dict, mutation_dict):
    num_mutations = len(mutation_dict)
    for n in names:
        print(n, end = " -> ")
        dict = tad_dict[n]

        mutation_sum = 0
        for id in mutation_dict:
            mutation_chr = mutation_dict[id][0]
            mutation_start = int(mutation_dict[id][1])
            mutation_end = int(mutation_dict[id][2])

            try:
                tads = dict[mutation_chr]
            
                for t in tads:
                    tad_start = t[0]
                    tad_end = t[1]
                    if mutation_start >= tad_start and mutation_end <= tad_end:
                        mutation_sum += 1
            except:
                continue

        print("{}, {:.2f}%".format(mutation_sum, mutation_sum/num_mutations*100))


def sample_mutation_count(names, tad_dict, sample_dict, mutation_dict):
    num_samples = len(sample_dict)
    for n in names:
        print(n)
        print("Mutations per sample", end = " -> ")
        dict = tad_dict[n]

        sample_count = 0
        total_mutations = 0
        loop_counter = 0
        for id in sample_dict:
            sample_exists = False
            mutation_count = 0
            mutations = sample_dict[id]
            for m_id in mutations:
                m_id = m_id.strip()
                mutation_chr = mutation_dict[m_id][0]
                mutation_start = mutation_dict[m_id][1]
                mutation_end = mutation_dict[m_id][2]

                try:
                    tads = dict[mutation_chr]
                    for t in tads:
                        tad_start = t[0]
                        tad_end = t[1]
                        if mutation_start >= tad_start and mutation_end <= tad_end:
                            mutation_count += 1
                            if not sample_exists:
                                sample_exists = True
                                sample_count += 1
                except:
                    pass
            loop_counter += 1
            total_mutations += mutation_count
            print(mutation_count, end = "")
            if loop_counter < num_samples:
                print("+", end = "")
            else:
                print("=", end = "")
        print(total_mutations)
        print("# samples with at least one mutation in TAD region: {}, {:.2f}%".format(sample_count, sample_count/num_samples*100))


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


def get_unique_cnv():
    cnv = []
    with open("copy_number_unique.csv") as csv_file:
        f_cnv = csv.reader(csv_file, delimiter = ",")
        next(f_cnv)

        for row in f_cnv:
            cnv.append(row)
    return cnv


def unique_cnv_count(names, tad_dict, cnv):
    num_cnv = len(cnv)
    for n in names:
        print(n)
        print("Unique CNVs", end = " -> ")
        
        dict = tad_dict[n]
        cnv_count = 0
        for c in cnv:
            seq = c[0].split(":")
            cnv_chr = re.sub("chr", "", seq[0])
            cnv_start = int(seq[1])
            cnv_end = int(seq[2])

            try:
                tads = dict[cnv_chr]
                for t in tads:
                    tad_start = t[0]
                    tad_end = t[1]

                    if cnv_start >= tad_start and cnv_end <= tad_end:
                        cnv_count += 1
            except:
                pass
        print("{}, {:.2f}%\n".format(cnv_count, cnv_count/num_cnv*100))


def cnv_sample_count(names, tad_dict, sample_dict):
    num_samples = len(sample_dict)
    for n in names:
        print(n)
        print("CNVs per sample", end = " -> ")
        dict = tad_dict[n]

        sample_count = 0
        total_cnvs = 0
        loop_counter = 0
        for id in sample_dict:
            sample_exists = False
            cnv_count = 0
            cnvs = sample_dict[id]
            for cnv in cnvs:
                seq = cnv.split(":")
                cnv_chr = re.sub("chr", "", seq[0]).strip()
                cnv_start = int(seq[1])
                cnv_end = int(seq[2])

                try:
                    tads = dict[cnv_chr]
                    for t in tads:
                        tad_start = t[0]
                        tad_end = t[1]
                        if cnv_start >= tad_start and cnv_end <= tad_end:
                            cnv_count += 1
                            if not sample_exists:
                                sample_exists = True
                                sample_count += 1
                except:
                    pass
            loop_counter += 1
            total_cnvs += cnv_count
            print(cnv_count, end = "")
            if loop_counter < num_samples:
                print("+", end = "")
            else:
                print("=", end = "")
        print(total_cnvs)
        print("# samples with at least one CNV in TAD region: {}, {:.2f}%".format(sample_count, sample_count/num_samples*100))


def main():
    mutation_dict = get_mutations()
    sample_dict = get_samples()
    unique_cnv = get_unique_cnv()
    cnv_samples = get_cnv_samples()

    overlap_names = [
        "HiCLNCaP_HiCPC3",
        "HiCPrEC_HiCLNCaP",
        "HiCPrEC_HiCPC3",
    ]
    tad_names = [
        "HiCLNCaP",
        "HiCPC3",
        "HiCPrEC"
    ]

    print("Total mutations in region:")
    overlap_tads = tad_dictionary(overlap_names, "overlap")
    total_mutation_count(overlap_names, overlap_tads, mutation_dict)
        
    regular_tads = tad_dictionary(tad_names, "regular")
    total_mutation_count(tad_names, regular_tads, mutation_dict)

    print("-----------------------------------")

    print("# Overlapped mutations per sample:")
    sample_mutation_count(overlap_names, overlap_tads, sample_dict, mutation_dict)
    sample_mutation_count(tad_names, regular_tads, sample_dict, mutation_dict)

    print("-----------------------------------")

    print("Unique CNVs")
    unique_cnv_count(overlap_names, overlap_tads, unique_cnv)
    unique_cnv_count(tad_names, regular_tads, unique_cnv)

    cnv_sample_count(overlap_names, overlap_tads, cnv_samples)
    cnv_sample_count(tad_names, regular_tads, cnv_samples)


if __name__ == "__main__":
    main()