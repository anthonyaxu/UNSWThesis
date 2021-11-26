#!/usr/bin/python3

import sys, csv, re, glob

csv.field_size_limit(sys.maxsize)

RESOLUTION = 40000
mut_path = "../data/mutations/"

def get_mutations():
    dict = {}
    # count_dict = {}
    with open(mut_path + "top_mutations_5.csv") as csv_file:
    # with open(mut_path + "samples_per_mutation.csv") as csv_file:
            f_m = csv.reader(csv_file, delimiter = ",")
            # Skip csv header
            next(f_m)

            for row in f_m:
                id = row[0]
                seq = row[1].split(":")
                chr = re.sub("chr", "", seq[0])
                if chr != "MT":
                    start = int(seq[1])
                    end = int(seq[2])
                    dict[id] = [chr, start, end]
                    # count = int(row[3])
                    # count_dict[id] = count

    # count_dict = dict(sorted(count_dict.items(), key = lambda x:x[1], reverse = True))
    # return dict, count_dict
    return dict


def get_samples():
    dict = {}
    with open(mut_path + "mutations_per_sample.csv") as csv_file:
            f_s = csv.reader(csv_file, delimiter = ",")
            # Skip csv header
            next(f_s)

            for row in f_s:
                sample_id = row[0]
                mutations = row[1].split(",")
                dict[sample_id] = mutations
    return dict


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

# Total number of mutations appearing in TAD regions
def tad_mutation_count(mutation_dict, tad_dict):
    num_mutations = len(mutation_dict)
    print("Mutations")
    for k in tad_dict:
        print(k, end = " -> ")
        dict = tad_dict[k]
        mutation_sum = 0
        for id in mutation_dict:
            m = mutation_dict[id]
            mut_chr = m[0]
            mut_start = int(m[1])
            mut_end = int(m[2])

            try:
                tds = dict[mut_chr]
                for t in tds:
                    t_start = int(t[0])
                    t_end = int(t[1])
                    if mut_start >= t_start and mut_end <= t_end:
                        mutation_sum += 1
            except KeyError:
                continue
        print("{}, {:.2f}%".format(mutation_sum, mutation_sum/num_mutations*100))


# How many samples had at least one mutation in each TAD region
def tad_sample_count(sample_dict, mutation_dict, tad_dict):
    num_samples = len(sample_dict)
    print("Samples")
    for k in tad_dict:
        print(k, end = " -> ")
        dict = tad_dict[k]
        sample_count = 0

        for id in sample_dict:
            muts = sample_dict[id]
            for m_id in muts:
                if m_id in mutation_dict:
                    m_id = m_id.strip()
                    m = mutation_dict[m_id]
                    mut_chr = m[0]
                    mut_start = int(m[1])
                    mut_end = int(m[2])

                    try:
                        tds = dict[mut_chr]
                        for t in tds:
                            t_start = t[0]
                            t_end = t[1]
                            if mut_start >= t_start and mut_end <= t_end:
                                sample_count += 1
                                continue
                    except KeyError:
                        continue
        print("{}, {:.2f}%".format(sample_count, sample_count/num_samples*100))
        


def main():
    # mutations, mutation_counts = get_mutations()
    mutations = get_mutations()
    # print(mutations)
    samples = get_samples()
    print(len(samples))

    # paths = [
    #     "../data/ontad/PrEC/{}/tads/*".format(RESOLUTION),
    #     "../data/ontad/LNCaP/{}/tads/*".format(RESOLUTION),
    #     "../data/ontad/output/*"
    # ]
    
    # tad_paths = []
    # for p in paths:
    #     for fname in glob.glob(p):
    #         if re.search(".*_(\d+)\.tad", fname) or re.search(".*intersect\.tad", fname):
    #             tad_paths.append(fname)

    # tads = tad_dictionary(tad_paths)
    # tad_mutation_count(mutations, tads)
    # tad_sample_count(samples, mutations, tads)



if __name__ == "__main__":
    main()