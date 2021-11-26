#!/usr/bin/python3

import subprocess
import statistics as stats

OnTAD = "../OnTAD/OnTAD"
hic_dir = "../data/hic/"
outdir = "../data/ontad/"

resolution = [
    "40000"
]

cells = [
    "PrEC",
    "LNCaP",
    "PC3"
]

chromosomes = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"
]

total_gene_length = 3088286401

def get_OnTADs():
    try:
        for c in cells:
            for r in resolution:
                chr_path = hic_dir + "{}/{}/".format(c, r)
                tad_path = outdir + "{}/{}/tads/".format(c, r)

                for ch in chromosomes:
                    matrix = chr_path + "{}_{}_chr{}_dense.matrix".format(c, r, ch)
                    out = tad_path + "{}_{}_chr{}_tad".format(c, r, ch)

                    subprocess.run([OnTAD, matrix, "-o", out])
        return True
    except:
        return False


def filter_levels():
    try:
        for c in cells:
            for r in resolution:
                tad_path = outdir + "{}/{}/tads/".format(c, r)

                for ch in chromosomes:
                    try:
                        readfile = tad_path + "{}_{}_chr{}_tad.tad".format(c, r, ch)
                        f = open(readfile)
                        content = f.readlines()
                        f.close()

                        writefile = tad_path + "{}_{}_chr{}_filtered.tad".format(c, r, ch)
                        f = open(writefile, "w+")
                        for line in content:
                            tad_level = line.split("\t")[2]
                            if tad_level == "1":
                                f.writelines(line)
                        f.close()
                    except:
                        continue
        return True
    except:
        return False


def combine():
    try:
        for c in cells:
            for r in resolution:
                tad_path = outdir + "{}/{}/tads/".format(c, r)

                res_tad_size = 0
                res_tad_list = []

                tad_out_name = tad_path + "{}_{}.tad".format(c, r)
                tad_out = open(tad_out_name, "w+")
                tad_out.writelines("chr\tstart\tend\n")

                print("-------------------------------------------")
                print(c + "-" + r + ":")
                for ch in chromosomes:
                    chr_tad_size = 0
                    chr_tad_list = []
                    try:
                        read_file = tad_path + "{}_{}_chr{}_filtered.tad".format(c, r, ch)
                        f = open(read_file)
                        content = f.readlines()
                        f.close()

                        for line in content:
                            l = line.split("\t")
                            start = int(l[0])*int(r)
                            end = int(l[1])*int(r)
                            mean = float(l[3])
                            score = float(l[4])
                            diff = (end - start)
                            chr_tad_list.append(diff)
                            chr_tad_size += diff
                            tad_out_line = "chr{}\t{}\t{}\n".format(ch, start, end)
                            # tad_out_line = "chr{}\t{}\t{}\t{}\t{}\t{}\n".format(ch, start, end, diff, mean, score)
                            tad_out.writelines(tad_out_line)
                        
                        chr_mean_tad = stats.mean(chr_tad_list)
                        chr_median_tad = stats.median(chr_tad_list)
                        print("CHR{} - NUMBER OF TADs: {}, MEAN: {:.0f}".format(ch, len(chr_tad_list), chr_mean_tad))

                        res_tad_size += chr_tad_size
                        res_tad_list.extend(chr_tad_list)
                    except:
                        continue
                
                tad_out.close()
                res_mean_tad = stats.mean(res_tad_list)
                print("OVERALL")
                print("TAD SIZE: {}, # TAD: {}, MEAN: {:.0f}".format(res_tad_size, len(res_tad_list), res_mean_tad))
                print("TAD COVERAGE: {:.2f}%".format(res_tad_size/total_gene_length*100))
        return True
    except:
        return False

def main():
    # print("TADs failed!") if not get_OnTADs() else print("TADs successful!")
    # print("Filtering failed!") if not filter_levels() else print("Filtering successful!")
    print("Combine failed!") if not combine() else print("Combine successful!")


if __name__ == "__main__":
    main()