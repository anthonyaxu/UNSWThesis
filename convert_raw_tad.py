#!/usr/bin/python3

import statistics as stats

cell_path = "../../MahdiehLabani/HiC_prostate/{}/output/hic_results/matrix/Dixon_2M/raw/{}/{}{}"
cell = [
    "HiCPrEC",
    "HiCLNCaP",
    "HiCPC3"
]

resolution = [
    # "5000",
    # "10000",
    # "25000",
    "40000",
    # "500000"
]

chr = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"
]  

# Filter TADlevel 1 from raw file and write to new file
# for c in cell:
#     for r in resolution:
#         tad_path = "data/{}/{}/TADs/".format(c, r)

#         for ch in chr:
#             try:
#                 read_file = tad_path + "{}_{}_chr{}_tad.tad".format(c, r, ch)
#                 f = open(read_file)
#                 content = f.readlines()
#                 f.close()

#                 write_file = tad_path + "{}_{}_chr{}_filtered.tad".format(c, r, ch)
#                 f = open(write_file, "w+")
#                 for line in content:
#                     tad_level = line.split("\t")[2]
#                     if tad_level == "1":
#                         f.writelines(line)

#                 f.close()
#                 print(write_file + " file written\n")
#             except:
#                 print("File does not exist. Continuing\n")
#                 continue

total_gene_length = 3000000000

# Calculate median, mean TAD sizes for each chromosome and cell
for c in cell:
    for r in resolution:
        tad_path = "data/{}/{}/TADs/".format(c, r)

        res_tad_size = 0
        res_tad_list = []

        # tad_out_name = tad_path + "{}_{}.tad".format(c, r)
        tad_out_name = "test_{}_40000.tad".format(c)
        tad_out = open(tad_out_name, "w+")

        print("-------------------------------------------")
        print(c + "-" + r + ":")
        for ch in chr:
            chr_tad_size = 0
            chr_tad_list = []
            try:
                read_file = tad_path + "{}_{}_chr{}_filtered.tad".format(c, r, ch)
                f = open(read_file)
                content = f.readlines()
                f.close()

                for line in content:
                    print(line)
                    l = line.split("\t")
                    start = int(l[0])
                    end = int(l[1])
                    level = int(l[3])
                    mean = float(l[4])
                    score = float(l[5])
                    diff = (end - start)*int(r)
                    chr_tad_list.append(diff)
                    chr_tad_size += diff
                    tad_out_line = "chr{}\t{}\t{}\n".format(ch, start, end)
                    # tad_out_line = "chr{}\t{}\t{}\t{}\t{}\n".format(ch, start, end, level, mean, score)
                    tad_out.writelines(tad_out_line)
                
                # chr_mean_tad = stats.mean(chr_tad_list)
                # chr_median_tad = stats.median(chr_tad_list)
                # print("CHR{} - NUMBER OF TADs: {}, MEAN: {:.0f}".format(ch, len(chr_tad_list), chr_mean_tad))

                # res_tad_size += chr_tad_size
                # res_tad_list.extend(chr_tad_list)
            except:
                continue
        
        tad_out.close()

        # res_mean_tad = stats.mean(res_tad_list)
        # res_median_tad = stats.median(res_tad_list)
        # print("OVERALL")
        # print("TAD SIZE: {}, # TAD: {}, MEAN: {:.0f}".format(res_tad_size, len(res_tad_list), res_mean_tad))
        # print("TAD COVERAGE: {:.2f}%".format(res_tad_size/total_gene_length*100))