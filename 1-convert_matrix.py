#!/usr/bin/python3

import subprocess
import os.path

sparseToDense = "utils/sparseToDense.py"
splitSparse = "utils/split_sparse.py"
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

def split2chromosome():
    try:
        for c in cells:
            for r in resolution:
                matrix = hic_dir + "{}/{}/dixon_2M_{}.matrix".format(c, r, r)
                bed = hic_dir + "{}/{}/dixon_2M_{}_abs.bed".format(c, r, r)
                output = hic_dir + "{}/{}/{}_{}".format(c, r, c, r)

                subprocess.run(["python3", splitSparse, "-b", bed, matrix, "-o", output])
        return True
    except:
        return False


def sparse2dense():
    try:
        for c in cells:
            for r in resolution:
                path = hic_dir + "{}/{}/".format(c, r)
                for ch in chromosomes:
                    matrix = path + "{}_{}_chr{}.matrix".format(c, r, ch)
                    bed = path + "{}_{}_chr{}_abs.bed".format(c, r, ch)
                    out = path + "{}_{}_chr{}_dense.matrix".format(c, r, ch)

                    subprocess.run(["python3", sparseToDense, "-b", bed, matrix, "-o", out])
        return True
    except:
        return False

def main():
    print("Splitting chromosome failed!") if not split2chromosome() else print("Splitting chromosome successful!")
    print("Sparse2Dense failed!") if not sparse2dense() else print("Sparse2Dense successful!")


if __name__ == "__main__":
    main()
