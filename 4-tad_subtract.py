#!/usr/bin/python3

import subprocess
import os
import re

RESOLUTION = "40000"

def subtract(a, b, outfile):
    subprocess.run("bedtools subtract -a {} -b {} -header -bed > {}".format(a, b, outfile), shell = True)


def main():
    cancer = "../data/ontad/LNCaP/{}/tads/LNCaP_{}.tad".format(RESOLUTION, RESOLUTION)
    healthy = "../data/ontad/PrEC/{}/tads/PrEC_{}.tad".format(RESOLUTION, RESOLUTION)
    subtract(healthy, cancer, "../data/ontad/output/ontad_healthy_{}.tad".format(RESOLUTION))
    subtract(cancer, healthy, "../data/ontad/output/ontad_cancer_{}.tad".format(RESOLUTION))

if __name__ == "__main__":
    main()