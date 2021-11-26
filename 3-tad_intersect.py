#!/usr/bin/python3

import subprocess
import os
import re

cells = "PrEC LNCaP"
ontad_path = "../data/ontad/"

RESOLUTION = "40000"

def intersect():
    c = cells.split(" ")
    a = c[0]
    b = c[1]

    a_path = ontad_path + "{}/{}/tads/{}_{}.tad".format(a, RESOLUTION, a, RESOLUTION)
    b_path = ontad_path + "{}/{}/tads/{}_{}.tad".format(b, RESOLUTION, b, RESOLUTION)

    out = ontad_path + "output/ontad_{}_{}_{}_intersect.tad".format(a, b, RESOLUTION)

    subprocess.run("bedtools intersect -a {} -b {} -header -bed > {}".format(a_path, b_path, out), shell = True)


def main():
    intersect()


if __name__ == "__main__":
    main()