#!/usr/bin/python3

file = "../data/ontad/PrEC/40000/tads/PrEC_40000.tad"

f = open(file, "r")
content = f.readlines()[1:]
f.close()

i = 0
x = True
while i < len(content):
    r = content[i].strip()
    r = r.split("\t")
    chr1 = r[0]
    start1 = int(r[1])
    end1 = int(r[2])

    j = i + 1
    if j < len(content):  
        while True:
            t = content[j].strip()
            t = t.split("\t")
            chr2 = t[0]
            start2 = int(t[1])
            end2 = int(t[2])

            if chr1 == chr2 and end1 == start2:
                j += 1
                end1 = end2
                if j >= len(content):
                    end = end2
                    x = False
                    break
            else:
                i = j - 1
                end = end1
                break

        print("{}\t{}\t{}".format(chr1, start1, end))
    else:
        if x:
            print("{}\t{}\t{}".format(chr1, start1, end1))
    i += 1