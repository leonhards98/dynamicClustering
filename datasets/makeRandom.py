# Randomly reorders all points in the given file
# 1. Input File
# 2. Output File


import sys
import random

lines_seen = []
counter = 0
outfile = open(sys.argv[2], "w")
for line in open(sys.argv[1], "r"):
        counter += 1
        lines_seen.append(line)
for i in range(0,counter,1):
    rand = random.randint(0,counter-i-1)
    outfile.write(lines_seen[rand])
    lines_seen.pop(rand)
outfile.close()
