# Removes any duplicate points (lines) from the input dataset.

# 1. Input file
# 2. Output file

import sys
import random

lines_seen = {} # holds lines already seen
outfile = open(sys.argv[2], "w")
for line in open(sys.argv[1], "r"):
    #line = line.split("\t")[1]
    if line not in lines_seen:
        lines_seen[line] = 1
    else:
        lines_seen[line] += 1

counter = 0
while len(lines_seen) > 0:
    (line,count) = lines_seen.popitem()
    outfile.write(str(counter) + " 1 "  + str(count) + " " + line)
    counter += 1

outfile.close()
