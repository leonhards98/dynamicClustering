# Constructs a snake window. First, points are insered without deletions until argv[2] points were inserted.
# Then, points are deleted with a probability of argv[3] - 1 until only 0.2 * argv[2] points remain.
# The insertion probability is increase to argv[3], until all argv[2] points are present in the dataset. 
# This pattern is repeated, until all points of the static dataset have been used.

# 1. Name of the input file without the .txt
# 2. Length of construction phase (how many points are inserted before any can be deleted)
# 3. Probabilty of a operation after the construction phase being an insertion (0.01 and 0.99)

import sys
import queue
import re
import random

currentLines = []
windowSize = int(sys.argv[2])
#prob = float(sys.argv[3])
prob = 0.9
infile = sys.argv[1]
count = 0

outfile = infile+"_snake"+str(windowSize)+".txt"
outfile = open(outfile, "w")
fstream = open(infile+".txt", "r")
while count <= windowSize:
    line = fstream.readline()
    parts = re.split(r' |\t,', line)
    count += 1
    outfile.write(' '.join(map(str,parts)))
    parts[1] = 0
    currentLines.append(' '.join(map(str,parts)))

inc = False
while len(currentLines) > 0:
    rand = random.uniform(0.01,0.99)
    if rand < prob:
        if inc:
            line = fstream.readline()
            if line == "":
                break
            parts = re.split(r' |\t,', line)
            outfile.write(' '.join(map(str,parts)))
            parts[1] = 0
            currentLines.append(' '.join(map(str,parts)))
        else:
            rand = random.randint(0,len(currentLines)-1)
            outfile.write(currentLines.pop(rand))
    else:
        if inc:
            rand = random.randint(0,len(currentLines)-1)
            outfile.write(currentLines.pop(rand))
        else:
            line = fstream.readline()
            if line == "":
                break
            parts = re.split(r' |\t,', line)
            outfile.write(' '.join(map(str,parts)))
            parts[1] = 0
            currentLines.append(' '.join(map(str,parts)))

    if inc:
        if len(currentLines) >= windowSize:
            inc = False
    else:
        if len(currentLines) <= windowSize*0.2:
            inc = True

outfile.close()
