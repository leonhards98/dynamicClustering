# Constructs a random const window. First, argv[2] points are inserted without deletions. Then, argv[4] points are inserted
# as a random sliding window with an insertion probability of argv[3].
# This is followed by 5000 operations with an insertion probability of 0.5 and another 5000 operations with an insertion probability of 0.1
 
# 1. Name of the input file without the .txt
# 2. Length of construction phase (how many points are inserted before any can be deleted)
# 3. Probabilty of a operation after the construction phase being an insertion (0.01 and 0.99)
# 4. Total number of points to insert

import sys
import queue
import re
import random

currentLines = []
windowSize = int(sys.argv[2])
prob = float(sys.argv[3])
infile = sys.argv[1]
totalSize = int(sys.argv[4])
count = 0

outfile = infile+"_randomConst"+str(totalSize)+"_"+str(prob)+".txt"
outfile = open(outfile, "w")
fstream = open(infile+".txt", "r")
while count <= windowSize:
    line = fstream.readline()
    parts = re.split(r' |\t,', line)
    count += 1
    outfile.write(' '.join(map(str,parts)))
    parts[1] = 0
    currentLines.append(' '.join(map(str,parts)))

while len(currentLines) <= totalSize:
    rand = random.uniform(0.01,0.99)
    if rand < prob:
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

for i in range(0,5000,1):
    if (i % 2) == 0:
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

for i in range(0,5000,1):
    rand = random.uniform(0.01,0.99)
    if rand >= prob:
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


outfile.close()
