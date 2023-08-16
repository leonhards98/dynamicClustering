# Constructs a random sliding window using all points of the input file

#1. Name of the input file without the .txt
#2. Length of construction phase (how many points are inserted before any can be deleted)
#3. Probabilty of a operation after the construction phase being an insertion (0.01 and 0.99)

import sys
import queue
import re
import random

currentLines = []
windowSize = int(sys.argv[2])
prob = float(sys.argv[3])
infile = sys.argv[1]
count = 0

outfile = infile+"_random"+str(windowSize)+"_"+str(prob)+".txt"
outfile = open(outfile, "w")
fstream = open(infile+".txt", "r")
while count <= windowSize:
    line = fstream.readline()
    parts = re.split(r' |\t,', line)
    count += 1
    outfile.write(' '.join(map(str,parts)))
    parts[1] = 0
    currentLines.append(' '.join(map(str,parts)))

while len(currentLines) > 0:
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
        #rand = random.randint(0,len(currentLines)-1)
        outfile.write(currentLines.pop(0))

outfile.close()
