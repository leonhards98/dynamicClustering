import sys

inFile = open("birch2.txt", "r")
outFile = open("birch2Skip.txt", "w")
skip = 3
clusterSize = 1000

lines = inFile.readlines()

counter = 0
counter2 = 0
for i,line in enumerate(lines):
    tokens = line.split(" ")
    if tokens[0] == '':
        tokens.remove('')
    if i >= counter*skip*clusterSize:
        if i < counter*skip*clusterSize+clusterSize:
            outFile.write(str(counter2) + " 1 1 " + line)
            counter2 += 1
    if  i == counter*skip*clusterSize+clusterSize:
        counter += 1
