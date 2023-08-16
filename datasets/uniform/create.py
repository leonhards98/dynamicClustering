import sys
import math

size = sys.argv[1]
outFile = open("uniform"+size+"_insertAll.txt","w")

width = math.ceil(math.sqrt(int(size)))

for i in range(0,int(size),1):
    outFile.write(str(i)+" 1 1 "+str(math.floor(i/width))+" "+str(i%width)+"\n")


