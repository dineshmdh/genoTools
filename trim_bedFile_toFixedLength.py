#!usr/local/bin/python
'''
 Created on 9/2/2015

Given a bed3+ format file, generate a bed3+ format which has the pkLengths fixed in relation to the pkMid. 
i.e Fix ss and es as 100bp (for instance) away from the pkMid. The remaining fields are kept constant. 
'''

import re
import argparse

parser = argparse.ArgumentParser(description="Given a bed3+ format file, generate a bed3+ format which has the pkLengths fixed in relation to the pkMid. i.e Fix ss and es as 100bp (for instance) away from the pkMid. The remaining fields are kept constant.")
parser.add_argument("-d", "--distance", default=100, help="Number of bases away from the pkMid to go on either end to fix ss and es.", type=int, required=True)
parser.add_argument("-f", "--inputFile", help="Input file (with full path)", required = True)
parser.add_argument("-o", "--outputDir", help="Full path of the output directory", required=True)
parser.add_argument("-t", "--hasHeader", help="The input file has a header. Use 0 if no header, 1 if there is one.", default=0, choices=[0,1], type=int, required=True)

args = parser.parse_args()
dist = args.distance
inputFile = args.inputFile
outputDir = args.outputDir
hasHeader = args.hasHeader

inputFileName = inputFile[(inputFile.rfind("/")+1):]
outputFileName = inputFileName[:-4]+"."+str(dist)+"bpFromMid.txt"

handleIn = open(inputFile, "r")
handleOut = open(outputDir+"/"+outputFileName, "w")
lineCount = 0
for line in handleIn:
    lineCount += 1
    if ( (hasHeader == 1) and (lineCount == 1) ):
	handleOut.write(line.strip()+"\n")
    else:
	vals = re.split("\s+", line.strip())
	ss, es = int(vals[1]), int(vals[2])
	mid = int( (ss+es)/2.0 )
	ss = mid - dist
	es = mid + dist
	
	toWrite = vals[0]+"\t"+str(ss)+"\t"+str(es)
	if (len(vals) > 3):
	    toWrite += "\t"+"\t".join(vals[3:])
	handleOut.write(toWrite.strip()+"\n")

handleOut.close()
handleIn.close()


