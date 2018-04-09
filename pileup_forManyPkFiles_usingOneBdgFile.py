'''
Created on 04/06/2014
Updated on 8/12/2015

@author: Dinesh

Given a list of pkFiles (with their directories and their names) and a pileUp.bdg file, calculate and output the pileUp
sum values in all the peak files for that single bdg file.

Notes:
1. NOTE THAT BOTH THE BDG FILE AND THE PEAK FILE ARE ASSUMED TO HAVE BEEN NAMED IN THIS FORMAT:
    Nhdf.H2az.Rep1.q0.05_treat_pileup.bdg
    cellLine.protein.Rep#.qVal_peaks*
2. chrM peaks are ignored for pileup calculation.


>>>>>> IMPORTANT <<<<<<<
1. If condition is "usePkSummit",IT IS IMPORTANT THAT THE INPUT PEAK FILE IS THE *_PEAKS.ENCODEPEAKS.
2. If the pkFileName has encodePeak, then the last column is assumed to be dist of summit from ss, and the next to last is assumed to be the qVal.
3. The code is not suited for when peakFile is *broadPeak or *gappedPeak and condition is 'usePkSummit'.

=== updated on Nov 17, 2017 ====
Removing this: "Both bdg and pkFiles are assumed to be named in this format: cellLine.[proteinName/someInfo].Rep#.[p_/q_]Val_[peaks/treat_pileup].[txt/bed/narrowPeak/encodePeak[s]/bdg]. For eg.: Nhdf.H2az.Rep1.q0.05_treat_pileup.bdg"
'''

import re
import numpy as np
import collections as c
import argparse

parser = argparse.ArgumentParser(description="Given a list of pkFiles (with their directories and their names) and a pileUp.bdg file, calculate and output the pileUp sum values in all the peak files for that single bdg file. Note: 1. chrM peaks are ignored for pileup calculation. 2. The code is not suited for when peakFile is *broadPeak or *gappedPeak and condition is 'usePkSummit'.")
parser.add_argument("-c", "--condition", default="useWholePk", help="Can be any of 'useWholePk' or 'usePkMid' or 'usePkSummit'; the last is only applicable for when *encodePeaks are used.", choices=['useWholePk', 'usePkMid', 'usePkSummit'], required=True)
parser.add_argument("-z", "--chromsize", default=None, help="Path with filename to chromosome size file.", required=True)
parser.add_argument("-s", "--scan", action='store_true', help="If this option is set, the region of interest is scanned using a window of size -w and step size of -d. IMPORTANT!! --> Currently, the code is written for 'usePkMid' as the condition, and with the assumption that the region is equally partitioned by a fixed number of non-overlapping windows. Another assumption is that the scanning window = 2 * step size. So, a good set of parameters to use is say, 300bp as flank, 200bp as window and 100bp as the stepsize..")
parser.add_argument("-f", "--flank", type=int, help="Number of flanking bps to go on each side of the pk. If 'condition' is 'useWholePk' and flanks are not desired, set this to 0.", default=0, required=True)
parser.add_argument("-w", "--window", type=int, help="If --scan is true, this window size (in bp) used to scan the region of interest around the peak. The region of interest is defined as [pkMid-flank, pkMid+flank] (if condition is 'usePkMid'), for instance, and along this is this window scanned using the step size option (-d) that also must be specified", default=None)
parser.add_argument("-d", "--stepsize", type=int, help="If --scan is true, this step size (in bp) used to move the window of size -w in the region of interest (defined by the flanks from the pkMid (if condition is 'usePkMid') or edges of peak-domain (if condition is 'useWholePk')).", default=None)
parser.add_argument("-p", "--pkFiles", help="Full path (with name) of file with lists of pkFiles (with full path)", required=True)
parser.add_argument("-b", "--bdgFile", help="Full path (with name) of the bdg file to be used.", required=True)
parser.add_argument("-o", "--outputDir", help="Full path of the output directory", required=True)


args = parser.parse_args()
chrSizeFile = args.chromsize
condition = args.condition
flank = args.flank
listOfPkFiles = args.pkFiles
bdgFile = args.bdgFile
outputDir = args.outputDir

if (args.scan):
    windowLen = args.window
    stepSize = args.stepsize
    assert condition == "usePkMid"
    assert (2 * flank) % windowLen == 0  # not checking a/b is an int. since int/int returned a float (while testing)
    assert (2 * stepSize) % stepSize == 0


def getChrSize(chrom):
    handle = open(chrSizeFile, "r")
    size = 0
    for line in handle:
        vals = re.split("\s+", line.strip())
        if (vals[0] == chrom):
            size = int(vals[1])
    handle.close()
    return size


def initArray(chrSize):
    D = np.zeros(chrSize, dtype=float)
    return D


'''
Make either of the two G dictionaries.
G[ (chrom, ss, es) ] = 0 for the sum (or max) pileup value at the peak
G[ (chrom, ss, es) ] = [0, 0,..,0] if we are 'scanning'
'''


def getPeakDict_forAChr(chrom, pkFile):
    handleGenesIn = open(pkFile, "r")
    G = c.OrderedDict()
    chrFound = False
    for line in handleGenesIn:
        vals = re.split("\s+", line.strip())
        ss = int(vals[1])
        es = int(vals[2])

        if (condition == "usePkSummit"):
            assert len(vals) > 6  # bed6 also doesn't have pksummit
            summitLoc = ss + int(vals[-1])  # ss + number of bps off ss that the summit is located
            ss = summitLoc - flank
            es = summitLoc + flank
        elif (condition == "usePkMid"):
            mid = int((ss + es) / 2)
            ss = mid - flank
            es = mid + flank
        elif (condition == "useWholePk"):
            ss = ss - flank
            es = es + flank
        else:
            raise Exception("condition can only be one of 'useWholePk' or 'usePkMid' or 'usePkSummit'")

        if (vals[0] == chrom):
            G = updateG_whileInitializing(G, vals)
            chrFound = True
        if ((chrFound) and (vals[0] != chrom)):
            break
    handleGenesIn.close()
    return G


'''
While initializing the G, update it by adding a new element. Just as a reminder/note,
we are making certain assumptions about the length of window, stepsize and the flanks. (see the argparser notes)
'''


def updateG_whileInitializing(G, vals):
    if (len(vals) < 3):  # can happen with EOF
        return G
    chrom, ss, es = vals[:3]
    if (args.scan):
        numOfScansToDo = ((2 * flank) / windowLen) - 1
        list_val = list(np.zeros(numOfScansToDo, dtype=int))
        G[(chrom, ss, es)] = list_val
    else:
        G[(chrom, ss, es)] = 0
    return G


'''
Update each of the G dicts in Gs to update them with the pileups.
'''


def updateGs_orig(Gs, D):
    for G in Gs:
        for k in G.keys():
            for i in range(0, k[2] - k[1]):
                G[k] += D[k[1] + i]
    return Gs


'''
Update each of the G dicts in Gs to update them with the pileups.
'''


def updateGs(Gs, D):
    for G in Gs:
        for k in G.keys():
            r = k[2] - k[1]  # region of interest; should be divisible by windowLen if scan is used (which also means we are using 'usePkMid')

            if (args.scan):
                pileups = G[k]
                window_ss = k[1]
                for ind in range(0, len(pileups)):
                    for i in range(0, windowLen):
                        pileups[ind] += D[window_ss + i]
                    window_ss = window_ss + stepSize
                G[k] = pileups

            else:
                for i in range(0, r):
                    try:
                        G[k] += D[k[1] + i]
                    except:
                        print("k", k)
                        print("i", i)
                        print("D", D)
                        print("len(D)", len(D))
                        raise Exception()
    return Gs


def getDforNewChr(currChr):
    chrSize = getChrSize(currChr)
    D = initArray(chrSize)
    return D


def getGsForNewChr(currChr):
    handlePks = open(listOfPkFiles, "r")
    Gs = []
    for pkFile in handlePks:
        G = getPeakDict_forAChr(currChr, pkFile.strip())
        Gs.append(G)
    return Gs


def getPileUpNameInfo(pk_or_bdgFileName, isPkFile):
    # assumed to be named this way: cellName.protein.rep.qVal_** where ** can be peaks.encodePeak or peaks.bed or peaks.txt
    # pkFile is a boolean variable
    pk_or_bdgFileName = pk_or_bdgFileName[pk_or_bdgFileName.rfind("/") + 1:]
    vals = re.split("\.", pk_or_bdgFileName)  # ['Fd', 'myoDwExtSize100', 'Rep5', 'p0', '001_treat_pileup', 'bdg']
    cellName = vals[0]
    protein = vals[1]
    rep = vals[2]
    cutOff_prefix = vals[3]
    cutOffInfo = vals[4]
    # fileExtension = vals[5]
    cutOff = cutOff_prefix + "." + re.split("_", cutOffInfo)[0]
    return cellName, protein, rep, cutOff


def writeGs(Gs, outFileNames):
    if (len(Gs) != len(outFileNames)):
        raise Exception()
    for i in range(0, len(Gs)):
        outFileName = outFileNames[i]
        G = Gs[i]
        handleOut = open(outputDir + "/" + outFileName, "a")
        for k, v in G.items():
            if (k[0] == "chrM"):
                continue
            else:
                if (args.scan):
                    pileups = ""
                    for item in v:
                        pileups += str(item) + ","
                    pileups = pileups[:-1]
                    toWrite = k[0] + "\t" + str(k[1]) + "\t" + str(k[2]) + "\t" + pileups + "\n"
                else:
                    if (len(k) == 3):
                        toWrite = k[0] + "\t" + str(k[1]) + "\t" + str(k[2]) + "\t" + str(v) + "\n"
                    else:
                        toWrite = k[0] + "\t" + str(k[1]) + "\t" + str(k[2]) + "\t" + k[3] + "\t" + str(v) + "\n"
                handleOut.write(toWrite)
        handleOut.close()


def getPkFileInfo(pkFile):
    # get the fileName from the directory
    pkFileName = pkFile[pkFile.rindex("/") + 1:]
    return pkFileName


def get_outputFileNames_old():
    outFileNames = []
    bdg_cellName, bdg_protein, bdg_rep, bdg_qVal = getPileUpNameInfo(bdgFile, isPkFile=False)  # bdg_qVal can be bdg_pVal
    handlePks = open(listOfPkFiles, "r")
    for pkFile in handlePks:
        pk_cellName, pk_protein, pk_rep, pk_qVal = getPileUpNameInfo(pkFile.strip(), isPkFile=True)  # pk_qVal can be pk_pVal or pk_idrVal
        outFileName = bdg_cellName + "_" + bdg_protein + "_" + bdg_rep + "_" + bdg_qVal + "Pileups.using" + condition[3:] + ".wFlankBp" + str(flank) + ".for." + pk_cellName + "_" + pk_protein + "Pks_" + pk_rep + "_" + pk_qVal + ".txt"
        outFileNames.append(outFileName)
    handlePks.close()
    return outFileNames


def get_outputFileNames():
    outFileNames = []
    bdgfile_info = ".".join(re.split("\.", bdgFile.strip())[:-1])  # just removing the file suffix
    handlePks = open(listOfPkFiles, "r")
    for apkfile in handlePks:
        pkfile_info = ".".join(re.split("\.", apkfile.strip())[:-1])  # just removing the file suffix
        outFileName = "{}_pileupOn_{}.txt".format(bdgfile_info, pkfile_info)
        outFileNames.append(outFileName)
    handlePks.close()
    return outFileNames


def parseBdg():
    handle = open(bdgFile, "r")
    print("working with bdgFile", bdgFile)
    currChr = ""
    chrFound = False
    outFileNames = get_outputFileNames()
    while 1:
        line = handle.readline()
        if (not line):
            break
        bdgVals = re.split("\s+", line.strip())
        ss, es, pileUp = int(bdgVals[1]), int(bdgVals[2]), float(bdgVals[3])
        if (bdgVals[0] == "chrM"):  # skip chrM
            continue
        if (bdgVals[0] != currChr):  # new chrom encountered
            currChr = bdgVals[0]
            if (not chrFound):  # for the first time
                D = getDforNewChr(currChr)
                Gs = getGsForNewChr(currChr)
                for i in range(0, es - ss):
                    D[ss + i] = pileUp  # update the array for the first line
                chrFound = True
            else:
                Gs = updateGs(Gs, D)
                writeGs(Gs, outFileNames)
                del D
                del Gs
                D = getDforNewChr(currChr)
                Gs = getGsForNewChr(currChr)
                for i in range(0, es - ss):
                    D[ss + i] = pileUp  # update the array for the first line
        else:  # same previous chromosome encountered
            for i in range(0, es - ss):
                try:
                    D[ss + i] = pileUp  # update the array
                except:
                    raise Exception()
    Gs = updateGs(Gs, D)
    writeGs(Gs, outFileNames)
    del D
    del Gs
    handle.close()


def main():
    parseBdg()


main()
