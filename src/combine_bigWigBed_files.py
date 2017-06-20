#!/usr/bin/env python2.7
"""
This program will combine multiple BED format outputs from the program bigWigAverageOverBed
(could be modified in the future for TAB format)

In addition to adding a header, the program will add a final column that reflects the file name from which the row of data were taken.
You can modify the first function of this program directly if you want to manipulate this name in any way.

usage: combine_bigWigBed_files.py inDir outFile
"""

import os
import sys

inDir = sys.argv[1] #will process any file that ends in BED
outFile = sys.argv[2]


def modify_fileName(fileName):
    #should be given a basename file name
    if fileName.find('clust') != -1:
        atoms = fileName.split('.')
        atoms = atoms[0].split('_')
        return atoms[-2] + '_' + atoms[-1]
    else:
        atoms = fileName.split('.')
        return atoms[1]

######################

def main(inDir, outFile):
    fileList = [k for k in os.listdir(inDir) if k.split('.')[-1] == 'bed']
    f_out = open(outFile, 'w')
    headerAtoms = ['chromosome', 'bed_start', 'bed_end', 'feature_name', 'score', 'strand', 'statistic', 'source']
    outLine = '\t'.join(headerAtoms) + '\n'
    f_out.write(outLine)
    for myFile in fileList:
        fileName = modify_fileName(myFile)
        f = open(os.path.join(inDir,myFile), 'r')
        for myLine in f:
            atoms = myLine.rstrip('\r\n').split('\t')
            if len(atoms) == len(headerAtoms) - 1:
                atoms.append(fileName)
                outLine = '\t'.join(atoms) + '\n'
                f_out.write(outLine)
        f.close()
    f_out.close()


main(inDir, outFile)
