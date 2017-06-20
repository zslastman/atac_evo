"""Reads in file names from stdin
Basically the same as the other program except that the column names are taken from the first input file

find . -name "*perm.results.txt" -print | python ../../src/compile_permutations.py myout
"""

import sys
import os
outFile = sys.argv[1]
toStrip = sys.argv[2] #if anything is given here, the program will use the file name to replace the dataID. What is written here will be stripped from the rihs of the file name. Otherwise give None.
fileNames = sys.stdin



nameVector = None
f_out = open(outFile, 'w')
for myFile in fileNames:
    f = open(myFile.rstrip('\n\r'), 'r')
    myLines = [k.rstrip('\n\r') for k in f if k.rstrip('\r\n') !='']
    if nameVector is None:
        nameVector = myLines[0]
        f_out.write(nameVector + '\n')
    myLines.pop(0)
    if toStrip != 'None' and toStrip != 'none' and toStrip != 'NONE':
        newName = os.path.basename(myFile.rstrip('\n\r')).rstrip(toStrip)
        for i in range(len(myLines)):
            atoms = myLines[i].split('\t')
            atoms[0] = newName
            myLines[i] = '\t'.join(atoms)
    for myLine in myLines:
        f_out.write(myLine + '\n')
f_out.close()

