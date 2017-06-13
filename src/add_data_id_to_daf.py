#!/usr/bin/env python2.7
"""
Just a quick script that adds the basename of the file to a DAF result.
To be safe, always use full path name
"""

import sys
import os


inFile = sys.argv[1]
fileName = os.path.basename(inFile)
f = open(inFile, 'r')
myLines = [k.rstrip('\r\n') for k in f if k.rstrip('\r\n') != '']
f.close()
headerAtoms = myLines.pop(0).split('\t')
if headerAtoms[0] != 'dataID':
    headerAtoms = ['dataID'] + headerAtoms
    f = open(inFile, 'w')
    outLine = '\t'.join(headerAtoms) + '\n'
    f.write(outLine)
    for myLine in myLines:
        atoms = [fileName] + myLine.split('\t')
        outLine = '\t'.join(atoms) + '\n'
        f.write(outLine)
    f.close()

