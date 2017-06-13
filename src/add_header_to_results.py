"""Reads in file names from stdin for INSIGHT results files.
Writes the results as a table


find . -name "*.ins.results.txt" -print | grep -v joint | grep -v permutations | grep -v not | python ../../src/compile_results.py myout
"""

import sys
outFile = sys.argv[1]
fileNames = sys.stdin


myLines = []
for myFile in fileNames:
    f = open(myFile.rstrip('\n'), 'r')
    myLines = myLines + [k.rstrip('\r\n') for k in f if k.rstrip('\r\n') != '']
    f.close()

nameVector=["dataID","thres","rho","rho_stderr","E.A.",
              "E.A._stderr","E.W.","E.W._stderr","alpha","alpha_stderr",
              "tau","tau_stderr","eta","eta_stderr","gamma","gamma_stderr",
              "lnLd","LRT.rho.0.","LRT.eta.0.","LRT.gamma.0.","em_status",
              "anal_sites","anal_blocks"]
minLen = len(nameVector) - 2
f = open(outFile, 'w')
f.write('\t'.join(nameVector) + '\n')
for myLine in myLines:
    atoms = myLine.split('\t')
    if len(atoms) < minLen:
        continue
    atoms[1] = str(float(atoms[1].split('=')[-1].rstrip('%'))/100)
    outLine = '\t'.join(atoms)
    f.write(outLine + '\n')
f.close()
