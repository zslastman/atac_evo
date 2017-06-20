#!/usr/bin/env python2.7
"""Will read in from stdin a set of sites in INSIGHT input format.
The result, written to stdout, will be the counts of polymorphic sites in two categories.
1) Rare derived allele (e.g. < 1.5%)
2) Common derived allele (e.g. > 5%)
"""

import argparse
import sys







def main(args):
    f_in = sys.stdin
    f_out = sys.stdout
    commonCount = 0
    rareCount = 0
    totalCount = 0
    nonAncestral = 0
    totalSites = 0
    totalPolymorphisms = 0
    totalFixed = 0
    for myLine in f_in:
        atoms = myLine.rstrip('\n').split('\t')
        if atoms[0] == 'block' or atoms[0].split(' ')[0] == 'samples':
            continue
        if atoms[0] == 'site':
            siteType = 2
            ancProb1 = 3
            ancProb2 = 4
            alleleCount1 = 5
            alleleCount2 = 6
        else:
            siteType = 5
            ancProb1 = 6
            ancProb2 = 7
            alleleCount1 = 8
            alleleCount2 = 9
        totalSites = totalSites + 1
        if atoms[siteType] == 'M':
            if float(atoms[ancProb1]) <= args.fixed_prob:
                totalFixed = totalFixed + 1
        else:
            totalPolymorphisms = totalPolymorphisms + 1
            #first off, we need to make sure that our probabilities are at least the ancestral minmum
            ancestralProbs = [float(atoms[ancProb1]), float(atoms[ancProb2])]
            totalCount = totalCount + 1
            if max(ancestralProbs) < args.min_ancestral_prob:
                nonAncestral = nonAncestral + 1
                continue
            alleleCounts = [int(atoms[alleleCount1]), int(atoms[alleleCount2])]
            if alleleCounts[1] > alleleCounts[0]:
                raise ValueError('Minor allele count is greater than major allele count!')
            daf = float(alleleCounts[ancestralProbs.index(min(ancestralProbs))])/sum(alleleCounts)
            if daf <= args.rare_threshold:
                rareCount = rareCount + 1
            elif daf >= args.common_threshold:
                commonCount = commonCount + 1
    outHeader = 'rareCount\tcommonCount\ttotalPolymorphisms\tnonAncestralPolymorphisms\tfixedDivergence\ttotalSites\n'
    f_out.write(outHeader)
    outLine = str(rareCount) + '\t' + str(commonCount) + '\t' + str(totalCount) + '\t' + str(nonAncestral) + '\t' + str(totalFixed) + '\t' + str(totalSites) + '\n'
    f_out.write(outLine)











parser = argparse.ArgumentParser(description='Get counts of sites with rare and common derived polymorphic alleles. The input (in INSIGHT format) comes from stdin. The output is printed to stdout')
parser.add_argument('--rare_threshold', '-r', required = False, default = 0.015, type = float, help = 'maximum threshold to be considered a rare allele [default 0.015]')
parser.add_argument('--common_threshold', '-c', required = False, default = 0.05, type = float, help = 'minimum threshold to be considered a common allele [default 0.05]')
parser.add_argument('--min_ancestral_prob', '-m', required = False, default = 0.5, type = float, help = 'minimum posterior probability that the that one of the segregating alleles is the ancestral state. If neither allele can be considered ancestral above this threshold, the site is skipped for the DAF analysis.')
parser.add_argument('--fixed_prob', '-f', required = False, default = 0.5, type = float, help = 'A fixed must have less than this probability of being ancestral to be considered a derived allele')




args = parser.parse_args()
main(args)
