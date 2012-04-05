#!/usr/bin/env python

from sys import argv
from random import shuffle

infile = argv[1]
samplesize = float(argv[2])

with open(infile) as Infile:
    count = 0
    for line in Infile:
        while count < 1:
            count += 1
            Locs = range(1, len(line.split()[1:]) + 1)
            shuffle(Locs)
            jacknife = int(samplesize * len(Locs))
            start = 0 - jacknife
            end = 0
            numfiles = int(len(Locs)) / int(jacknife)

for i in range(numfiles):
    outfile = '%s_%s.txt' % (infile.split('.')[0], str(i))
    print '\nwriting %s' % (outfile)
    with open(outfile, 'w') as Outfile:
        start += jacknife
        end += jacknife
        with open(infile) as Infile:
            for line in Infile:
                specimen = line.rstrip('\n').split()
                Outfile.write(specimen[0])
                for locus in Locs[start:end]:
                    Outfile.write('\t %s' % (specimen[locus]))
                Outfile.write('\n')

print '\nDONE\n'
