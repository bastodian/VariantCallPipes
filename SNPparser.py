#!/usr/bin/env python
'''
    Script to convert vcf file into Structure input file.

    As Infile specify the vcf file from which to call genotypes. MissingData 
    specifies the amount of missing data allowed for each locus.
'''

import sys
from os import system as bash

Infile = sys.argv[1]
CombinedOutfile = sys.argv[2]
MissingRowData = float(sys.argv[3])
MissingColumnData = float(sys.argv[4])
OutfilesA = []
OutfilesB = []

with open(Infile) as File:
    for line in File:
        if '#' in line:
            for individual in line.split('\t')[9:]:
                OutfilesA.append(individual.rstrip('\n') + '_A')
                OutfilesB.append(individual.rstrip('\n') + '_B')
            FileCount = -1
            for FileA in OutfilesA:
                FileCount += 1
                fA = open(OutfilesA[FileCount], 'w')
                fA.write(OutfilesA[FileCount][:-2] + '\t')
                fA.close()
            FileCount = -1
            for FileB in OutfilesB:
                FileCount += 1
                fB = open(OutfilesB[FileCount], 'w')
                fB.write(OutfilesB[FileCount][:-2] + '\t')
                fB.close()
        elif '#' not in line and line.split('\t')[9:].count('./.') < \
                MissingRowData * float(len(line.split('\t')[9:])):
            Genotypes = line.split('\t')[9:]
            FileCount = -1
            for genotype in Genotypes:
                FileCount += 1
                fA = open(OutfilesA[FileCount], 'a')
                fB = open(OutfilesB[FileCount], 'a')
                fA.write(genotype[0:1] + '\t')
                fB.write(genotype[2:3] + '\t')
                fA.close()
                fB.close()
        else:
            continue
               
FileCount = -1
for newline in OutfilesA:
    FileCount += 1
    fA = open(OutfilesA[FileCount], 'r+')       
    for line in fA:
        newline = line.rstrip('\t')
    fA.seek(0)
    fA.write(newline + '\n')
    fA.close()
FileCount = -1
for newline in OutfilesB:
    FileCount += 1
    fB = open(OutfilesB[FileCount], 'r+')
    for line in fB:
        newline = line.rstrip('\t')
    fB.seek(0)
    fB.write(newline + '\n')
    fB.close()

### COMBINE outfiles
FileCount = -1
Outfiles = OutfilesA
for file in OutfilesB:
    Outfiles.append(file)

CombinedGenotypes = open(CombinedOutfile, 'w')
for infile in sorted(Outfiles):
    with open(infile) as Combine:
        for line in Combine:
            if line.count('.') < MissingColumnData * float(len(line.split('\t')[1:])):
                line = line.replace('.', '-9')
                CombinedGenotypes.write(line)
            else:
                continue
        bash('rm ' + infile)
CombinedGenotypes.close()
