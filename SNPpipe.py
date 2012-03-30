#!/usr/bin/env python
'''
This script maps paired-end Illumina reads to a reference genome using bowtie2. It then converts
the SAM output from Bowtie mapping to BAM. The BAM files are then sorted and indexed and an
identifier is added, identifying the specimens the reads were generated from. All individual
BAM files are then merged, mis-aligned regions are identified and re-aligned using GATK's
re-aligner (i.e., a full Smith-Waterman re-alignment is performed). Lastly, SNPs are called 
from the re-aligned BAM file using GATK's unified genotyper.

Java version => 1.6 required ($ java -version).

You can fine tune Java's memory handling by adding -XX:MaxPermSize=4096m

Input: 1: directory containing input fastq files (name must contain R1 for mate1, R2 for 
          mate2 before file extension; use .fastq file extension)
       2: output directory 
       3: Reference genome in fasta format (extension must be .fa) 
       4: number of threads allocated 
       5: Maximum Java heap size in GB (initial heap size 1GB)

author: Bastian Bentlage
email: bastodian@gmail.com
license: Creative Commons Attribution
'''

import os, sys, glob

inDir = os.path.abspath(sys.argv[1])
outDir = os.path.abspath(sys.argv[2])
ReferenceGenome = os.path.abspath(sys.argv[3]).split('.')[0]
Threads = sys.argv[4]
JavaMem = sys.argv[5]

### Dependencies -- edit these to the paths on your system!

bowtie2 = '/home/bastodian/SNP/Tools/bowtie2-2.0.0-beta5/bowtie2'
samtools ='/home/bastodian/SNP/Tools/SAMtools/samtools-0.1.18/samtools'
picardAddOrReplaceGroups = '/home/bastodian/SNP/Tools/Picard/picard-tools-1.63/AddOrReplaceReadGroups.jar'
picardMergeSamFiles = '/home/bastodian/SNP/Tools/Picard/picard-tools-1.63/MergeSamFiles.jar'
GATK = '/home/bastodian/SNP/Tools/GATK/GenomeAnalysisTK-1.4-37-g0b29d54/GenomeAnalysisTK.jar'

with open(os.path.join(outDir,'log.txt'), 'w') as log:

    ### Check that all the required software is where I expect it

    if (os.path.exists(bowtie2) is False) or (os.path.exists(samtools) is False) \
                    or (os.path.exists(picardAddOrReplaceGroups) is False) or \
                    (os.path.exists(picardMergeSamFiles) is False) or (os.path.exists(GATK) is False):
        log.write('One or more required programs missing\nBowtie2 present? %s\nSamTools present? %s\n\
                    Picard AddOrReplaceGroups present? %s\nPicard MergeSamFiles present? %s\nGATK present? %s' \
                    % (os.path.exists(bowtie2), os.path.exists(samtools), os.path.exists(picardAddOrReplaceGroups), \
                    os.path.exists(picardMergeSamFiles), os.path.exists(GATK)))
        sys.exit()
    
    ### If software is present create pairs of files and run bowtie2 against the reference

    else:
        fastqFileList = sorted(os.listdir(inDir))
        # Make sure there are files in the inDirectory        
        if len(fastqFileList) == 0:
            log.write('Input directory is empty!')
            sys.exit()
        # Make sure there are fastq files in the inDirectory
        elif 'fastq' not in " ".join(fastqFileList):
            log.write('Did not find fastq in file names!')
            sys.exit()
        else:
            pos = 0
            for filename in fastqFileList:
                pos+=1
                myfile = filename.split('_')
                if 'R1' in myfile[-1]:
                    file1 = '%s/%s' % (inDir, filename)
                    SamFileName = '%s.sam' % ('_'.join(myfile[0:3]))
                    BowtieOutFile = os.path.join(outDir, SamFileName)
                    match = '%s' % ('_'.join(myfile[1:3]))
                    for anotherfile in fastqFileList[pos:]:
                        if match in anotherfile:
                            file2 = '%s/%s' % (inDir, anotherfile)
                            log.write('\nBowtie input file 1: %s\nBowtie input file 2: %s\n\nRunning \
                                        Bowtie with %s threads\nReference: %s\nOutput: %s\n' \
                                        % (file1, file2, Threads, ReferenceGenome, BowtieOutFile))
                            bowtieRun = ('%s -p %s -x %s -1 %s -2 %s --no-overlap --no-contain \
                                            --sensitive-local --phred33 -S %s 2>&1' % (bowtie2, \
                                            Threads, ReferenceGenome, file1, file2, BowtieOutFile))
                            os.system(bowtieRun)
                            break
                else:
                    continue
            log.write('\nAll Bowtie runs completed\n')

            ### Convert SAM files to BAM files using SamTools
            
            if len(glob.glob(os.path.join(outDir, '*sam'))) != 0:
                for samfile in glob.glob(os.path.join(outDir, '*sam')):
                    bamfile = '%s.bam' % (samfile.split('.')[0])
                    log.write('\nConverting %s to %s\n\n' % (samfile, bamfile))
                    SamConvert = os.popen('%s import %s.fai %s %s 2>&1' % (samtools, ReferenceGenome, samfile, bamfile)) 
                    log.write("".join(SamConvert.readlines()))
                    if os.path.exists(bamfile) is True:
                        rmSAMfile = 'rm %s' % (samfile)
                        os.system(rmSAMfile)
                    else:
                        log.write('\n error writing %s' % (bamfile))
                        sys.exit()
            else:
                log.write('\nSAM files missing!')
                sys.exit()

            ### Add information to BAM file and sort BAM and index BAm for Realigner

            if len(glob.glob(os.path.join(outDir, '*bam'))) != 0:
                for bamfile in glob.glob(os.path.join(outDir, '*.bam')):
                    ID = os.path.splitext(os.path.split(bamfile)[-1])[0]
                    outbamfile = '%s.sorted.fixed.bam' % (os.path.splitext(bamfile)[0])
                    log.write('\nSorting and indexing: %s\nCreating index: %s\nWriting to: %s\n\n' % (bamfile, ID, outbamfile)) 
                    AddOrReplace = os.popen('java -Xms1g -Xmx%sg -jar %s I=%s O=%s SORT_ORDER=coordinate RGPL=Illumina RGID=%s \
                                            RGLB=%s RGPU=%s RGSM=%s CREATE_INDEX=True -nt %s 2>&1' \
                                            % (JavaMem, picardAddOrReplaceGroups, bamfile, outbamfile, ID, ID, ID, ID, Threads))
                    log.write("".join(AddOrReplace.readlines()))
            else:
                log.write('\nBAM files missing!\nIndexing and sorting failed.')
                sys.exit()

            ### Merge individual BAM files using Picard

            if len(glob.glob(os.path.join(outDir, '*.sorted.fixed.bam'))) != 0:
                bamlist = []
                for sortedFixedBamfile in glob.glob(os.path.join(outDir, '*.sorted.fixed.bam')):
                    bamlist.append(sortedFixedBamfile)
                for infile in range(len(bamlist)):
                    bamlist[infile] = 'I=%s' % (bamlist[infile])
                I = ' '.join(bamlist)
                mergedBAM = '%s.merged.bam' % (os.path.join(outDir, myfile[0]))
                log.write('\nFiles to be merged: %s\nOtputfile: %s' % (I, mergedBAM))
                PicardMerge = os.popen('java -Xms1g -Xmx%sg -jar %s USE_THREADING=yes %s O=%s \
                                MERGE_SEQUENCE_DICTIONARIES=yes -nt %s 2>&1' \
                                % (JavaMem, picardMergeSamFiles, I, mergedBAM, Threads)) 
                log.write("".join(PicardMerge.readlines()))
            else:
                log.write('\nSorted, indexed BAM files missing!\nMerging failed.')
                sys.exit()

            ### After merging insert index using samtools index

            if os.path.exists(mergedBAM) is True:
                # insert index using SAMtools
                log.write('\nIndexing %s\n' % (mergedBAM))
                indexMerge = os.popen('%s index %s 2>&1' % (samtools, mergedBAM))
                log.write("".join(indexMerge.readlines()))
                # Create realignment targets using GATK
                realignTargets =  '%s.realign.targets.list' % (os.path.join(outDir, myfile[0]))
                log.write('\nTarget creation for realignment.\nTargets: %s\n' % (realignTargets))
                RealignerTarget = os.popen('java -Xms1g -Xmx%sg -jar %s -T RealignerTargetCreator -I %s \
                                            -R %s.fa -o %s  -nt %s 2>&1' \
                                            % (JavaMem, GATK, mergedBAM, ReferenceGenome, realignTargets, Threads))
                log.write("".join(RealignerTarget.readlines()))
            else:
                log.write('\nMerged BAM file for indexing is missing!!!\n')
                sys.exit()

            ### Realign using Smith-Waterman algorithm, as implemented in GATK's IndelRealigner

            if os.path.exists(realignTargets) is True:
                realignedBAM = '%s.realigned.bam' % (os.path.join(outDir, myfile[0]))
                log.write('\nRealigning %s\nOutput %s\n\n' % (mergedBAM, realignedBAM))
                Realign = os.popen('java -Xms1g -Xmx%sg -jar %s -I %s -R %s.fa -T IndelRealigner \
                                    -targetIntervals %s -o %s 2>&1' \
                                    % (JavaMem, GATK, mergedBAM, ReferenceGenome, realignTargets, realignedBAM))
                log.write("".join(Realign.readlines()))
            else:
                log.write('\nTarget for alignment missing!\n')
                sys.exit()

            ### Call SNPs using GATK's Unified Genotyper
            
            if os.path.exists(realignedBAM) is True:
                SNPvcfFile = '%s.vcf' % (os.path.join(outDir, myfile[0]))
                log.write('\nCalling SNPs from %s\nOutput %s\n\n' % (realignedBAM, SNPvcfFile))
                SNPcall = os.popen('java -Xms1g -Xmx%sg -jar %s -R %s.fa -T UnifiedGenotyper -I %s -o %s \
                                    -stand_call_conf 50.0 -stand_emit_conf 30.0 -nt %s 2>&1' \
                                    % (JavaMem, GATK, ReferenceGenome, realignedBAM, SNPvcfFile, Threads))
                log.write("".join(SNPcall.readlines()))
            else:
                log.write('\nRealigned file missing!\n')
                sys.exit()
