#!/usr/bin/env python 

###############################################################################
#   jpdFilter5.py by Jessilyn Dunn on 2015-10-08
#   Dependencies: python/2.7
#      
#   Use: module load python/2.7;  ./jpdFilter7.py -f inputfile -o output file -r (tab-delimited list of cases) -c (tab-delimited list of ctrls)
#   Batch use: run with betaThalSTMP.sh wrapper
#   
###############################################################################

###############
#Import Modules
from __future__ import print_function
import argparse
import re
import operator
from collections import OrderedDict
import os


###################
# Argument parsing

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

descr = """filter a VCF or STMP file based on user-defined criteria to generate list of candidate genes."""
Arguments = argparse.ArgumentParser(description=descr, add_help=True)

Arguments.add_argument('-f',
	'--filename',
    ### check that filename ends in .vcf or .tsv. if not, print "Filename does not end in .vcf or .tsv"
	help='VCF or STMP file (.vcf or .tsv)')
Arguments.add_argument('-o',
		'--outfile', type=str,
		help='output filename', default='outfile.genvar')
Arguments.add_argument('-t',
		'--totals', type=str,
		help='total gene SNP count filename', default='totals.genvar')
Arguments.add_argument('--cases',
        nargs='*', type=str, help='list of case sample names')
Arguments.add_argument('--controls',
		nargs='*', type=str, help='list of control sample names')
Arguments.add_argument('--mtype',
		nargs='+', type=str, help='list of types of mutations; default is for potentially deleterious mutations: nonsynonymous SNV, frameshift deletion, frameshift insertion, stopgain SNV, stoploss SNV, unknown',
        default=['nonsynonymous SNV', 'frameshift deletion', 'frameshift insertion', 'stopgain SNV', 'stoploss SNV', 'unknown'])
Arguments.add_argument('--threshold',
        type=int, choices=range(0,101), metavar="[0-101]",
		help='threshold for minimum percent difference between SNP occurrence in cases and controls (value between 0 and 100)', default=100)
Arguments.add_argument('--occurrence',
        type=int, choices=range(0,101), metavar="[0-101]",
		help='minimum percent of samples in which the SNP appears', default=100)
Arguments.add_argument('--maf',
        type=float,help='MAF < threshold (value between 0 and 1); default is 0.1 for variants', default=0.1)
Arguments.add_argument('--sd',
        type=float,help='segDup < threshold (value between 0 and 1); default is 0.5', default=0.5)
ParsedArgs = Arguments.parse_args()

############
# Functions

def count_samples(filename):
    """A function to define the columns that contain sample genotype data"""
    #   A message that we are going through filename
    print('Reading ' + filename + ' to define columns')
    global chrom, pos, id, format, firstSample, samples, funct, knowngen, exonFcn, hapmap, segdup
    with open(filename, 'r') as f:
        first_line = f.readline()
        p = (first_line.split())
        chrom = p.index('CHROM')
        pos = p.index('POS')
        id = p.index('ID')
        format = p.index('FORMAT')
        funct = p.index('function')
        knowngen = p.index('gene_knowngene')
        exonFcn = p.index('exonicFunction_knowngene')
        hapmap = p.index('hapmap2and3_ASW') # covers all hapmap and 1000genomes and cg/epi values (25 total, excluding dbSNP38 (#21))
        segdup = p.index('segDup')
        firstSample = format + 1 #first sample in filename
        samples = p[firstSample:len(p)] #list of all samples in filename
        print('sample names are' + str(len(samples)) + '.')

def filter_stmp(filename, outfile, mtype, maf, totals, cases, controls, threshold, occurrence, sd):
    """A function to filter vcf or stmp files for candidate homozygous SNPs and/or Genes of Interest (GOI) based on user-defined crtieria """

    print('Reading ' + filename + ' ...' ) #   A message that we are going through filename
    genotype = [None] * len(samples)
    # define counters to keep track of SNPs of interest
    hom=0 # total number of homozygous mutations in all samples listed in file
    het=0 # total number of heterozygous mutations in all samples listed in file
    het_alt=0 # total number of heterozygous (with alt. allele) mutations in all samples listed in file
    case=0 # number of cases that contain this homozygous mutation
    ctrl=0 # number of controls that contain this homozygous mutation
    j=0
    goi_total = {} # goi_total keeps track of # of SNPs in group for gene.

    #begin reading from filename
    with open(filename, 'r') as f:
        casesoutfile = str(outfile + '.cases')
        controlsoutfile = str(outfile + '.ctrls')
        r = open(casesoutfile, 'a')
        c = open(controlsoutfile, 'a')
        print("chrom",'\t',"pos",'\t',"id",'\t',"genes",'\t',"Cases",'\t',"Controls",'\t'.join(str(e) for e in samples),'\t',"homozygous 1/1",'\t',"heterozygous 1/0",'\t',"heterozygous_alt 1/2", file=r)
        print("chrom",'\t',"pos",'\t',"id",'\t',"genes",'\t',"Cases",'\t',"Controls",'\t'.join(str(e) for e in samples),'\t',"homozygous 1/1",'\t',"heterozygous 1/0",'\t',"heterozygous_alt 1/2", file=c)

        next(f)
        for line in f:
            idx = firstSample-1
            fields = line.split('\t')
            curr_gen = fields[knowngen] # keep track of genes associated with each SNP
            curr_genes = re.sub("[(\[].*?[\)\]]","",curr_gen)
            gen_list = curr_genes.split(',')

            # filters for funct, exonFcn, hapmap (MAF compares to all the columns in both hapmap and onekgenomes)
            if (fields[funct] == "exonic" and fields[exonFcn] in mtype):
                for t in range(0,25): #check that the maf is below our threshold for all hapmap and 1000genomes, and cg46    cg69    dbSNP138        esp6500si_ALL   esp6500si_EA    esp6500si_AA    esp_pi columns
                    if t != 21 and (fields[hapmap+t] != '' and float(fields[hapmap+t]) >= maf): # maf is out of range, go to next line in f.
                        j = 1
                if fields[segdup] != '' and fields[segdup] >= sd: #filter out SNPs that fall in repeats and low complexity DNA
                    j = 1
                if j == 0:
                    # once the SNP meets all our criteria, check genotypes for the samples and keep count of occurrences
                    for i in xrange(1, len(samples)+1): # look at the genotype for each sample at the genomic location
                         idx += 1
                         ###  make option for phased genotypes to replace / with | in ./. etc
                         if "./." in fields[idx]:
                             genotype[i-1] = "unk"
                         if "0/0" in fields[idx]:
                             genotype[i-1] = "ref"
                         if ("1/1" in fields[idx]) or ("2/2" in fields[idx]):
                             genotype[i-1] = "hom"
                             hom+=1
                             #if len(samples)>1: #only count cases v controls this for filename with more than one sample
                             for y in cases:
                                 if y == samples[i-1]:
                                     case+=1 # number of times this hom. SNP appears in cases
                             for x in controls:
                                 if x == samples[i-1]:
                                         ctrl+=1 # number of times this hom. SNP appears in controls
                             for gen in gen_list:
                                 if gen in goi_total:
                                     goi_total[gen] += 1 # number of times a particular gene has a hom. SNP appear in any sample in filename
                                 else:
                                     goi_total[gen] = 1 # if we haven't seen a SNP in this gene before, add that gene to the dictionary
                         if ("1/0" in fields[idx]) or ("0/1" in fields[idx]):
                             genotype[i-1] = "het"
                             het+=1
                         if ("1/2" in fields[idx]) or ("2/1" in fields[idx]):
                             genotype[i-1] = "het_alt"
                             het_alt+=1
                    if 'NULL' in controls: # take care of files that have families containing only cases
                        print(fields[chrom],'\t',fields[pos],'\t',fields[id],'\t',';'.join(str(e) for e in gen_list),'\t',case,'\t',ctrl,'\t'.join(str(e) for e in genotype),'\t',hom,'\t',het,'\t',het_alt, file=r)
                        c.close()

                    if 'NULL' in cases: # take care of files that have families containing only controls
                        print(fields[chrom],'\t',fields[pos],'\t',fields[id],'\t',';'.join(str(e) for e in gen_list),'\t',case,'\t',ctrl,'\t'.join(str(e) for e in genotype),'\t',hom,'\t',het,'\t',het_alt, file=c)
                        r.close()

                    if (case + ctrl) > (occurrence/100)*(len(samples)) and case + ctrl != 0: # if SNP appears in sufficient % of our samples
                        if case + ctrl != 0 and abs(100*( (case - ctrl) / (case + ctrl) ) ) > threshold: # if the SNP difference between cases and controls meets a minimum threshold, then print it to the file that we will use for R plotting and GSEA
                            if case>ctrl:
                                print(fields[chrom],'\t',fields[pos],'\t',fields[id],'\t',';'.join(str(e) for e in gen_list),'\t',case,'\t',ctrl,'\t'.join(str(e) for e in genotype),'\t',hom,'\t',het,'\t',het_alt, file=r)
                            if ctrl>case:
                                print(fields[chrom],'\t',fields[pos],'\t',fields[id],'\t',';'.join(str(e) for e in gen_list),'\t',case,'\t',ctrl,'\t'.join(str(e) for e in genotype),'\t',hom,'\t',het,'\t',het_alt, file=c)

            # reset counts and genotypes for next line in filename
            hom=0
            het=0
            het_alt=0
            case=0
            ctrl=0
            genotype = [None] * len(samples)
            j=0
    if not 'NULL' in cases:
        r.close()
    if not 'NULL' in controls:
        c.close()
    if 'NULL' in controls:
        os.remove(controlsoutfile)
    if 'NULL' in cases:
        os.remove(casesoutfile)
    sortedgoi_total = OrderedDict(reversed(sorted(goi_total.items(), key=lambda t: t[1])))
    j = open(totals, 'a')
    print("Gene",'\t',"Total # of Homozygous SNPs in group", file=j)
    for key, value in sortedgoi_total.iteritems() :
        print(key, '\t', value, file=j)
    print(len(sortedgoi_total))
    j.close()

# Carry out commands

count_samples(ParsedArgs.filename)
filter_stmp(ParsedArgs.filename, ParsedArgs.outfile, ParsedArgs.mtype, ParsedArgs.maf, ParsedArgs.totals, ParsedArgs.cases, ParsedArgs.controls, ParsedArgs.threshold, ParsedArgs.occurrence, ParsedArgs.sd)
