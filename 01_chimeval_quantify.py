#!/usr/bin/env python

"""

:Author: Max Cocca
:Contact: massimiliano.cocca [@] burlo.trieste.it
:Date: *09.03.2019

:Description:

Script to extract counts and perform "sequence-typing" based on amplicon frequency on R100 and R0

"""

import re 
import sys 
import time
import collections
import itertools
import argparse


parser = argparse.ArgumentParser(description='Script to extract counts and perform "sequence-typing" based on amplicon frequency on R100 or R0')

parser.add_argument('-i', action="store",dest="input_file", help="Alignment counts")
# parser.add_argument('-a_list', action="store",dest="a_list", help="List of amplicons", nargs="*")
parser.add_argument('-o', action="store",dest="out_file", help="Alignment frequencies output file prefix")


all_parameters = parser.parse_args()
input_file = all_parameters.input_file
# a_list = all_parameters.a_list
out_file = all_parameters.out_file+"_most_freq.counts"
out_file_no_quant=all_parameters.out_file+"_not_typed.counts"
out_file_no_quant_list=all_parameters.out_file+"_not_typed.list"
out_file_exclude_ampl=all_parameters.out_file+"_not_typed_exclude.list"

# input_file ="/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_24072019181711/STRICT_NOCUSTOM/F_H_L/DNA_3_r0.fq/DNA_3_r0.fq_align.counts"
# out_file_prefix = "/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_24072019181711/STRICT_NOCUSTOM/F_H_L/DNA_3_r0.fq/DNA_3_r0.fq_align"
# out_file = out_file_prefix+"_most_freq.counts_test"
# out_file_no_quant=out_file_prefix+"_not_typed.counts_test"
# out_file_no_quant_list=out_file_prefix+"_not_typed.list_test"
# out_file_exclude_ampl=out_file_prefix+"_not_typed_exclude.list_test"

#we need to add 
#1) open the sample file
sample_align=open('%s' %(input_file),'r')
out_file_name=open('%s' %(out_file),'w')
out_file_no_quant_name=open('%s' %(out_file_no_quant),'w')
out_file_no_quant_list_name=open('%s' %(out_file_no_quant_list),'w')
out_file_exclude_ampl_name=open('%s' %(out_file_exclude_ampl),'w')
#2) generate a dictionary of dictionaries: outer key will be the amplicon name, innner keys will be the amplicon allele name, values will be counts
# sample_counts=collections.defaultdict(lambda: collections.defaultdict(list))
sample_counts=collections.defaultdict(dict)

#read all amplicon with alignment for the current sample(in sample_align)
# for aln in a_list:
for line in sample_align:
	#enter the aligned sample counts
	# print(aln)
	ampl = line.strip().split(" ")[1].split("_")[0]
	# if re.match(aln + "_", line.strip().split(" ")[1]):
	# print(line.strip())
	c_align=line.strip().split(" ")
	sample_counts[ampl][c_align[1]]=int(c_align[0])

# 3) Calculate for each amplicon the proportion reads aligned with respect to the overall amplicon reads, to proceed with the quantification
ampl_prop=collections.defaultdict(dict)
ampl_no_quant=collections.defaultdict(dict)
ampl_quant=collections.defaultdict(lambda: collections.defaultdict(dict))
# ampl_quant=collections.defaultdict(collections.defaultdict(dict))

for aln in sample_counts:
	# aln='PCR15'
	#get the sum of all reads aligned for the current apmplicon in the current sample (we are not counting unmapped reads)
	c_aln_sum=sum(sample_counts[aln].values())
	#cicle on the aligned amplicons to sum get the proportion of each allele
	for allele in sample_counts[aln]:
		ampl_prop[aln][allele]=float(sample_counts[aln][allele])/float(c_aln_sum)
		if ampl_prop[aln][allele] >= 0.9:
			#we can say that this is ok, we can skip all other alignment for that amplicon
			ampl_quant["hom"][aln][allele]=ampl_prop[aln][allele]
		elif ampl_prop[aln][allele] >= 0.3 and ampl_prop[aln][allele] <= 0.6:
			ampl_quant["het"][aln][allele]=ampl_prop[aln][allele]
		else:
			ampl_no_quant[allele]=ampl_prop[aln][allele]

# 4) now we can print the normal output
for geno in ("hom","het"):
	for aln in ampl_quant[geno]:
		for ampl in ampl_quant[geno][aln].keys():
			print >> out_file_name, '%s\t%s' %(ampl,ampl_quant[geno][aln][ampl])

out_file_name.close()
#5) we can print also the not typed amplicons alleles
for allele in ampl_no_quant:
	print >> out_file_no_quant_name, '%s\t%s' %(allele,ampl_no_quant[allele])
	print >> out_file_no_quant_list_name, '%s' %(allele)

out_file_no_quant_name.close()
out_file_no_quant_list_name.close()

#6) we need to check if the quantification worked for all amplicons in the correct way, otherwise we need to selcet a list of amplicons to realign the sample
#6a) if we have something in the het lists with only one allele, we need to check all the other alleles and select everithing that is above 10% and realign to it
to_filter = [ampl_het.keys() for ampl_het in ampl_quant['het'].values() if len(ampl_het) < 2]
# if we have something here, it means that we need to extract which amplicons to realign, filtering less represented alleles
to_filter_list = [k.split("_")[0] for j in to_filter for k in j]
# we will add to this list also those amplicon wich are completely untyped
#6b) we check if there are totally unmapped amplicons: those have to be realigned
#get all quantified
all_quant_ampl = list(set([ampl.split("_")[0] for gen in ampl_quant.keys() for ampl_all in ampl_quant[gen].values() for ampl in ampl_all.keys()]))

#compare with the amplicon counts and get all those amplicons which didn't pass the thresholds for the quantification
not_quant_ampl = [x for x in ampl_prop.keys() if x not in all_quant_ampl]
#merge this list to the tofilter list previously calculated
to_extract=list(set(to_filter_list + not_quant_ampl))


# for include_ampl in [k for ampli in to_extract for k,x in ampl_prop[ampli].items() if x >= 0.01]:
for exclude_ampl in [k for ampli in to_extract for k,x in ampl_prop[ampli].items() if x < 0.01]:
	# print >> out_file_include_ampl_name, '%s' %(include_ampl)
	print >> out_file_exclude_ampl_name, '%s' %(exclude_ampl)

out_file_exclude_ampl_name.close()

print("DONE writing alignment count informations!")


