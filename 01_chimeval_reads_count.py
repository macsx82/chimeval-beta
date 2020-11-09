#!/usr/bin/env python

"""

:Author: Max Cocca
:Contact: massimiliano.cocca [@] burlo.trieste.it
:Date: *09.03.2019

:Description:

Script to extract counts and perform "sequence-typing" based on amplicon frequency on the chimaeric sample

"""

import re 
import sys 
import time
import collections
import itertools
import argparse


parser = argparse.ArgumentParser(description='Script to extract counts and perform "sequence-typing" based on amplicon frequency on the chimaeric sample')

parser.add_argument('-i', action="store",dest="input_file", help="Alignment counts")
parser.add_argument('-a_list', action="store",dest="a_list", help="List of amplicons", nargs="*")
parser.add_argument('-d_list', action="store",dest="dupe_list", default="", nargs="*" , help="List file of duplicate sequence alignment")
parser.add_argument('-r_list', action="store",dest="r_list", help="List of R100 aligned sequences", nargs="*")
parser.add_argument('-r100_aln_only', action="store",dest="r100_only_list", help="List of R100 only aligned sequences (not mapped in R0)", nargs="*")

parser.add_argument('-o', action="store",dest="out_file", help="Output file")


all_parameters = parser.parse_args()
input_file = all_parameters.input_file
a_list = all_parameters.a_list

dupe_list = all_parameters.dupe_list
r_list = all_parameters.r_list
r100_only=all_parameters.r100_only_list
out_file = all_parameters.out_file

# input_file ="/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_26062019173657/STRICT_NOCUSTOM/F_G_F/DNA_3_r0.fq_K/DNA_3_r0.fq_align.counts"
# out_file = "/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_26062019173657/STRICT_NOCUSTOM/F_G_F/DNA_3_r0.fq_K/R100_test"
# a_list = ["PCR10","PCR11","PCR12","PCR13","PCR15","PCR16","PCR19","PCR20","PCR21","PCR22","PCR23","PCR25","PCR26","PCR28","PCR29","PCR3","PCR5","PCR6","PCR8","PCR9"]
# dupe_list = ["PCR10_6","PCR11_1","PCR11_3","PCR12_2","PCR12_5","PCR13_4","PCR15_3","PCR16_2","PCR16_3","PCR20_3","PCR21_4","PCR23_1","PCR26_3","PCR28_1","PCR3_5","PCR6_1","PCR8_1","PCR8_4","PCR9_2","PCR9_3"]
# r_list = ["PCR10_6","PCR11_1","PCR11_3","PCR12_2","PCR12_5","PCR13_3","PCR13_4","PCR15_2","PCR15_3","PCR16_2","PCR16_3","PCR19_1","PCR20_3","PCR21_3","PCR21_4","PCR23_1","PCR23_5","PCR26_2","PCR26_3","PCR28_1","PCR28_4","PCR3_5","PCR5_3","PCR5_5","PCR6_1","PCR6_3","PCR8_1","PCR8_4","PCR9_2","PCR9_3"]
###### r_list = ["PCR10_6","PCR11_1","PCR11_3","PCR12_2","PCR12_5","PCR13_3","PCR13_4","PCR15_3","PCR16_2","PCR16_3","PCR20_3","PCR21_4","PCR23_1","PCR26_3","PCR28_1","PCR3_5","PCR6_1","PCR8_1","PCR8_4","PCR9_2","PCR9_3"]
# r100_only=["PCR13_1","PCR13_2","PCR13_3","PCR13_6","PCR22_55","PCR22_60","PCR22_68","PCR9_1"]

# in r_list we need to add the r100_only stuff, removing duplicates
r_all_list = list(set(r_list + r100_only))
################################################################
#we need to understand if this approach is the right one or not!
################################################################

# out_file ="/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_26062019173657/STRICT_NOCUSTOM/F_G_F/DNA_3_r0.fq_K/DNA_3_r0.fq_align_ric_100_res.counts_test"

#we need to add 
#1) open the sample file
sample_align=open('%s' %(input_file),'r')
out_file_name=open('%s' %(out_file),'w')
#2) generate a dictionary of dictionaries: outer key will be the amplicon name, innner keys will be the amplicon allele name, values will be counts
# sample_counts=collections.defaultdict(lambda: collections.defaultdict(list))
sample_counts=collections.defaultdict(dict)

#read all amplicon with alignment for the current sample(in sample_align)
for line in sample_align:
	for aln in a_list:
		#enter the aligned sample counts
		# print(aln)
		if re.match(aln + "_", line.strip().split(" ")[1]):
			# print(line.strip())
			c_align=line.strip().split(" ")
			sample_counts[aln][c_align[1]]=int(c_align[0])

# 3) calculate, for each amplicon the proportion of R100 in the chimaera, checking also if the amplicon allele is in the duplicate list
#Remove all duplicates from the overall sum, and define 2 different sums
for aln in sample_counts:
	#get the sum of all reads aligned for the current apmplicon in the current sample (we are not counting unmapped reads)
	c_aln_sum=sum(sample_counts[aln].values())
	#cicle on the R100 aligned amplicons to sum them up in the chimaera (this is the stuff we have in r_list)
	c_c_r_aln_sum=0
	dupe_f=False
	# for r_aln in r_list:
	for r_aln in r_all_list:
		#check if the amplicon allele is present in the chimaera
		if r_aln in sample_counts[aln].keys():
			#check if this amplicon is in teh dupe list: we need to remove it if this is true
			c_c_r_aln=sample_counts[aln][r_aln]
			# c_c_r_aln_sum=c_c_r_aln_sum+c_c_r_aln
			if r_aln not in dupe_list:
				c_c_r_aln_sum=c_c_r_aln_sum+c_c_r_aln
			else:
				print(r_aln)
				dupe_f=True
			# print(str(c_c_r_aln)+" "+r_aln+" "+str(c_c_r_aln_sum) + " " + str(c_aln_sum))
	#we need to keep track of the amplicons with dupe alignemnts, because we will need to multiply the obtained ratio by 2
	if dupe_f:
		#TODO: check if we need to remove the *2 multiplier
		c_r_prop_c=(float(c_c_r_aln_sum)/float(c_aln_sum))*2*100
		print(str(c_r_prop_c) + " " + aln +" " + str(c_aln_sum) +" " + str(c_c_r_aln_sum) + " DUPE")
	else:
		c_r_prop_c=(float(c_c_r_aln_sum)/float(c_aln_sum))*100
		print(str(c_r_prop_c) + " " + aln +" " + str(c_aln_sum) +" " + str(c_c_r_aln_sum) + " NODUPE")
	print >> out_file_name, '%s\t%s\t%s\t%s' %(aln,c_r_prop_c,c_aln_sum,c_c_r_aln_sum)

out_file_name.close()

print("DONE writing chimaeric informations!")

# #bash form
# for aln in ${align_list}
# do
#     # sum all alignments for the current amplicons
#     total_c_align=$(awk -v c_aln=${aln}_ '$2~c_aln' ${sample_out}/${sample_name}_align.counts | awk '{sum+=$1}END{print sum}')
    
#     #get the count for the most frequent alignments in RIC_100
#     if [[ -z ${align_list_to_extract_dup} ]]; then
#         # echo "vuoto"
#         ric_100_res_count=$(for ric_100_res in ${align_list_ric_100}
#         do
#             #perform the "sequence-typing"
#             awk -v c_aln=${aln}_ '$2~c_aln' ${sample_out}/${sample_name}_align.counts | awk -v ric_100_r=${ric_100_res} '{if($2==ric_100_r) print $1}'
#         done | tr " " "\n"|awk '{sum+=$1}END{print sum}')
#         #get the percentage over the total of the RIC_100 residual in the CHIMERA SAMPLE
#         ric_100_res_perc=$(echo "${total_c_align} ${ric_100_res_count}" | awk '{print ($2*100)/$1}')
#     else
#         ric_100_res_count=$(for ric_100_res in ${align_list_ric_100}
#         do
#             #perform the "sequence-typing"
#             fgrep -v -w "${align_list_to_extract_dup}" ${sample_out}/${sample_name}_align.counts |awk -v c_aln=${aln}_ '$2~c_aln' | awk -v ric_100_r=${ric_100_res} '{if($2==ric_100_r) print $1}'
#         done | tr " " "\n"|awk '{sum+=$1}END{print sum}')        
#         #get the percentage over the total of the RIC_100 residual in the CHIMERA SAMPLE
#         if list_include "${align_list_to_extract_dup}" "${aln}";then
#             ric_100_res_perc=$(echo "${total_c_align} ${ric_100_res_count}" | awk '{print (($2*100)/$1)*2}')
#         else
#             ric_100_res_perc=$(echo "${total_c_align} ${ric_100_res_count}" | awk '{print ($2*100)/$1}')
#         fi
#     fi

#     echo ${aln} ${total_c_align} ${ric_100_res_count} ${ric_100_res_perc}

# done > ${sample_out}/${sample_name}_align_ric_100_res.counts

