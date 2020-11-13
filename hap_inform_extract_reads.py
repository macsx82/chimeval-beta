#!/usr/bin/env python

"""

:Author: Max Cocca
:Contact: massimiliano.cocca [@] burlo.trieste.it
:Date: *14.02.2018

:Description:

Script to extract informative haplotypes from different kind of dataset

"""

import gzip 
import re 
import sys 
import time
import collections
import itertools
import argparse

# input_file=sys.argv[1]
# # input_file="/home/cocca/analyses/michelangelo/09102017/EUR/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.EUR.vcf.gz"
# n_snp=sys.argv[2]
# # n_snp=5
# chrom=sys.argv[3]
# # chrom=1
# pop=sys.argv[4]
# # pop="EUR"
# out_dir=sys.argv[5]
# # out_dir="/home/cocca/analyses/michelangelo/19102017/EUR"
# block_mode=sys.argv[6]
# # operation_mode=sys.argv[6]
# amp_size=sys.argv[7]
# #amp_size=100 should be the default

parser = argparse.ArgumentParser(description='Script to extract informative haplotypes from different kind of dataset.', version='%(prog)s 02142018.v1.0')

parser.add_argument('-i', action="store",dest="input_file", help="Input file")
parser.add_argument('-n_snp', action="store",dest="n_snp",default=5, help="Minimum number of SNPs in the block",type=int)
parser.add_argument('-chr', action="store",dest="chrom", help="Cromosome")
parser.add_argument('-pop', action="store",dest="pop", help="Population name")
parser.add_argument('-o_path', action="store",dest="out_dir", help="Output data path")
parser.add_argument('-b_mode', action="store",dest="block_mode",choices=('ALL', 'LIST'), help="Define how to generate the blocks. ALL: extract blocks from all data; LIST: extract blocks using a bed file.")
parser.add_argument('-b_list', action="store",dest="block_list", default="", help="List file for block creation")
parser.add_argument('-amp_s', action="store",dest="amp_size", default=100, help="Define the amplicon size", type=int)

all_parameters = parser.parse_args()
input_file = all_parameters.input_file
n_snp = all_parameters.n_snp
chrom = all_parameters.chrom
pop = all_parameters.pop
out_dir = all_parameters.out_dir
block_mode = all_parameters.block_mode
block_list = all_parameters.block_list
amp_size = all_parameters.amp_size


# lines= [if not re.match('#',i): i.split("\t") for i in current_file.readlines()]
#Read gzipped vcf files
current_file=gzip.open('%s' %(input_file), 'r')
print "Reading VCF file...."
lines=[]
for line in current_file:
	if not re.match('#',line):
		site=line.rstrip().split("\t")
		geno=site[9:]
		# need to take care of multiallelic sites
		geno=[w.replace('0', site[3]) for w in geno]
		if re.search(',',site[4]):
		# if len(site[4]) > 1:
			print site[4]
			multi=site[4].split(',')
			print len(multi)
			for n_alt in xrange(1,len(multi)):
				geno=[w.replace(str(n_alt), multi[n_alt-1]) for w in geno]
				# pass
			# geno=[w.replace('2', multi[1]) for w in geno]
			# multi_geno=site
		else:
			geno=[w.replace('1', site[4]) for w in geno]
		# geno=[w.replace('1', site[4]) for w in geno]
		# variant=[site[0],site[1],site[2], site[3],site[4],site[9:]]
		variant=[site[0],site[1],site[2], site[3],site[4],geno]
		lines.append(variant)
print "Vcf file read!"



print "Creating blocks list...."
if block_mode == 'ALL':
	# Create a list of blocks
	#Mode 1:
	all_rs_list=[]
	rs_list=[]
	i=0
	while i <= len(lines):
	# for i in range(0, len(lines)):
		print i
		j =i+int(n_snp)
		if j < len(lines):
			if int(lines[j][1]) - int(lines[i][1]) > int(amp_size): 
				if len(rs_list) != 0:
					# all_rs_list.append(list(set(rs_list)))
					rs_list.sort()
					all_rs_list.append(list(rs_list for rs_list,_ in itertools.groupby(rs_list)))
				rs_list=[]
				i +=1
				continue
			else: 
				for k in range(0, int(n_snp)):
					#print(lines[i+k][2], lines[i+k][4], lines[i+k][9])
					# rs_list.append(lines[i+k][2])
					rs_list.append(lines[i+k])
				while int(lines[j][1]) -int(lines[i][1]) <= int(amp_size):
					#print(lines[j][2], lines[j][4], lines[j][9])
					rs_list.append(lines[j])
					j +=1
				if len(rs_list) != 0:
					# all_rs_list.append(list(set(rs_list)))
					rs_list.sort()
					all_rs_list.append(list(rs_list for rs_list,_ in itertools.groupby(rs_list)))
					rs_list=[]
				i=j
		else:
			break
					# rs_list.append(lines[j][2])
elif block_mode == "LIST":
	#Read tab separated block region
	# block_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_MEA.tab"
	current_list=open('%s' %(block_list), 'r')
	print "Reading list file...."
	r_lines=[]
	for r_line in current_list:
		if re.search("\\b"+str(chrom)+"\\b",r_line):
			block_region=r_line.rstrip().split("\t")
			r_lines.append([block_region[0],block_region[1],block_region[2]])
	print "Block file read!"
	#Mode 2
	#Read User provided file with different block regions:
	# cycle through defined regions and extract alla snps in those boundaries
	r_lines.sort()

	all_rs_list=[]
	for region in r_lines:
		# print region[1]
		rs_list=[]
		start=int(region[1])
		end=int(region[2])
		for b_line in lines:
			if (int(b_line[1]) >= start and int(b_line[1]) <= end):
				rs_list.append(b_line)
				continue
		rs_list.sort()
		all_rs_list.append(list(rs_list for rs_list,_ in itertools.groupby(rs_list)))
		print len(rs_list)

print "Blocks list created!"


print "Write block list file!"
# ${pop}_chr${chr}_${n_snp}_snps_block.txt
#write a file with selected blocks
haps_blocks=open('%s/%s_chr%s_%s_snps_block.txt' %(out_dir,pop,chrom,n_snp), 'w')
for snp_set in all_rs_list:
	# snp_set=all_rs_list[0]
	print >> haps_blocks, 'block\t%s' %(all_rs_list.index(snp_set))
	for snp in snp_set:
		#print the haplotypes
		print >> haps_blocks, '%s\t%s\t%s\t%s\t%s\t%s' %(snp[0],snp[1],snp[2],snp[3],snp[4],"\t".join(snp[5]))

haps_blocks.close()

#for each haplotype we need to get the frequency
# we need a new structure to take in account each sample's hap and check for direct duplicates
# plus duplicated hap which arise from sites with different phasing
#fist we need to get how many samples we have..but to get that number , we have to access the all_rs_list element!
#but I need to do it only once, since the samples' number is always the same
samples=len(all_rs_list[0][0][5])

haps=collections.defaultdict(lambda: collections.defaultdict(list))

# haps={}
#we need to get, for each block, the haplotypes for each sample
#select sample
for sample in xrange(0,samples):
	#enter the hap block
	for block in all_rs_list:
		#enter the snp
		for snp in block:
			#create a dictionary with the hap for that sample for that block
			haps[all_rs_list.index(block)][sample].append(snp[5][sample])
# We need to check if the hap are the same but only with different phase
#but with het sites, they're not the same!!


print "Calculate haps frequencies..."
#now, for each block, we have to count the hap frequency
haps_freq=collections.defaultdict(lambda: collections.defaultdict(list))
for block in haps:
	for sample in haps[block]:
		# print(sample)
		haps_freq[block][tuple(haps[block][sample])].append(sample)

print "Done Calculate haps frequencies!"

#now , lets print out the blocks with their haps and haps frequecies
haps_blocks_freq=open('%s/%s_chr%s_%s_snps_block_freq.txt' %(out_dir,pop,chrom,n_snp), 'w')
print >> haps_blocks_freq,'block\thaplotype\toccurrences\tfreq'
for b in haps_freq:
	for hap in haps_freq[b]:
		hap_fr=len(haps_freq[b][hap])/float(samples)
		print >> haps_blocks_freq, '%s\t%s\t%s\t%s' %(b,hap,len(haps_freq[b][hap]),hap_fr)

haps_blocks_freq.close()
print "Done writing haps frequencies file!"


#now, for each block, I want to get only the most informative haps: aka all those haps within a block 
#that have at least 5 informative snps aka changes, as described here:
# AA vs BB
# BB vs AA
# AB vs BB
# AB vs AA
# BA vs BB
# BA vs AA

#function to extract informativeness of each haplotype in a block
def info_hap(h1,h2):
	info_h1=0
	# info_r1_h1_r1_h2=0
	# info_r1_h1_r2_h2=0
	# info_r2_h1_r1_h2=0
	# info_r2_h1_r2_h2=0
	if len(h1) == len(h2):
		#first extract the first read from h1
		# r1_h1=[x[0] for x in h1]
		# r2_h1=[x[2] for x in h1]
		# r1_h2=[x[0] for x in h2]
		# r2_h2=[x[2] for x in h2]
		for snp_h in range(0,len(h1)):
			# if r1_h1[snp_h] != r1_h2[snp_h]:
			# 	info_r1_h1_r1_h2+=1
			# if r1_h1[snp_h] != r2_h2[snp_h]:
			# 	info_r1_h1_r2_h2+=1
			# if r2_h1[snp_h] != r1_h2[snp_h]:
			# 	info_r2_h1_r1_h2+=1
			# if r2_h1[snp_h] != r2_h2[snp_h]:
			# 	info_r2_h1_r2_h2+=1
			if len(list(set(h1[snp_h].split("|")))) == 1:
				if len(list(set(h2[snp_h].split("|")))) == 1:
					#here we are in the hom -> hom state, informative only if the alleles are different
					if h1[snp_h] != h2[snp_h]:
						info_h1+=1
				# if len(list(set(h2[snp_h].split("|")))) == 2:
					#here we don't add informativeness
			elif len(list(set(h1[snp_h].split("|")))) == 2:
				if len(list(set(h2[snp_h].split("|")))) == 1:
					info_h1+=1
				# if len(list(set(h2[snp_h].split("|")))) == 2:
					#here we don't add informativeness
	#return data for informativenss and the reads correlated
	# return info_r1_h1_r1_h2,info_r1_h1_r2_h2,info_r2_h1_r1_h2,info_r2_h1_r2_h2,r1_h1,r2_h1,r1_h2,r2_h2
	return info_h1


def info_hap_reads(h1,h2):
	info_r1_h1_r1_h2=0
	info_r1_h1_r2_h2=0
	info_r1_h1_r2_h1=0
	info_r2_h1_r1_h2=0
	info_r2_h1_r2_h2=0
	info_r2_h1_r1_h1=0

	if len(h1) == len(h2):
		#first extract the first read from h1
		r1_h1=[x[0] for x in h1]
		r2_h1=[x[2] for x in h1]
		r1_h2=[x[0] for x in h2]
		r2_h2=[x[2] for x in h2]
		for snp_h in range(0,len(h1)):
			if r1_h1[snp_h] != r1_h2[snp_h]:
				info_r1_h1_r1_h2+=1
			if r1_h1[snp_h] != r2_h2[snp_h]:
				info_r1_h1_r2_h2+=1
			if r1_h1[snp_h] != r2_h1[snp_h]:
				info_r1_h1_r2_h1+=1
			if r2_h1[snp_h] != r1_h2[snp_h]:
				info_r2_h1_r1_h2+=1
			if r2_h1[snp_h] != r2_h2[snp_h]:
				info_r2_h1_r2_h2+=1
			if r2_h1[snp_h] != r1_h1[snp_h]:
				info_r2_h1_r1_h1+=1
			# if len(list(set(h1[snp_h].split("|")))) == 1:
			# 	if len(list(set(h2[snp_h].split("|")))) == 1:
			# 		#here we are in the hom -> hom state, informative only if the alleles are different
			# 		if h1[snp_h] != h2[snp_h]:
			# 			info_h1+=1
			# 	# if len(list(set(h2[snp_h].split("|")))) == 2:
			# 		#here we don't add informativeness
			# elif len(list(set(h1[snp_h].split("|")))) == 2:
			# 	if len(list(set(h2[snp_h].split("|")))) == 1:
			# 		info_h1+=1
				# if len(list(set(h2[snp_h].split("|")))) == 2:
					#here we don't add informativeness
	#return data for informativenss and the reads correlated
	return info_r1_h1_r1_h2,info_r1_h1_r2_h2,info_r2_h1_r1_h2,info_r2_h1_r2_h2,r1_h1,r2_h1,r1_h2,r2_h2,info_r1_h1_r2_h1,info_r2_h1_r1_h1
	# return info_h1

print "Calculate haps informativeness..."
haps_sig=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
# haps_sig=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
# haps_selected=collections.defaultdict(lambda: collections.defaultdict(list))

for block in haps_freq.keys():
	for h in range(0,len(haps_freq[block].keys())):
		# h_freqs=[]
		#take the first haplotype and compare it with the other
		for j in range(0,len(haps_freq[block].keys())):
			if haps_freq[block].keys()[h] != haps_freq[block].keys()[j]:
				# print len(haps_freq[block].keys()[h])
				# print len(haps_freq[block].keys()[j])
				hap_info_r1_h1_r1_h2,hap_info_r1_h1_r2_h2,hap_info_r2_h1_r1_h2,hap_info_r2_h1_r2_h2,h_r1,h_r2,h2_r1,h2_r2,hap_info_r1_h1_r2_h1,hap_info_r2_h1_r1_h1 =info_hap_reads(haps_freq[block].keys()[h],haps_freq[block].keys()[j])
				# hap_info_r1_h1_r1_h2,hap_info_r1_h1_r2_h2,hap_info_r2_h1_r1_h2,hap_info_r2_h1_r2_h2,h_r1,h_r2=info_hap(haps_freq[block].keys()[h],haps_freq[block].keys()[j])
				# hap_info=info_hap(haps_freq[block].keys()[h],haps_freq[block].keys()[j])
				#calculate frequencies if hap_info is  >= 3 
				h1_fr=len(haps_freq[block][haps_freq[block].keys()[h]])/float(samples)
				h2_fr=len(haps_freq[block][haps_freq[block].keys()[j]])/float(samples)					
				h_freq= h1_fr * h2_fr
				# h_freqs.append(h_freq)
				# if hap_info_r1_h1_r1_h2 >= 3 | hap_info_r1_h1_r2_h2 >= 3:
				# if hap_info >= 5:
					# haps_sig[block][haps_freq[block].keys()[h]].append(h_freq)
					# haps_sig[block][haps_freq[block].keys()[h]][tuple(haps_freq[block].keys()[j])].append(h_freq)
				# if hap_info_r1_h1_r1_h2 >= 3:
				# 	haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r1)].append(h_freq)
			 # 	else:
			 # 		if hap_info_r1_h1_r2_h2 >= 3:
				# 		haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r1)].append(h_freq)
				# if hap_info_r1_h1_r1_h2 >= 3:
				# 	haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r1)].append(h_freq)
			 # 	else:
			 # 		if hap_info_r1_h1_r2_h2 >= 3:
				# 		haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r1)].append(h_freq)
				# 	else:
				# 		# if a read from the receiver is informative with respect to both the 
				# 		# elif hap_info_r2_h1_r1_h2 >= 3 | hap_info_r2_h1_r2_h2 >= 3:
				# 		if hap_info_r2_h1_r1_h2 >= 3:
				# 			haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r2)].append(h_freq)
				# 		else:
				# 			if hap_info_r2_h1_r2_h2 >= 3:
				# 				haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r2)].append(h_freq)
				# if hap_info_r1_h1_r1_h2 >= 3 and hap_info_r1_h1_r2_h2 >= 3 :
				# 	haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r1)].append(h_freq)
				# else:
				# 	if hap_info_r2_h1_r1_h2 >= 3 and hap_info_r2_h1_r2_h2 >= 3:
				# 		haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r2)].append(h_freq)

				#19/01/2018: modified how to select the informative read! two cases: 1) both reads are informative 2) only one is informative
				# in the second case, we'll compare the two ric reads and we'll keep the read only if it's informative against the other one
				if hap_info_r1_h1_r1_h2 >= 3 and hap_info_r1_h1_r2_h2 >= 3 :
					if hap_info_r2_h1_r1_h2 >= 3 and hap_info_r2_h1_r2_h2 >= 3:
						haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r1)].append(h_freq/2)
						haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r2)].append(h_freq/2)
					else:
						if hap_info_r1_h1_r2_h1 >= 3:
							haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r1)].append(h_freq)
				else:
					if hap_info_r2_h1_r1_h2 >= 3 and hap_info_r2_h1_r2_h2 >= 3:
						if hap_info_r2_h1_r1_h1 >= 3:
							haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r2)].append(h_freq)

				# if hap_info_r1 >= 3:
				# 	haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r1)].append(h1_fr)
				# elif hap_info_r2 >= 3:
				# 	haps_sig[block][haps_freq[block].keys()[h]][tuple(h_r2)].append(h1_fr)
				# haps_sig[block][haps_freq[block].keys()[h]].append(hap_info)
			#now for each haplotype, remove duplicate info
			# if max(list(set(haps_sig[block][haps_freq[block].keys()[h]]))) >= 5:
				#we'll select this haplotype in this block and we'll print out the info
				# haps_selected[block][haps_freq[block].keys()[h]] = list(set(haps_sig[block][haps_freq[block].keys()[h]]))

print "DONE Calculate haps informativeness!"
#CALCULATE CUMULATIVE FREQUENCIES: GET FREQUENCY of each haplotype during comparison, and if haps_info >= 3,
 # keep the freqm than multiply with the second hap, the num all together
#now , lets print out the blocks with their haps and haps frequecies
haps_blocks_sig=open('%s/%s_chr%s_%s_snps_block_sig.txt' %(out_dir,pop,chrom,n_snp), 'w')
# print >> haps_blocks_sig,'block\thaplotype\thaplotype_2\tsig_values\tblock_sig'
# for b in haps_sig:
# 	block_sig=[s for w in haps_sig[b] for k in haps_sig[b][w] for s in haps_sig[b][w][k]]
# 	for hap in haps_sig[b]:
# 		for h_read in haps_sig[b][hap]:
# 			print >> haps_blocks_sig, '%s\t%s\t%s\t%s\t%s' %(b,hap,h_read,sum(haps_sig[b][hap][h_read]),sum(block_sig))
print >> haps_blocks_sig,'block\thaplotype\tread\tsig_values\tblock_sig'
for b in haps_sig:
	block_sig=[z for w in haps_sig[b] for k in haps_sig[b][w] for z in haps_sig[b][w][k] ]
	for hap in haps_sig[b]:
		for hap_read in haps_sig[b][hap]:
			print >> haps_blocks_sig, '%s\t%s\t%s\t%s\t%s' %(b,hap,hap_read,sum(haps_sig[b][hap][hap_read]),sum(block_sig))

haps_blocks_sig.close()

print "DONE writing haps informativeness files!"
# print(all_rs_list)
