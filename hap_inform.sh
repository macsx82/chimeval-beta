#!/usr/bin/env bash
#20/01/2017
# Pipeline to extract most informative haplotypes from 1000G phase 3 data

TGP_input=$1
base_out=$2
mode=$3

case $mode in
	LIST )
		#Extract snp list based on input file
		region_list=$4
		subset_sample=$5
		pop=$6
		mkdir -p ${base_out}
		basename_out=`basename ${TGP_input}`
		#only on biallelic snps
		# bcftools view -R ${region_list} -m2 -M2 -v snps -S ${subset_sample} ${TGP_input} -O z -o ${base_out}/${basename_out}.reg_list.${pop}.vcf.gz 
		# bcftools view -R ${region_list} -m2 -M2 -S ${subset_sample} ${TGP_input} -O z -o ${base_out}/${basename_out}.reg_list.${pop}.vcf.gz 
		bcftools view -R ${region_list} -S ${subset_sample} ${TGP_input} -O z -o ${base_out}/${basename_out}.reg_list.${pop}.vcf.gz 
	
	;;
	FREQ )
		maf=$4
		# sd=$5
		subset_sample=$5
		pop=$6
		# TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
		# base_out="/home/cocca/analyses/michelangelo"
		# maf=0.5
		# sd=0.05
		# subset_sample="/home/cocca/analyses/michelangelo/EUR_samples_phase3.txt"
		#add subset Asian and whole TGP3
		mkdir -p ${base_out}
		basename_out=`basename ${TGP_input}`

		# ad_maf1=`echo "$maf $sd" | awk '{print $1-$2}'`
		# ad_maf2=$[maf + sd]
		#extract MAF filtered data only on biallelic snps
		# bcftools view -i"MAF>=${ad_maf1}" -m2 -M2 -v snps -S ${subset_sample} ${TGP_input} -O z -o ${base_out}/${basename_out}.${maf}.${pop}.vcf.gz 
		#23/02/2018: working syntax in v1.7 for maf filtering
		bcftools view -q ${maf}:minor -m2 -M2 -v snps -S ${subset_sample} ${TGP_input} -O z -o ${base_out}/${basename_out}.${maf}.${pop}.vcf.gz 
	;;
	LINKAGE )
		r2=$4
		subset_sample=$5
		pop=$6
		# TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
		# base_out="/home/cocca/analyses/michelangelo"
		# maf=0.5
		# sd=0.05
		# subset_sample="/home/cocca/analyses/michelangelo/EUR_samples_phase3.txt"
		#add subset Asian and whole TGP3
		mkdir -p ${base_out}
		basename_out=`basename ${TGP_input}`

		#extract R2 filtered data only on biallelic snps
		bcftools view -m2 -M2 -v snps -S ${subset_sample} ${TGP_input} | bcftools +prune -l ${r2} -O z -o ${base_out}/${basename_out}.${r2}.${pop}.vcf.gz 
	;;
esac

