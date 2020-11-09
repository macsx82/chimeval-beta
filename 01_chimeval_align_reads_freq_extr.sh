#!/usr/bin/env bash

#script to perform alignemt on custom designed amplicons and select most frequent alignments from a sample

# base_folder=/home/cocca/analyses/michelangelo/chimerismo/06092018
# sample=RIC_100
# amplicon_base=/home/cocca/analyses/michelangelo/chimerismo/06092018/BED_FILES/SNPS
# triplet=R_D_D
# block=PCR_29
base_folder=$1
sample=$2
ref_seq=$3
base_out=$4
MODE=$5
realign_mode=$6

sample_name=`basename ${sample}`
sample_out=${base_out}/${sample_name}

mkdir -p ${sample_out}

#align using bowtie2
# --very-sensitive-local = -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
#we use Local alignment : When the –local option is specified, Bowtie 2 performs local read alignment.
# In this mode, Bowtie 2 might “trim” or “clip” some read characters from one or both ends of the 
#alignment if doing so maximizes the alignment score.
# Check with Michelangelo or Mea the adapters size to be sure the -3 20 -5 20 parameters are correct

bowtie2 --local --very-sensitive-local --mm -x ${ref_seq} -U ${sample} -3 20 -5 20 --n-ceil L,0,0.5 --no-unal -S ${sample_out}/${sample_name}_test_noK.sam --met-file ${sample_out}/${sample_name}_test_noK.metrics

#count alignments for this sample
samtools view ${sample_out}/${sample_name}_test_noK.sam | cut -f 3 | sort | uniq -c > ${sample_out}/${sample_name}_align.counts

#For each alignment, select the most frequent

#get the list of all alignments to select each group
# sample=RIC_100
# sample_out=~/analyses/michelangelo/chimerismo/NON_HLA/ALL_seq
#Select the HET mode:
# case ${het_mode} in
# 	HET_strict)
# 		het_lb=0.3
# 		het_ub=0.6
# 	;;
# 	HET_slack)
# 		het_lb=0.2
# 		het_ub=0.8
# 	;;
# esac

#Select the operational mode
case ${MODE} in
	HLA)
		#WE need to treat all exons separately
		align_list=$(awk '{print $2}' ${sample_out}/${sample_name}_align.counts| cut -f 1 -d "*"|sort|uniq)

		for aln in ${align_list}
		do
			for exon in EXON_2 EXON_3
			do
			    # sum all alignments for the current amplicons in the current exon
			    total_c_align=$(fgrep "${exon}" ${sample_out}/${sample_name}_align.counts | awk -v c_aln=${aln} '$2~c_aln' | awk '{sum+=$1}END{print sum}')
			    # echo "${aln} ${total_c_align}"
			    #perform the "sequence-typing"
			    fgrep "${exon}" ${sample_out}/${sample_name}_align.counts | awk -v c_aln=${aln} '$2~c_aln' | awk -v tot_c_aln=${total_c_align} '{perc_type=($1/tot_c_aln)}{if(perc_type >= 0.9) print $2,perc_type;else if(perc_type >= 0.3 && perc_type <= 0.6) print $2,perc_type }'
			done >> ${sample_out}/${sample_name}_align_${aln}_most_freq.counts
		done 
		cat ${sample_out}/${sample_name}_align_*_most_freq.counts > ${sample_out}/${sample_name}_align_most_freq_ALL.counts

		for aln in ${align_list}
		do
			for exon in EXON_2 EXON_3
			do
			    # sum all alignments for the current amplicons in the current exon
			    total_c_align=$(fgrep "${exon}" ${sample_out}/${sample_name}_align.counts | awk -v c_aln=${aln} '$2~c_aln' | awk '{sum+=$1}END{print sum}')
			    # echo "${aln} ${total_c_align}"
			    #perform the "sequence-typing"
			    fgrep "${exon}" ${sample_out}/${sample_name}_align.counts | awk -v c_aln=${aln} '$2~c_aln' | awk -v tot_c_aln=${total_c_align} '{perc_type=($1/tot_c_aln)}{if(perc_type < 0.9 && perc_type > 0.6) print $2,perc_type;else if(perc_type < 0.3) print $2,perc_type }'
			done >> ${sample_out}/${sample_name}_align_${aln}_not_typed.counts
		done 
		cat ${sample_out}/${sample_name}_align_*_not_typed.counts > ${sample_out}/${sample_name}_align_most_freq_ALL_not_typed.counts
		cut -f 1 -d " " ${sample_out}/${sample_name}_align_most_freq_ALL_not_typed.counts > ${sample_out}/${sample_name}_align_most_freq_ALL_not_typed.list

	;;
	NON_HLA)
		#get the list of all aligned amplicons
		#call correct the python environment of the script
		source activate py27
		${MY_BASH_SCRIPTS}/01_chimeval_quantify.py -i ${sample_out}/${sample_name}_align.counts -o ${sample_out}/${sample_name}_align

case ${realign_mode} in
	RECURSIVE)
		#the previous script return a lot of files, we need to check for the existence/emptyness of the amplicon list to exclude
		${MY_BASH_SCRIPTS}/01_chimeval_custom_seq_generator.sh -e ${sample_out}/${sample_name}_align_not_typed_exclude.list -r ${ref_seq}.fasta -o ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta -m ON
	;;
	STANDARD)
		touch ${sample_out}/STANDARD_alignment.done
	;;
esac

		# if [[ -s ${sample_out}/${sample_name}_align_not_typed_exclude.list  ]]; then
		# 	#fix this to recursively realign stuff 
		# 	#we will perform only two rounds of alignment
		# 	align_list_to_extract=$(fgrep -v -w -f ${sample_out}/${sample_name}_align_not_typed_exclude.list <(egrep "^>" ${ref_seq}.fasta | sed 's/>//g'))
		# 	case ${stringent} in
		# 		ON )
		# 		for aln in ${align_list_to_extract}
		# 		do
		# 		#need a little tuning since this will work only if the ref seq is on one line
		# 		fgrep -w "${aln}" -A 1 ${ref_seq}.fasta
		# 		done > ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta
		# 		;;
		# 		OFF )
		# 		#OPT2: we select all the available sequences from the most frequent amplicons for the chimaera:
		# 		for aln in ${align_list_to_extract_ampl}
		# 		do
		# 		#need a little tuning since this will work only if the ref seq is on one line
		# 		fgrep "${aln}_" -A 1 ${ref_seq}.fasta
		# 		done > ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta
		# 		;;
		# 		NOCUSTOM )
		# 		#this option in used to align the chimera to all the amplicons alleles
		# 		cp ${ref_seq}.fasta ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta
		# 		;;
		# 		esac
		# fi

		# align_list=$(awk '{print $2}' ${sample_out}/${sample_name}_align.counts| cut -f 1 -d "_"|sort|uniq)

		# for aln in ${align_list}
		# do
		#     # we need to get the sum of all alignments for the current amplicon
		#     total_c_align=$(awk -v c_aln=${aln}_ '$2~c_aln' ${sample_out}/${sample_name}_align.counts | awk '{sum+=$1}END{print sum}')
		#     # echo "${aln} ${total_c_align}"
		#     #perform the "sequence-typing": get the amplicons with the most aligned reads in homozygous or heterozygous state
		#     # awk -v c_aln=${aln}_ '$2~c_aln' ${sample_out}/${sample_name}_align.counts | awk -v tot_c_aln=${total_c_align} '{perc_type=($1/tot_c_aln)}{if(perc_type >= 0.9) print $2,perc_type;else if(perc_type >= 0.3 && perc_type <= 0.6) print $2,perc_type }'
		#     awk -v c_aln=${aln}_ '$2~c_aln' ${sample_out}/${sample_name}_align.counts | awk -v tot_c_aln=${total_c_align} -v het_low=${het_lb} -v het_high=${het_ub} '{perc_type=($1/tot_c_aln)}{if(perc_type >= 0.9) print $2,perc_type;else if(perc_type >= het_low && perc_type <= het_high) print $2,perc_type }'
		# done > ${sample_out}/${sample_name}_align_most_freq.counts

		# for aln in ${align_list}
		# do
		#     # sum all alignments for the current amplicons
		#     total_c_align=$(awk -v c_aln=${aln}_ '$2~c_aln' ${sample_out}/${sample_name}_align.counts | awk '{sum+=$1}END{print sum}')
		#     # echo "${aln} ${total_c_align}"
		#     #perform the "sequence-typing": get the amplicons with are not within the threshold defined before for typing
		#     # awk -v c_aln=${aln}_ '$2~c_aln' ${sample_out}/${sample_name}_align.counts | awk -v tot_c_aln=${total_c_align} '{perc_type=($1/tot_c_aln)}{if(perc_type < 0.9 && perc_type > 0.6) print $2,perc_type;else if(perc_type < 0.3) print $2,perc_type }'
		#     awk -v c_aln=${aln}_ '$2~c_aln' ${sample_out}/${sample_name}_align.counts | awk -v tot_c_aln=${total_c_align} -v het_low=${het_lb} -v het_high=${het_ub} '{perc_type=($1/tot_c_aln)}{if(perc_type < 0.9 && perc_type > het_high) print $2,perc_type;else if(perc_type < het_low) print $2,perc_type }'
		# done > ${sample_out}/${sample_name}_align_not_typed.counts
		# cut -f 1 -d " " ${sample_out}/${sample_name}_align_not_typed.counts > ${sample_out}/${sample_name}_align_not_typed.list
	;;
esac

