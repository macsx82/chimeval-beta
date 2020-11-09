#!/usr/bin/env bash

#script to run the chimeval pipeline

triplet=$1
base_folder=$2
ref_seq=$3
base_out_path=$4
MODE=$5
het_mode=$6


echo ${triplet}
s_name=(${triplets[${triplet}]})
# base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
# ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
s_time=$(date +"%d%m%Y%H%M%S")
# base_out=~/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_${s_time}/STRICT_${strict}/${triplet}
base_out=${base_out_path}/test_${s_time}/STRICT_${strict}/${triplet}
mkdir -p ${base_out}
#the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
for i in 0 1
do
    sample=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[$i]}_r0.fq
    echo "~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA HET_strict" |qsub -N ${s_name[$i]}_chim_s01_${s_time} -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.log -e ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.e -V -l h_vmem=${m} -q all.q,fast
done

#check if there is the need of a custom realignment:

#after the second alignment step, we can go on and work on the chimera
#now the third element will be processed with the script 02
sample_c=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[2]}_r0.fq
ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta

ric_100_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
ric_0_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

#we also need a list of stuff that is only aligned in R100 but not in R0
ric_100_aln_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align.counts
ric_0_aln_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align.counts


echo "~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict} ${ric_100_aln_counts} ${ric_0_aln_counts}" | qsub -N ${s_name[2]}_chim_s02 -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[2]}_chim_s02.log -e ${base_out}/\$JOB_ID_${s_name[2]}_chim_s02.e -V -l h_vmem=${m} -q fast -hold_jid ${s_name[0]}_chim_s01_${s_time},${s_name[1]}_chim_s01_${s_time}
sleep 2




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
case ${het_mode} in
	HET_strict)
		het_lb=0.3
		het_ub=0.6
	;;
	HET_slack)
		het_lb=0.2
		het_ub=0.8
	;;
esac

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

		#the previous script return a lot of files, we need to check for the existence/emptyness of the amplicon list to exclude
		if [[ -s ${sample_out}/${sample_name}_align_not_typed_exclude.list  ]]; then
			#fix this to recursively realign stuff 
			#we will perform only two rounds of alignment
				case ${stringent} in
				ON )
				for aln in ${align_list_to_extract}
				do
				#need a little tuning since this will work only if the ref seq is on one line
				fgrep -w "${aln}" -A 1 ${ref_seq}
				done > ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta
				;;
				OFF )
				#OPT2: we select all the available sequences from the most frequent amplicons for the chimaera:
				for aln in ${align_list_to_extract_ampl}
				do
				#need a little tuning since this will work only if the ref seq is on one line
				fgrep "${aln}_" -A 1 ${ref_seq}
				done > ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta
				;;
				NOCUSTOM )
				#this option in used to align the chimera to all the amplicons alleles
				cp ${ref_seq} ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta
				;;
				esac
		fi

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

