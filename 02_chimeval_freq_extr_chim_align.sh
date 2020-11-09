#!/usr/bin/env bash

#script to extract from the fasta file the more frequent alignment to use for the chimera alignment and count the RIC_100 residual

############# funbction definition
list_include() {
    local list="$1"
    local item="$2"
    if [[ $list =~ (^|[[:space:]])"$item"($|[[:space:]]) ]] ; then
    # yes, list include item
        result=0
    else
        result=1
    fi
    return $result
}
#############

base_folder=$1
sample=$2
ref_seq=$3
# ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/custom_amplicon_hand_made.fasta
base_out=$4
ric_0_counts=$5
ric_100_counts=$6
not_typed_donor=$7
# ric_0_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_14022019/RIC_0%_r0.fq/RIC_0%_r0.fq_align.counts
# ric_100_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_14022019/RIC_100_r0.fq/RIC_100_r0.fq_align.counts
#this flag defines whether we'll use a more stringent amplicon ref sequences selection or not
stringent=$8

ric_100_aln_counts=$9
ric_0_aln_counts=${10}

sample_name=`basename ${sample}`
sample_out=${base_out}/${sample_name}_K

mkdir -p ${sample_out}/REF_SEQ

#read both ric_0 and ric_100 more frequent alignment files to extract the fasta data to wich align to
#########################################################################################################################################
#we need also a list of stuff that is not genotyped from the DONOR that we need to remove from the align_list_to_extract list -> WHY??? 
#only from the ric_100!!!
#########################################################################################################################################
# align_list_ric_0=$(awk '{print $1}' ${ric_0_counts} | fgrep -v -w -f ${not_typed_donor} )
# align_list_ric_100=$(awk '{print $1}' ${ric_100_counts} | fgrep -v -w -f ${not_typed_donor} )
#check what happens if we include also the amplicons not genotyped in the donor, wich, supposedly belongs only to the r100
align_list_ric_100=$(awk '{print $1}' ${ric_100_counts} )
align_list_ric_0=$(awk '{print $1}' ${ric_0_counts} )
# align_list_ric_100=$(awk '{print $1}' ${ric_100_counts} )


align_list_to_extract=$(echo "${align_list_ric_0} ${align_list_ric_100}"|tr " " "\n"|sort|uniq)
align_list_to_extract_ampl=$(echo "${align_list_ric_0} ${align_list_ric_100}"|tr " " "\n"|cut -f 1 -d "_"|sort|uniq)

#get align present in both RIC_100 and RIC_0: we need to remove them from the calculation of the residuals
align_list_to_extract_dup=$(echo "${align_list_ric_0} ${align_list_ric_100}"|tr " " "\n"|sort|uniq -c| awk '$1>1 {print $2}')

#I want to extract also a list of amplicons that are aligned in the R100 but not in R0
ric_100_only_aln=$(fgrep -v -w -f <(awk '{print $2}' ${ric_0_aln_counts}) <(awk '{print $2}' ${ric_100_aln_counts}))


#OPT1: we only use the most frequent amplicon sequences as reference for the chimaera:
#extract from the main fasta only the selected alignments
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


#
#index the refseq
bowtie2-build ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta ${sample_out}/REF_SEQ/${sample_name}_custom_f_seq

c_ref_seq=${sample_out}/REF_SEQ/${sample_name}_custom_f_seq

#align using bowtie2
bowtie2 --local --very-sensitive-local --mm -x ${c_ref_seq} -U ${sample} -3 20 -5 20 --n-ceil L,0,0.5 --no-unal -S ${sample_out}/${sample_name}_test_noK.sam --met-file ${sample_out}/${sample_name}_test_noK.metrics

#count alignments for the chimera
samtools view ${sample_out}/${sample_name}_test_noK.sam | cut -f 3 | sort | uniq -c > ${sample_out}/${sample_name}_align.counts

#get the list of all alignments to count residuals of RIC_100

#TODO: refactor with python
align_list=$(awk '{print $2}' ${sample_out}/${sample_name}_align.counts| cut -f 1 -d "_"|sort|uniq)

#call correct the python environment of the script
source activate py27

${MY_BASH_SCRIPTS}/01_chimeval_reads_count.py -i ${sample_out}/${sample_name}_align.counts -a_list ${align_list} -d_list ${align_list_to_extract_dup} -r_list ${align_list_ric_100} -r100_aln_only ${ric_100_only_aln} -o ${sample_out}/${sample_name}_align_ric_100_res.counts

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

