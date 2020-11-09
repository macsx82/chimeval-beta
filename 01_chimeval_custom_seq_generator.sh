#!/usr/bin/env bash
#
#
# Custom fasta file generator from a complete reference and an exlcusio amplicon list
# exclusion_list=$1
# ref_seq=$2
# fasta_out=$3
# stringent=$4
# fasta_out=${sample_out}/REF_SEQ/${sample_name}_custom_f_seq.fasta
# exclusion_list=${sample_out}/${sample_name}_align_not_typed_exclude.list
while getopts ":e:r:o:m:h" opt ${@}; do
  case $opt in
    e)
    # specify exclusion list
		echo "Exclusion list provided"
	    exclusion_list=${OPTARG}
	    ;;
    r)
    #specify reference in fasta format
	    echo "Reference file: ${OPTARG}"
	    ref_seq=${OPTARG}
	    ;;
    o)
      
      echo "Custom reference file : ${OPTARG}"
      fasta_out=${OPTARG}
      ;;
    m)
      echo "Stringency mode selected: ${OPTARG}"
      stringent=${OPTARG}
      ;;
    h)
      echo "#########################"
      echo "Usage:"
      echo "01_chimeval_custom_seq_generator.sh -r <reference_file_path> -o <custom_fasta_output>"
      echo "Execution options: "
      echo "                   -e: Exclusion list file "
      echo "                   -m: Stringency mode selection: ON/OFF/NOCUSTOM "
      echo "                   -h: this help message "
      echo "#########################"
      exit 1
      ;;
    *)
      echo $opt
    ;;
  esac

done


#the previous script return a lot of files, we need to check for the existence/emptyness of the amplicon list to exclude
if [[ -s ${exclusion_list} ]]; then
	#fix this to recursively realign stuff 
	#we will perform only two rounds of alignment
	align_list_to_extract=$(fgrep -v -w -f ${exclusion_list} <(egrep "^>" ${ref_seq} | sed 's/>//g'))
	case ${stringent} in
		ON )
			for aln in ${align_list_to_extract}
			do
			#need a little tuning since this will work only if the ref seq is on one line
				fgrep -w "${aln}" -A 1 ${ref_seq}
			done > ${fasta_out}
		;;
		OFF )
		#OPT2: we select all the available sequences from the most frequent amplicons for the chimaera (the list provided has to be consistent with this option)
			for aln in ${align_list_to_extract}
			do
			#need a little tuning since this will work only if the ref seq is on one line
				fgrep "${aln}_" -A 1 ${ref_seq}
			done > ${fasta_out}
		;;
		NOCUSTOM )
		#this option in used to align the chimera to all the amplicons alleles
			cp ${ref_seq} ${fasta_out}
		;;
		esac
else
	#trigger the NOCUSTOM option
	cp ${ref_seq} ${fasta_out}
fi

touch ${fasta_out}.custom_done