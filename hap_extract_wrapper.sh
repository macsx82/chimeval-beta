#!/usr/bin/env bash
hap_TGP_input=$1
n_snp=$2
chr=$3
pop=$4
hap_out_d=$5
hap_block_list=$6


eval "$(conda shell.bash hook)"

conda activate py27

/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 