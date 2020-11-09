#7/03/2017

Run the extraction for all the different population subset

```bash
for chr in {1..22}
do

for pop in EUR EAS ALL
do

TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
maf="0.5"
sd="0.05"
subset_sample="/home/cocca/analyses/michelangelo/${pop}_samples_phase3.txt"
out_d="/home/cocca/analyses/michelangelo/${pop}"

echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} ${maf} ${sd} ${subset_sample} ${pop}"| qsub -N ${pop}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G

done

done
```

Now extract blocks with genotypes

```bash
for n_snp in 4 5 6

for chr in {1..22}
do

for pop in EUR EAS ALL
do

TGP_input="/home/cocca/analyses/michelangelo/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.5.${pop}.vcf.gz"
out_d=/home/cocca/analyses/michelangelo/21042017/${pop}/
mkdir -p ${out_d}

for n_snp in 5
do

echo "/home/cocca/scripts/bash_scripts/hap_inform_extract.py ${TGP_input} ${n_snp} ${chr} ${pop} ${out_d}"| qsub -N ${pop}_${chr}_${n_snp} -o ${out_d}/\$JOB_ID_${pop}_${chr}_${n_snp}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}_${n_snp}.e -V -l h_vmem=4G

done
done
done

#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract.py ${TGP_input} ${n_snp} > ${out_d}/${pop}_chr${chr}_${n_snp}_snps_block.txt"| qsub -N ${pop}_${chr}_${n_snp} -o ${out_d}/\$JOB_ID_${pop}_${chr}_${n_snp}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}_${n_snp}.e -V -l h_vmem=4G
```
---

#08/06/2017

Rerun the scripts using MAF threshold >= 0 and at least 5 informative SNPs per block

Run the extraction for all the different population subset

```bash
for chr in {1..22}
do

for pop in EUR EAS ALL
do

TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
maf="0"
sd="0"
subset_sample="/home/cocca/analyses/michelangelo/${pop}_samples_phase3.txt"
out_d="/home/cocca/analyses/michelangelo/08062017/${pop}"

mkdir -p ${out_d}
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} ${maf} ${sd} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N ${pop}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G

done

done
```

Now extract blocks with genotypes

```bash
for n_snp in 4 5 6

for chr in {1..22}
do

for pop in EUR EAS ALL
do

#TGP_input="/home/cocca/analyses/michelangelo/08062017/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.5.${pop}.vcf.gz"
TGP_input="/home/cocca/analyses/michelangelo/08062017/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.${pop}.vcf.gz"
out_d=/home/cocca/analyses/michelangelo/08062017/${pop}/
mkdir -p ${out_d}

for n_snp in 5
do

echo "/home/cocca/scripts/bash_scripts/hap_inform_extract.py ${TGP_input} ${n_snp} ${chr} ${pop} ${out_d}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N ${pop}_${chr}_${n_snp} -o ${out_d}/\$JOB_ID_${pop}_${chr}_${n_snp}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}_${n_snp}.e -V -l h_vmem=4G

done
done
done


prefix=`date +"%d%m%Y"`
for chr in {1..22}
do

for pop in ALL EAS EUR
do

TGP_input="/home/cocca/analyses/michelangelo/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.5.${pop}.vcf.gz"
#TGP_input="/home/cocca/analyses/michelangelo/08062017/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.${pop}.vcf.gz"

out_d=/home/cocca/analyses/michelangelo/${prefix}/${pop}/
mkdir -p ${out_d}

for n_snp in 5
do

/home/cocca/scripts/bash_scripts/hap_inform_extract.py ${TGP_input} ${n_snp} ${chr} ${pop} ${out_d}

done
done
done
```

---
#25/09/2017

- Lista blocchi
- Di questi blocchi, calcolare informatività
	- tenendo tutti gli SNPs
	- vedendo almeno due differenze
	- vedendo almeno tre differenze
	- deve ricalcolare l'informatività usando le reads

```bash
for chr in {1..22}
do

for pop in EUR EAS ALL
do

TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/${pop}_samples_phase3.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati.tab"
out_d="/home/cocca/analyses/michelangelo/09102017/${pop}"

mkdir -p ${out_d}
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LIST ${region_list} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N ${pop}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G

done

done
```

---
#15/10/2017

```bash
prefix=`date +"%d%m%Y"`
for chr in {1..22}
do

for pop in ALL EAS EUR
do

#TGP_input="/home/cocca/analyses/michelangelo/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.5.${pop}.vcf.gz"
TGP_input=/home/cocca/analyses/michelangelo/09102017/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz
block_list=/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati.tab
out_d=/home/cocca/analyses/michelangelo/${prefix}/${pop}/
mkdir -p ${out_d}

for n_snp in 5
do

/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py ${TGP_input} ${n_snp} ${chr} ${pop} ${out_d} ${block_list}

done
done
done
```
---
#04/02/2018

Extract new amplicon list

First extract the needed regions:

```bash
for chr in {1..22}
do

for pop in EUR EAS ALL
do

TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab"
out_d="/home/cocca/analyses/michelangelo/22012018/${pop}"

mkdir -p ${out_d}
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LIST ${region_list} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G

#extract info haps
prefix=`date +"%d%m%Y"`
hap_TGP_input=/home/cocca/analyses/michelangelo/22012018/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz
hap_block_list=/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do

echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py ${hap_TGP_input} ${n_snp} ${chr} ${pop} ${hap_out_d} ${hap_block_list}" | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_${pop}_${chr}

done
done
done
```

---

#calculate blocks on the data filtered by linkage

```bash
for chr in {1..22}
do

#for pop in EUR EAS
for pop in ALL
do

TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab"
r2=0.2
out_d="/home/cocca/analyses/michelangelo/23012018_${r2}/${pop}"

#mkdir -p ${out_d}
#echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LINKAGE ${r2} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G

#extract info haps
prefix=`date +"%d%m%Y"`
hap_TGP_input=/home/cocca/analyses/michelangelo/23012018_${r2}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${r2}.${pop}.vcf.gz

hap_block_list=/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${r2}_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"

#/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py ${hap_TGP_input} ${n_snp} ${chr} ${pop} ${hap_out_d} ALL
echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=15G

#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py ${hap_TGP_input} ${n_snp} ${chr} ${pop} ${hap_out_d} ALL" | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=20G -hold_jid extract_${pop}_${chr}

done
done
done
```
---
#14/02/2014
#calculate blocks on the data filtered by linkage on CBM
```bash
for chr in 22
for chr in {1..22}
do

for pop in ALL
do

TGP_input="/home/cocca/analyses/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab"
r2=0.2
out_d="/home/cocca/analyses/michelangelo/14022018_${r2}/${pop}"

mkdir -p ${out_d}
#echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LINKAGE ${r2} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -q burlo

#extract info haps
prefix=`date +"%d%m%Y"`
hap_TGP_input=/home/cocca/analyses/michelangelo/14022018_${r2}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${r2}.${pop}.vcf.gz

hap_block_list=/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${r2}_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"

#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=15G -hold_jid extract_${pop}_${chr} -q burlo
echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=25G -q burlo


done
done
done
```

---
#19/02/2014
#calculate blocks on the data filtered by MAF=30% on CBM/GENEMONSTER

```bash
for maf in 0.1 0.3
do

for pop in EUR EAS ALL
do

for chr in {1..22}
do

TGP_input="/home/cocca/analyses/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab"


prefix_filt=`date +"%d%m%Y%H"`
out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_${maf}/${pop}"

mkdir -p ${out_d}
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} FREQ ${maf} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${maf}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -q burlo

#extract info haps
prefix=`date +"%d%m%Y%H"`
hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz

hap_block_list=/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${maf}_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"

echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=10G -hold_jid extract_${pop}_${maf}_${chr} -q burlo
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=15G -q burlo


done
done
done
done
```
---
#28/02/2018
Test amp size fix

```bash
for maf in 0.1
do

for pop in EAS
do

for chr in {1..22}
do

TGP_input="/home/cocca/analyses/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab"


prefix_filt=`date +"%d%m%Y%H"`
out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_${maf}/${pop}"

mkdir -p ${out_d}
#echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} FREQ ${maf} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${maf}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -q burlo

#extract info haps
prefix=`date +"%d%m%Y%H"`
#hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz
hap_TGP_input=/home/cocca/analyses/michelangelo/2302201803_MAF_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz

hap_block_list=/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${maf}_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"

#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=10G -hold_jid extract_${pop}_${maf}_${chr} -q burlo
echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=10G


done
done
done
done
```
---
#2/03/2018
Extract data with MEA's list

```bash
maf=0.45
pop=EUR
chr=22

for maf in 0.45
do

for pop in EUR EAS ALL
do

for chr in {1..22}
do

TGP_input="/home/cocca/analyses/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
#region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_MEA.tab"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_mea_08032018.txt"


prefix_filt=`date +"%d%m%Y%H"`
#out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_${maf}/${pop}"
#out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MEALIST_${maf}/${pop}"
out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MEALIST_20bl/${pop}"

mkdir -p ${out_d}
#echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} FREQ ${maf} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${maf}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -q burlo
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LIST ${region_list} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_LIST_${pop}_20bl_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G


#extract info haps
prefix=`date +"%d%m%Y%H"`
hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MEALIST_20bl/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MEALIST_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/0203201803_MAF_LIST_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/2302201803_MAF_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz

#hap_block_list=/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab
hap_block_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_mea_08032018.txt"
#hap_block_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_MEA.tab"
#hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${maf}_blocks/${pop}/
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${maf}_MEALIST_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"

#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 100 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_${pop}_${maf}_${chr} -q burlo
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=10G -hold_jid extract_${pop}_${maf}_${chr}
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=10G
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_LIST_${pop}_${maf}_${chr} -q burlo
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -q burlo
echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G


done
done
done
done




for pop in EUR EAS ALL
do

for chr in {1..22}
do

TGP_input="/home/cocca/analyses/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
#region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_MEA.tab"
#egion_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_mea_08032018.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_mea_test_chr6.txt"


prefix_filt=`date +"%d%m%Y%H"`
#out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_${maf}/${pop}"
#out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MEALIST_${maf}/${pop}"
#out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MEALIST_20bl/${pop}"
out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_CHR6_TEST_MEALIST/${pop}"

mkdir -p ${out_d}
#echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} FREQ ${maf} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${maf}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LIST ${region_list} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_LIST_${pop}_20bl_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G


#extract info haps
prefix=`date +"%d%m%Y%H"`
hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_CHR6_TEST_MEALIST/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MEALIST_20bl/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MEALIST_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/0203201803_MAF_LIST_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz
#hap_TGP_input=/home/cocca/analyses/michelangelo/2302201803_MAF_${maf}/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.${maf}.${pop}.vcf.gz

#hap_block_list=/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_04012018.tab
hap_block_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_mea_test_chr6.txt"
#hap_block_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_mea_08032018.txt"
#hap_block_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_tesi_MEA.tab"
#hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${maf}_blocks/${pop}/
#hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${maf}_MEALIST_blocks/${pop}/
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_CHR6_TEST_MEALIST_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"

#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 100 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_${pop}_${maf}_${chr} -q burlo
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=10G -hold_jid extract_${pop}_${maf}_${chr}
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=10G
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_LIST_${pop}_${maf}_${chr} -q burlo
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${maf}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -q burlo
echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_LIST_${pop}_20bl_${chr}


done
done
done
```

---
Reads trimming:

first, go back to FASTA:

```bash
for bam_file in `ls *.bam`
do 

out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/fastq
mkdir -p ${out_path}

echo "samtools fastq -1 ${out_path}/${bam_file}_R1.fq.gz -2 ${out_path}/${bam_file}_R2.fq.gz -O -t ${bam_file} -0 ${out_path}/${bam_file}_R0.fq.gz"| qsub -N bam_to_fasta_${bam_file} -o ${out_path}/\$JOB_ID_bam_to_fasta_${bam_file}.log -e ${out_path}/\$JOB_ID_bam_to_fasta_${bam_file}.e -V -l h_vmem=5G -cwd

done
```


PL=Illumina			#piattaforma
LB=Libx				#libreria
cl=5				#compression level
thr=16				#number of thread
ip1=100				#interval_padding
ip5=500				#interval_padding
cnt=0				#contamination
maa=3				#max alternate alleles
bs=50				#batch size
rt=1				#read thread

```bash
fastq_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/fastq

for fastq in `ls ${fastq_path}/*.fq.gz`
do

ph1=/home/pio/bin/trim_galore

#quality
q=20
#error rate
e=0.1


out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/fastq_trimmed
mkdir -p ${out_path}

echo "${ph1} --gzip -q ${q} -e ${e} -o ${out_path}/ -clip_R1 ${cR1} -three_prime_clip_R1 ${tpcR1} ${fastq}" |  qsub -N trim_fasta_${bam_file} -o ${out_path}/\$JOB_ID_trim_fasta_${bam_file}.log -e ${out_path}/\$JOB_ID_trim_fasta_${bam_file}.e -V -l h_vmem=5G -cwd

done
```
---
Reads trimming:

Using bamUtil trimBam

```bash
ph1=/home/shared/softwares/bin/bam

for bam_file in `ls /home/cocca/analyses/michelangelo/chimerismo/20AMP/*.bam`
do 

#cut 5'
cR1=10

#cut 3'
tpcR1=10

out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/trimmed_bam
final_out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam
read_filt_out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam/100bp_length
mkdir -p ${out_path}
mkdir -p ${final_out_path}
mkdir -p ${read_filt_out_path}

outname=`basename ${bam_file}`
#first trim
#echo "${ph1} trimBam ${bam_file} - -L ${cR1} -R ${tpcR1} --clip | samtools sort -n -O bam -o ${out_path}/${outname}_tempSort " | qsub -N bam_trim_${outname} -o ${out_path}/\$JOB_ID_bam_trim_${outname}.log -e ${out_path}/\$JOB_ID_bam_trim_${outname}.e -V -l h_vmem=5G -cwd

#then fixmate, then sort
#echo "samtools fixmate ${out_path}/${outname}_tempSort - | samtools sort - -O bam -o ${final_out_path}/${outname}" | qsub -N bam_sort_fix_${outname} -o ${final_out_path}/\$JOB_ID_bam_sort_fix_${outname}.log -e ${final_out_path}/\$JOB_ID_bam_sort_fix_${outname}.e -V -l h_vmem=5G -cwd -hold_jid bam_trim_${outname}
#echo "samtools fixmate ${out_path}/${outname}_tempSort - | samtools sort - -O bam -o ${final_out_path}/${outname}" | qsub -N bam_sort_fix_${outname} -o ${final_out_path}/\$JOB_ID_bam_sort_fix_${outname}.log -e ${final_out_path}/\$JOB_ID_bam_sort_fix_${outname}.e -V -l h_vmem=5G -cwd

#echo "samtools view -h -m 100 ${final_out_path}/${outname} -b -o ${read_filt_out_path}/${outname};samtools index ${read_filt_out_path}/${outname}" | qsub -N bam_read_filt_${outname} -o ${read_filt_out_path}/\$JOB_ID_bam_read_filt_${outname}.log -e ${read_filt_out_path}/\$JOB_ID_bam_read_filt_${outname}.e -V -l h_vmem=5G -cwd -hold_jid bam_sort_fix_${outname}
echo "samtools view -h -m 100 ${final_out_path}/${outname} -b -o ${read_filt_out_path}/${outname};samtools index ${read_filt_out_path}/${outname}" | qsub -N bam_read_filt_${outname} -o ${read_filt_out_path}/\$JOB_ID_bam_read_filt_${outname}.log -e ${read_filt_out_path}/\$JOB_ID_bam_read_filt_${outname}.e -V -l h_vmem=5G -cwd

done
```
---
Run chimeval on the filtered data

```bash
out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam/100bp_length/chimeval_blocks
mkdir -p ${out_path}

for chim_set in R_D_R R_D_C1 R_D_C01 R_D_D D_R_D D_R_C1 D_R_C01 D_R_R
do

echo "R CMD BATCH '--args '${chim_set}'' ~/scripts/r_scripts/chimeval_blocks_1.r ${out_path}/chimeval_blocks_1_${chim_set}.Rout" | qsub -N chimeval_blocks_${chim_set} -o ${out_path}/\$JOB_ID_chimeval_blocks_${chim_set}.log -e ${out_path}/\$JOB_ID_chimeval_blocks_${chim_set}.e -V -l h_vmem=5G -cwd

done
```


---
create start end amplicon region bed files

```bash
while read ampl
do

ampl_file=4HotSpot_IAD132502_181_40amp_mod1_${ampl}_sorted.bed

start=`head -1 ${ampl_file}| awk '{print $2}'` 
end=`tail -1 ${ampl_file}| awk '{print $3}'`
chr=`head -1 ${ampl_file} | cut -f 1`

echo -e "${chr}\t${start}\t${end}" > 4HotSpot_IAD132502_181_40amp_mod1_${ampl}_sorted_start_end.bed

done < 4HotSpot_IAD132502_181_40amp.list

```
---
phase bams

```bash
out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam/100bp_length/phased_bam
mkdir -p ${out_path}

for bam_file in *.bam
do

sample_name=${bam_file%.bam}

echo "samtools phase -b ${out_path}/${sample_name} ${bam_file}"| qsub -N phase_bam_${sample_name} -o ${out_path}/\$JOB_ID_phase_bam_${sample_name}.log -e ${out_path}/\$JOB_ID_phase_bam_${sample_name}.e -V -l h_vmem=5G -cwd

done
```

---
Extract regions from bam files

```bash
for ampl in `cat /home/cocca/analyses/michelangelo/chimerismo/4HotSpot_IAD132502_181_40amp.list`
do

echo ${ampl}
base_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam/100bp_length
out_path=${base_path}/extract_regions/${ampl}

mkdir -p ${out_path}

for bam_file in `ls ${base_path}/*.bam`
do

bam_base_name=`basename ${bam_file}`
bam_name=${bam_base_name%.bam}

ampl_file=/home/cocca/analyses/michelangelo/chimerismo/BED_FILES/SNPS/4HotSpot_IAD132502_181_40amp_mod1_${ampl}_sorted.bed

echo "~/scripts/bash_scripts/ampl_region_extract.sh ${bam_file} ${out_path} ${ampl_file}" | qsub -N extract_${ampl}_${bam_name} -o ${out_path}/\$JOB_ID_extract_${ampl}_${bam_name}.log -e ${out_path}/\$JOB_ID_extract_${ampl}_${bam_name}.e -V -l h_vmem=5G -cwd

done
done
```

##########################################
Run chimeval on the amplicon filtered data with and without exclusion of non hom blocks

```bash
for mode_set in HOM
do

out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam/100bp_length/extract_regions/chimeval_blocks_${mode_set}
mkdir -p ${out_path}

for chim_set in R_D_R R_D_C1 R_D_C01 R_D_D D_R_D D_R_C1 D_R_C01 D_R_R
do

echo "R CMD BATCH '--args '${chim_set}' '${mode_set}'' ~/scripts/r_scripts/chimeval_blocks_1_ampl_filt.r ${out_path}/chimeval_blocks_1_ampl_filtr_${chim_set}.Rout" | qsub -N chimeval_blocks_${chim_set}_${mode_set} -o ${out_path}/\$JOB_ID_chimeval_blocks_${chim_set}_${mode_set}.log -e ${out_path}/\$JOB_ID_chimeval_blocks_${chim_set}_${mode_set}.e -V -l h_vmem=5G -cwd

done
done
```
####################################################################
Run chimeval on the unfiltered data but with the exclusion of non hom blocks

```bash
for mode_set in HOM ALL
do

out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam/100bp_length/chimeval_blocks_${mode_set}
mkdir -p ${out_path}

for chim_set in R_D_R R_D_C1 R_D_C01 R_D_D D_R_D D_R_C1 D_R_C01 D_R_R
do

echo "R CMD BATCH '--args '${chim_set}' '${mode_set}'' ~/scripts/r_scripts/chimeval_blocks_1.r ${out_path}/chimeval_blocks_1_${chim_set}.Rout" | qsub -N chimeval_blocks_${chim_set}_${mode_set} -o ${out_path}/\$JOB_ID_chimeval_blocks_${chim_set}_${mode_set}.log -e ${out_path}/\$JOB_ID_chimeval_blocks_${chim_set}_${mode_set}.e -V -l h_vmem=5G -cwd

done
done
```
####################################################################
```bash
for pop in EUR EAS ALL
do
for chr in {1..12}
do

TGP_input="/home/cocca/analyses/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
prefix_filt=`date +"%d%m%Y%H"`
out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_150bp/${pop}"
maf=0.1
mkdir -p ${out_d}
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} FREQ ${maf} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${maf}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=2G -q burlo
prefix=`date +"%d%m%Y%H"`

hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_150bp/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.1.${pop}.vcf.gz
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${maf}_blocks/${pop}

mkdir -p ${hap_out_d}

for n_snp in 5
do echo "${pop} ${chr}"
echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 150 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=15G -hold_jid extract_${pop}_${maf}_${chr} -q burlo 
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 150 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=10G

done
done
done

#####################################################################################
for pop in ALL
do
for chr in {13..22}
do

#TGP_input="/home/cocca/analyses/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
#subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
#prefix_filt=`date +"%d%m%Y%H"`
#out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_150bp/${pop}"
maf=0.1
#mkdir -p ${out_d}
#echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} FREQ ${maf} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_${pop}_${maf}_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=2G -q burlo
prefix=`date +"%d%m%Y%H"`

hap_TGP_input=/home/cocca/analyses/michelangelo/1304201817_150bp/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.1.${pop}.vcf.gz
hap_out_d=/home/cocca/analyses/michelangelo/1304201817_${maf}_blocks/${pop}

mkdir -p ${hap_out_d}

for n_snp in 5
do echo "${pop} ${chr}"
#echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 150 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=15G -hold_jid extract_${pop}_${maf}_${chr} -q burlo 
echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode ALL -amp_s 150 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=15G

done
done
done

```
#####################################################################################
#12/07/2018

phase bams

```bash
out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam/100bp_length/12072018_phased_bam
bam_file=chimera_01.bam

out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/15072018_phased_bam


#m=20G -> good value till 75
q=fast
m=100G

for k in 13 20 50 75 95 100
for k in 90 95 100
do

#k=90

#out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/15072018_phased_bam_NO_A_${k}
out_path=/home/cocca/analyses/michelangelo/chimerismo/20AMP/15072018_phased_bam_${k}

mkdir -p ${out_path}

for bam_file in *.bam
do

sample_name=${bam_file%.bam}


#echo "samtools phase -A -b ${out_path}/${sample_name} ${bam_file}"| qsub -N phase_bam_${sample_name} -o ${out_path}/\$JOB_ID_phase_bam_${sample_name}.log -e ${out_path}/\$JOB_ID_phase_bam_${sample_name}.e -V -l h_vmem=${m} -cwd -q ${q}
#echo "samtools phase -k ${k} -q 30 -b ${out_path}/${sample_name} ${bam_file}"| qsub -N phase_bam_${sample_name}_${k} -o ${out_path}/\$JOB_ID_phase_bam_${sample_name}.log -e ${out_path}/\$JOB_ID_phase_bam_${sample_name}.e -V -l h_vmem=${m} -cwd -q ${q}
#echo "samtools phase -A -k ${k} -q 30 -b ${out_path}/${sample_name} ${bam_file}"| qsub -N phase_bam_${sample_name}_${k} -o ${out_path}/\$JOB_ID_phase_bam_${sample_name}.log -e ${out_path}/\$JOB_ID_phase_bam_${sample_name}.e -V -cwd -q ${q}
echo "samtools phase -A -k ${k} -q 30 -b ${out_path}/${sample_name} ${bam_file}"| qsub -N phase_bam_${sample_name}_${k} -o ${out_path}/\$JOB_ID_phase_bam_${sample_name}.log -e ${out_path}/\$JOB_ID_phase_bam_${sample_name}.e -V -l h_vmem=${m} -cwd -q ${q}

done
done
```

---
file bam check:

```bash
for bam_file in *.bam
do

samtools quickcheck -vvv ${bam_file}

done

index bams:

for bam_f in *.bam;do samtools index ${bam_f};done
```
###############################

Extract data with Michelangelo's list

```bash
maf=0.45
pop=EUR
chr=22

for maf in 0.01
do

for pop in ALL
do

for chr in {1..22}
do

TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/coordinate_ampliconi_PCR_MICHI_31082018.tab"


prefix_filt=`date +"%d%m%Y%H"`
out_d="/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MICHILIST/${pop}"

mkdir -p ${out_d}
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LIST ${region_list} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_LIST_${pop}_20bl_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G


#extract info haps
prefix=`date +"%d%m%Y%H"`
hap_TGP_input=/home/cocca/analyses/michelangelo/${prefix_filt}_MAF_MICHILIST/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz

hap_block_list="/home/cocca/analyses/michelangelo/REG_TAB/coordinate_ampliconi_PCR_MICHI_31082018.tab"
hap_out_d=/home/cocca/analyses/michelangelo/${prefix}_${maf}_MAF_MICHILIST_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"

echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_LIST_${pop}_20bl_${chr}


done
done
done
done
```

---

Extract from Chimera bam file the reads overlapping to the informative snps (or the complete amplicon)

```bash
bam writeRegion --in RIC_0%.bam --out region_test --bamIndex RIC_0%.bai --bed /home/cocca/analyses/michelangelo/chimerismo/06092018/BED_FILES/SNPS/test_region.bed
bam writeRegion --in RIC_100.bam --out region_test_RIC100 --bamIndex RIC_100.bai --bed /home/cocca/analyses/michelangelo/chimerismo/06092018/BED_FILES/SNPS/test_region.bed

current_pos=(31405285 31405378 31405436 31405439)

samtools view region_test_RIC100| cut -f 4,10 | awk -v c_pos="${current_pos[*]}" 'BEGIN{split(c_pos,c_pos_arr," ");}{
start_pos=$1;split($2, seq, "");for(i=1;i<=length(c_pos_arr);i++){    
c_ind=(c_pos_arr[i]-start_pos+1); printf("%s|",seq[c_ind]) };print "";
}' | sort | uniq -c |sort -g -r |more


samtools view region_test_RIC100| cut -f 4,10 | awk -v c_pos="${current_pos[*]}" 'BEGIN{split(c_pos,c_pos_arr," ");}{
start_pos=$1;split($2, seq, ""); i=1;c_ind=(c_pos_arr[i]-start_pos+1); printf("%s %s %s",seq[c_ind],c_pos_arr[i],c_ind);print ""
}' | cut -f 1 -d " "|sort | uniq -c |sort -g -r |more


samtools view region_test_RIC100_samtools| cut -f 4,10 | awk -v c_pos="${current_pos[*]}" 'BEGIN{split(c_pos,c_pos_arr," ");}{
start_pos=$1;split($2, seq, ""); i=1;c_ind=(c_pos_arr[i]-start_pos+1); printf("%s %s %s",seq[c_ind],c_pos_arr[i],c_ind);print ""
}' | cut -f 1 -d " "|sort | uniq -c |sort -g -r |more


# }' | cut -f 1 -d " "|sort | uniq -c |sort -g -r |more

# chr15   31405284        31405285        rs7170825       A G 
# chr15   31405377        31405378        rs12907809      C T 
# chr15   31405435        31405436        rs2338861       C T 
# chr15   31405438        31405439        rs2879276       G A 

# rs7170825  A G
# rs12907809 C T
# rs2338861  C T
# rs2879276  G A



echo "here is a string" | awk '
{ 
  split($0, chars, "")
  for (i=1; i <= length($0); i++) {
    printf("%s\n", chars[i])
  }
}'


print i,c_pos_arr[i],c_ind,seq[c_ind]

```

---
#2/11/2018


Extract data with Michelangelo's list for the thesis:

######USING CBM cluster!##############
```bash
for maf in 0.01
do

for pop in EUR EAS ALL
do

for chr in {1..21}
do

TGP_input="/home/cocca/analyses/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/cocca/analyses/michelangelo/TGP_POPS/${pop}_samples_phase3.txt"
region_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_tesi_michi_02112018.tab"


prefix_filt=`date +"%d%m%Y%H"`
out_d="/home/cocca/analyses/michelangelo/phd_thesis/${prefix_filt}_MAF_MICHILIST/${pop}"

mkdir -p ${out_d}

q=burlo

echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LIST ${region_list} ${subset_sample} ${pop}"| qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N extract_LIST_${pop}_20bl_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -q ${q}


#extract info haps
prefix=`date +"%d%m%Y%H"`
hap_TGP_input=/home/cocca/analyses/michelangelo/phd_thesis/${prefix_filt}_MAF_MICHILIST/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz

hap_block_list="/home/cocca/analyses/michelangelo/REG_TAB/blocchi_sequenziati_tesi_michi_02112018.tab"
hap_out_d=/home/cocca/analyses/michelangelo/phd_thesis/${prefix}_${maf}_MAF_MICHILIST_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"

echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_LIST_${pop}_20bl_${chr} -q ${q}


done
done
done
done
```


---
#3/11/2018

Steps to replicate chimeval counts and genotyping :

- First use chimeval on the selected amplicon
    - use the chimeval_clocks_1.r script
        - use the chimeval_blocks_1_functions.r for function loading (mode_set <- "ALL")
        - set data_set <- c("R_D_D")
        - set j counter for aplicon selection j <- 22



-           

---

Extract only read overlapping each position 

for configuration R_D_D
    - use the appropriate bed file in /home/cocca/analyses/michelangelo/chimerismo/06092018/BED_FILES/SNPS/R_D_D.bed

```bash
amplicon_base=/home/cocca/analyses/michelangelo/chimerismo/06092018/BED_FILES/SNPS

for triplet in R_D_D D_R_R
do

mkdir -p ${amplicon_base}/${triplet}
for block in HLA-AE2 HLA-AE3 HLA-BE2 HLA-BE3 HLA-CE2 HLA-CE3 HLA-DPQ1E2 HLA-DPQ1E3 HLA-DRB1E2 HLA-DRB1E3 HLA-DRP1E2 HLA-DRP1E3 PCR_10 PCR_11 PCR_12 PCR_13 PCR_15 PCR_16 PCR_19 PCR_20 PCR_21 PCR_22 PCR_23 PCR_25 PCR_26 PCR_28 PCR_29 PCR_3 PCR_5 PCR_6 PCR_8 PCR_9
do

amplicon_file=${amplicon_base}/chimeval_HLAandNOHLA_752snps_${block}.bed
bedtools intersect -a ${amplicon_file} -b ${amplicon_base}/${triplet}.bed -wa -wb | cut -f -4,8| sort -g -k2,2 >  ${amplicon_base}/${triplet}/chimeval_HLAandNOHLA_INFORMATIVE_${block}.bed

done
done
```

- Now get the reads for each sample in wich we have all informative snps

```bash

# sample=RIC_100
# triplet=R_D_D
# block=PCR_29

base_folder=/home/cocca/analyses/michelangelo/chimerismo/06092018
amplicon_base=/home/cocca/analyses/michelangelo/chimerismo/06092018/BED_FILES/SNPS
log_folder=${base_folder}/logs

mkdir -p ${log_folder}

for sample in "RIC_01%" "RIC_0%" "RIC_100" "RIC_1%"
do

for triplet in R_D_D D_R_R
do

for block in HLA-AE2 HLA-AE3 HLA-BE2 HLA-BE3 HLA-CE2 HLA-CE3 HLA-DPQ1E2 HLA-DPQ1E3 HLA-DRB1E2 HLA-DRB1E3 HLA-DRP1E2 HLA-DRP1E3 PCR_10 PCR_11 PCR_12 PCR_13 PCR_15 PCR_16 PCR_19 PCR_20 PCR_21 PCR_22 PCR_23 PCR_25 PCR_26 PCR_28 PCR_29 PCR_3 PCR_5 PCR_6 PCR_8 PCR_9
do

echo "/home/cocca/scripts/bash_scripts/chimeval_read_selection.sh ${base_folder} ${sample} ${amplicon_base} ${triplet} ${block}" | qsub -m ea -M massimiliano.cocca@burlo.trieste.it -N read_extr_${block}_${triplet}_${sample} -o ${log_folder}/\$JOB_ID_read_extr_${block}_${triplet}_${sample}.log -e ${log_folder}/\$JOB_ID_read_extr_${block}_${triplet}_${sample}.e -V -l h_vmem=5G -q fast

done
done
done

#bedtools intersect -a ${sample_file} -b <( echo ${line}) -F 1 | samtools view | less -S
#echo ${line}| bedtools intersect -a ${sample_file} -b stdin -F 1 | samtools view | less -S
#bedtools intersect -a ${sample_file} -b <(head -1 ${amplicon_base}/${triplet}/chimeval_HLAandNOHLA_INFORMATIVE_${block}.bed) -F 1| bedtools intersect -a stdin -b <(head -2 ${amplicon_base}/${triplet}/chimeval_HLAandNOHLA_INFORMATIVE_${block}.bed|tail -n 1) -F 1|bedtools intersect -a stdin -b <(tail -n 1 ${amplicon_base}/${triplet}/chimeval_HLAandNOHLA_INFORMATIVE_${block}.bed) -F 1|samtools view |wc -l

samtools view -L ${amplicon_base}/${triplet}/chimeval_HLAandNOHLA_INFORMATIVE_${block}.bed -M ${sample_file} | wc -l
samtools bedcov ${amplicon_base}/${triplet}/chimeval_HLAandNOHLA_INFORMATIVE_${block}.bed ${sample_file}

samtools view ${sample_file} chr21:25405327-25405328 -b | samtools view /dev/stdin chr21:25405349-25405350 -b | samtools view /dev/stdin chr21:25405362-25405363 | wc -l

samtools view -L <(head -1 ${amplicon_base}/${triplet}/chimeval_HLAandNOHLA_INFORMATIVE_${block}.bed) ${sample_file} | wc -l

```

---

#20/01/2019

Test HLAminer

```bash
input_bam=/home/cocca/analyses/michelangelo/chimerismo/06092018/RIC_100.bam


Need to get back to fastq

for sample in RIC_0% RIC_01% RIC_1%
do 

input_bam=/home/cocca/analyses/michelangelo/chimerismo/06092018/${sample}.bam
samtools fastq -1 /home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r1.fq -2 /home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r2.fq -0 /home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r0.fq ${input_bam}

done

out_dir= ~/analyses/michelangelo/HLA

```
---
#25/01/2019

Test HLA-HD

```bash
/home/cocca/softwares/hlahd.1.2.0.1/bin/hlahd.sh -m 100 -c 0.95 -f ~/softwares/hlahd.1.2.0.1/freq_data/ ~/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_100_r0.fq ~/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_100_r0.fq /home/cocca/softwares/hlahd.1.2.0.1/HLA_gene.split.txt /home/cocca/softwares/hlahd.1.2.0.1/dictionary/ RIC_100 ~/analyses/michelangelo/chimerismo/HLA/test_hla_hd/ric_100_test/
```
Questo però sembra volere per forza paired end sequencing...not sure funzioni bene con SE...forse ci saranno dei duplicati o comunque dei results stronzi

```bash
/home/cocca/softwares/HLAminer_v1.4/bin/HPRAwgs_classI-II_SE.sh ~/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_100_r0.fq ~/analyses/michelangelo/chimerismo/HLA/test_fq/RIC_100_r0_30012019
```

Need to get back to fastq

```bash
for sample in RIC_0% RIC_01% RIC_1%
do 

input_bam=/home/cocca/analyses/michelangelo/chimerismo/06092018/${sample}.bam
samtools fastq -1 /home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r1.fq -2 /home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r2.fq -0 /home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r0.fq ${input_bam}

done


out_dir= ~/analyses/michelangelo/HLA


# Test hla-HD in single end mode

echo "/home/cocca/softwares/hlahd.1.2.0.1/bin/hlahd_SE.sh -t 16 -m 100 -c 0.95 -f ~/softwares/hlahd.1.2.0.1/freq_data/ ~/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_100_r0.fq /home/cocca/softwares/hlahd.1.2.0.1/HLA_gene.split.txt /home/cocca/softwares/hlahd.1.2.0.1/dictionary/ RIC_100 /home/cocca/analyses/michelangelo/chimerismo/HLA/test_hla_hd_30012019/" | qsub -N hla_hd_test -V -cwd -l h_vmem=20G -q all.q -o ./\$JOB_ID_hla_hd_test.log -e ./\$JOB_ID_hla_hd_test.err -V -cwd -pe orte 16

input_file=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_100_chr6_r0.fq

echo "/home/cocca/softwares/hlahd.1.2.0.1/bin/hlahd_SE.sh -t 32 -m 100 -c 0.95 -f ~/softwares/hlahd.1.2.0.1/freq_data/ ${input_file} /home/cocca/analyses/michelangelo/chimerismo/HLA/custom_align_data/HLA_gene.split.txt /home/cocca/softwares/hlahd.1.2.0.1/dictionary/ RIC_100 /home/cocca/analyses/michelangelo/chimerismo/HLA/test_hla_hd_1022019/" | qsub -N hla_hd_test -V -cwd -l h_vmem=20G -q all.q -o ./\$JOB_ID_hla_hd_test.log -e ./\$JOB_ID_hla_hd_test.err -V -cwd -pe orte 32
```
---
#11/2/1019

Align all using bowtie2

First, create a bowtie index of the reference

```bash
bowtie2-build custom_fasta_PCR11_19.fasta custom_fasta_PCR11_19
```
Now, align the test file to the data:

Test 1 - no k option set:

```bash
	bowtie2 --local --very-sensitive-local --mm -x ../REF_FASTA/custom_fasta_PCR11_19 -U ../06092018/FASTQ/RIC_100_r0.fq -3 20 -5 20 --n-ceil L,0,0.5 --no-unal -S RIC_100_test_noK.sam
```
1098094 reads; of these:
*  1098094 (100.00%) were unpaired; of these:
*    1032939 (94.07%) aligned 0 times
*    17 (0.00%) aligned exactly 1 time
*    65138 (5.93%) aligned >1 times
*5.93% overall alignment rate


Test 2 - k option set to 1:

```bash
    bowtie2 --local --very-sensitive-local --mm -x ../REF_FASTA/custom_fasta_PCR11_19 -U ../06092018/FASTQ/RIC_100_r0.fq -3 20 -5 20 --n-ceil L,0,0.5 -k 1 --no-unal -S RIC_100_test_k1.sam
    ```

1098094 reads; of these:
*  1098094 (100.00%) were unpaired; of these:
*    1032939 (94.07%) aligned 0 times
*    65155 (5.93%) aligned exactly 1 time
*    0 (0.00%) aligned >1 times
*5.93% overall alignment rate

Run all alignment:

```bash
for sample in RIC_100 RIC_0% RIC_1% RIC_01%
do

bowtie2 --local --very-sensitive-local --mm -x ~/analyses/michelangelo/chimerismo/REF_FASTA/custom_amplicon_hand_made -U ~/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r0.fq -3 20 -5 20 --n-ceil L,0,0.5 --no-unal -S ${sample}_test_noK.sam --met-file ${sample}_test_noK.metrics
bowtie2 --local --very-sensitive-local --mm -x ~/analyses/michelangelo/chimerismo/REF_FASTA/custom_amplicon_hand_made -U ~/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r0.fq -3 20 -5 20 --n-ceil L,0,0.5 -k 1 --no-unal -S ${sample}_test_k1.sam --met-file ${sample}_test_k1.metrics

done

Test only PCR19 most freqent alignments:

for sample in RIC_01% RIC_0% RIC_100
do

bowtie2 --local --very-sensitive-local --mm -x ~/analyses/michelangelo/chimerismo/REF_FASTA/custom_fasta_PCR19 -U ~/analyses/michelangelo/chimerismo/06092018/FASTQ/${sample}_r0.fq -3 20 -5 20 --n-ceil L,0,0.5 --no-unal -S ${sample}_PCR19_test_noK.sam --met-file ${sample}_PCR19_test_noK.metrics

done
```

---

We need to:

- Align using all existing alignment (fasta from defined informative blocks)
- count alignments for RIC_0 (donor) and RIC_100, to get the most frequent
	- we need to use the same criteria used in chimeval:
		- 30-60% of the total -> HET
		- 90% -> HOM
- with the most frequent aligment from RIC_0 and RIC_100, align to CHIMERA
	- calculate the percentage of RIC_100 alignments in CHIMERA



Test script 01: 

```bash
for s_name in 
for s_name in RIC_0% RIC_100 RIC_0%
do

base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name}_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_14022019


~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out}

done
```
Test script 02:

```bash
base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_01%_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_11022019
ric_0_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_11022019/RIC_0%_r0.fq/RIC_0%_r0.fq_align_most_freq.counts
ric_100_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_11022019/RIC_100_r0.fq/RIC_100_r0.fq_align_most_freq.counts


~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample} ${ref_seq} ${base_out} ${ric_0_counts} ${ric_100_counts}
```

---

#Test 14/2/2019

```bash
for s_name in RIC_0% RIC_100
do

base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name}_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_14022019


~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out}

done
```
Test script 02:

```bash
base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_01%_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_14022019
ric_0_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_14022019/RIC_0%_r0.fq/RIC_0%_r0.fq_align_most_freq.counts
ric_100_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_14022019/RIC_100_r0.fq/RIC_100_r0.fq_align_most_freq.counts


~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample} ${ref_seq} ${base_out} ${ric_0_counts} ${ric_100_counts}
```
---
#22/02/2019

Check HLA outcome:

```bash
bowtie2-build HLA_ALL_ex2_ex3_nodup.fasta HLA_ALL_ex2_ex3_nodup
```

Run the first step with RIC_0 and RIC_100 on HLA data

```bash
for s_name in RIC_0% RIC_100
do

base_folder=~/analyses/michelangelo/chimerismo/HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name}_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/HLA_ALL_ex2_ex3_nodup
base_out=~/analyses/michelangelo/chimerismo/HLA/22022019


/usr/bin/time ~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} HLA

done
```

---
#25/2/2019

```bash
for s_name in RIC_0%
do

base_folder=~/analyses/michelangelo/chimerismo/HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name}_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/HLA_ALL_ex2_ex3_nodup
base_out=~/analyses/michelangelo/chimerismo/HLA/25022019


/usr/bin/time ~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} HLA

done
```
```bash
for s_name in RIC_0% RIC_100 
do

base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name}_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_090232019

~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA

done
```

Test script 02:

```bash
base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_01%_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_090232019
ric_0_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_090232019/RIC_0%_r0.fq/RIC_0%_r0.fq_align_most_freq.counts
ric_100_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_090232019/RIC_100_r0.fq/RIC_100_r0.fq_align_most_freq.counts
not_typed_donor=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_090232019/RIC_0%_r0.fq/RIC_0%_r0.fq_align_not_typed.list


~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample} ${ref_seq} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} OFF

~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample} ${ref_seq} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ON
```

---
#09/03/2019

We need to test different combinations

```bash
declare -A triplets=()
# Set 1
triplets[R_D_R]=$(echo "RIC_100 RIC_0% RIC_100")
triplets[R_D_C1]=$(echo "RIC_100 RIC_0% RIC_1%")
triplets[R_D_C01]=$(echo "RIC_100 RIC_0% RIC_01%")
triplets[R_D_D]=$(echo "RIC_100 RIC_0% RIC_0%")

# Set 2
triplets[D_R_D]=$(echo "RIC_0% RIC_100 RIC_0%")
triplets[D_R_C1]=$(echo "RIC_0% RIC_100 RIC_1%")
triplets[D_R_C01]=$(echo "RIC_0% RIC_100 RIC_01%")
triplets[D_R_R]=$(echo "RIC_0% RIC_100 RIC_100")

for strict in ON OFF
do
    for triplet in R_D_R R_D_C1 R_D_C01 R_D_D D_R_D D_R_C1 D_R_C01 D_R_R
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_100232019/STRICT_${strict}/${triplet}

        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name[$i]}_r0.fq
            ~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA
        done

        #now the third element will be processed with the script 02
        sample_c=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name[2]}_r0.fq
        ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
        ric_0_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
        ric_100_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
        not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

        ~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict}
    done
done
```

---

#10/03/2019

Test with align rerun on R_100 after first alignment

```bash
for s_name in RIC_0% RIC_100 
do

base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name}_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_100232019_night

~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA

done

```

Test script 02:

```bash
base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/RIC_01%_r0.fq
ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_090232019
ric_0_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_090232019/RIC_0%_r0.fq/RIC_0%_r0.fq_align_most_freq.counts
ric_100_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_090232019/RIC_100_r0.fq/RIC_100_r0.fq_align_most_freq.counts
not_typed_donor=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/test_090232019/RIC_0%_r0.fq/RIC_0%_r0.fq_align_not_typed.list


~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample} ${ref_seq} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} OFF

~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample} ${ref_seq} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ON
```

Try to downsample at 100x R_0 and R_100:
*R_0 ~ 18000 RD
*R_100 ~ 28000 RD

```bash
samtools view -s 0.006 -b -o RIC_0%.bam ../RIC_0%.bam
samtools view -s 0.004 -b -o RIC_100.bam ../RIC_100.bam
```
now generate fastq data from the downsampled samples:

```bash
for sample in RIC_0% RIC_100
do 

input_bam=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample/${sample}.bam
samtools fastq -1 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample/FASTQ/${sample}_r1.fq -2 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample/FASTQ/${sample}_r2.fq -0 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample/FASTQ/${sample}_r0.fq ${input_bam}

done
```

We need to test different combinations

```bash
declare -A triplets=()
# Set 1
triplets[R_D_R]=$(echo "RIC_100 RIC_0% RIC_100")
triplets[R_D_C1]=$(echo "RIC_100 RIC_0% RIC_1%")
triplets[R_D_C01]=$(echo "RIC_100 RIC_0% RIC_01%")
triplets[R_D_D]=$(echo "RIC_100 RIC_0% RIC_0%")

# Set 2
triplets[D_R_D]=$(echo "RIC_0% RIC_100 RIC_0%")
triplets[D_R_C1]=$(echo "RIC_0% RIC_100 RIC_1%")
triplets[D_R_C01]=$(echo "RIC_0% RIC_100 RIC_01%")
triplets[D_R_R]=$(echo "RIC_0% RIC_100 RIC_100")


for strict in ON
do
    for triplet in R_D_R D_R_D D_R_C1 D_R_C01
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_100232019_night_lowcov/STRICT_${strict}/${triplet}

        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample/FASTQ/${s_name[$i]}_r0.fq
            ~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA
        done

        #now the third element will be processed with the script 02
        sample_c=/home/cocca/analyses/michelangelo/chimerismo/06092018/FASTQ/${s_name[2]}_r0.fq
        ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
        ric_0_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
        ric_100_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
        not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

        ~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict}
    done
done
```
---
```bash
#We need to test different combinations on all LOW coverage

declare -A triplets=()
# Set 1
triplets[R_D_R]=$(echo "RIC_100 RIC_0% RIC_100")
triplets[R_D_C1]=$(echo "RIC_100 RIC_0% RIC_1%")
triplets[R_D_C01]=$(echo "RIC_100 RIC_0% RIC_01%")
triplets[R_D_D]=$(echo "RIC_100 RIC_0% RIC_0%")

# Set 2
triplets[D_R_D]=$(echo "RIC_0% RIC_100 RIC_0%")
triplets[D_R_C1]=$(echo "RIC_0% RIC_100 RIC_1%")
triplets[D_R_C01]=$(echo "RIC_0% RIC_100 RIC_01%")
triplets[D_R_R]=$(echo "RIC_0% RIC_100 RIC_100")


for strict in ON
do
    for triplet in R_D_R D_R_D D_R_C1 D_R_C01
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_100232019_night_lowcov_all/STRICT_${strict}/${triplet}

        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample/FASTQ/${s_name[$i]}_r0.fq
            ~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA
        done

        #now the third element will be processed with the script 02
        sample_c=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample/FASTQ/${s_name[2]}_r0.fq
        ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
        ric_0_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
        ric_100_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
        not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

        ~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict}
    done
done
```

---

Try to downsample at 1000x R_0 and R_100:
*R_0 ~ 18000 RD
*R_100 ~ 28000 RD

```bash
samtools view -s 0.06 -b -o RIC_0%.bam ../RIC_0%.bam
samtools view -s 0.04 -b -o RIC_100.bam ../RIC_100.bam
```
now generate fastq data from the downsampled samples:

```bash
for sample in RIC_0% RIC_100
do 

input_bam=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_1000/${sample}.bam
samtools fastq -1 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_1000/FASTQ/${sample}_r1.fq -2 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_1000/FASTQ/${sample}_r2.fq -0 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_1000/FASTQ/${sample}_r0.fq ${input_bam}

done
```
Now run the scripts:

```bash
declare -A triplets=()
# Set 1
triplets[R_D_R]=$(echo "RIC_100 RIC_0% RIC_100")
triplets[R_D_C1]=$(echo "RIC_100 RIC_0% RIC_1%")
triplets[R_D_C01]=$(echo "RIC_100 RIC_0% RIC_01%")
triplets[R_D_D]=$(echo "RIC_100 RIC_0% RIC_0%")

# Set 2
triplets[D_R_D]=$(echo "RIC_0% RIC_100 RIC_0%")
triplets[D_R_C1]=$(echo "RIC_0% RIC_100 RIC_1%")
triplets[D_R_C01]=$(echo "RIC_0% RIC_100 RIC_01%")
triplets[D_R_R]=$(echo "RIC_0% RIC_100 RIC_100")


for strict in ON
do
    for triplet in R_D_R D_R_D D_R_C1 D_R_C01
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_100232019_night_lowcov_1000/STRICT_${strict}/${triplet}

        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_1000/FASTQ/${s_name[$i]}_r0.fq
            ~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA
        done

        #now the third element will be processed with the script 02
        sample_c=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_1000/FASTQ/${s_name[2]}_r0.fq
        ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
        ric_0_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
        ric_100_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
        not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

        ~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict}
    done
done
```

---
Try to downsample at 5000x R_0 and R_100:
*R_0 ~ 18000 RD
*R_100 ~ 28000 RD

```bash
samtools view -s 0.28 -b -o RIC_0%.bam ../RIC_0%.bam
samtools view -s 0.18 -b -o RIC_100.bam ../RIC_100.bam
```
now generate fastq data from the downsampled samples:

```bash
for sample in RIC_0% RIC_100
do 

input_bam=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_5000/${sample}.bam
samtools fastq -1 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_5000/FASTQ/${sample}_r1.fq -2 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_5000/FASTQ/${sample}_r2.fq -0 /home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_5000/FASTQ/${sample}_r0.fq ${input_bam}

done
```
Now run the scripts:

```bash
declare -A triplets=()
# Set 1
triplets[R_D_R]=$(echo "RIC_100 RIC_0% RIC_100")
triplets[R_D_C1]=$(echo "RIC_100 RIC_0% RIC_1%")
triplets[R_D_C01]=$(echo "RIC_100 RIC_0% RIC_01%")
triplets[R_D_D]=$(echo "RIC_100 RIC_0% RIC_0%")

# Set 2
triplets[D_R_D]=$(echo "RIC_0% RIC_100 RIC_0%")
triplets[D_R_C1]=$(echo "RIC_0% RIC_100 RIC_1%")
triplets[D_R_C01]=$(echo "RIC_0% RIC_100 RIC_01%")
triplets[D_R_R]=$(echo "RIC_0% RIC_100 RIC_100")


for strict in ON
do
    for triplet in R_D_R D_R_D D_R_C1 D_R_C01
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_100232019_night_lowcov_5000/STRICT_${strict}/${triplet}

        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_5000/FASTQ/${s_name[$i]}_r0.fq
            ~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA
        done

        #now the third element will be processed with the script 02
        sample_c=/home/cocca/analyses/michelangelo/chimerismo/06092018/downsample_5000/FASTQ/${s_name[2]}_r0.fq
        ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
        ric_0_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
        ric_100_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
        not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

        ~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict}
    done
done
```
---
#15/05/2019

Downloaded new data from Torrent site

/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019
rename files:

```bash
mv IonXpress_002_R_2019_05_09_09_30_00_user_SN2-758-Chimerismo_HLA_chip4_Auto_user_SN2-758-Chimerismo_HLA_chip4_1253.bam.bai DNA_13.bam.bai
mv IonXpress_003_R_2019_05_09_09_30_00_user_SN2-758-Chimerismo_HLA_chip4_Auto_user_SN2-758-Chimerismo_HLA_chip4_1253.bam.bai DNA_11_PADRE.bam.bai
mv IonXpress_004_R_2019_05_09_09_30_00_user_SN2-758-Chimerismo_HLA_chip4_Auto_user_SN2-758-Chimerismo_HLA_chip4_1253.bam.bai DNA_12_MADRE.bam.bai
mv IonXpress_005_R_2019_05_09_09_30_00_user_SN2-758-Chimerismo_HLA_chip4_Auto_user_SN2-758-Chimerismo_HLA_chip4_1253.bam.bai Sc_I_0.1.bam.bai
mv IonXpress_001_R_2019_05_09_09_30_00_user_SN2-758-Chimerismo_HLA_chip4_Auto_user_SN2-758-Chimerismo_HLA_chip4_1253.bam.bai Sc_L_0.1.bam.bai
```

To get:

* DNA_13.bam
* DNA_11_PADRE.bam
* DNA_12_MADRE.bam
* Sc_I_0.1.bam
* Sc_L_0.1.bam

Need to get back to fastq:

```bash
mkdir -p /home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ

for sample in DNA_13 DNA_11_PADRE DNA_12_MADRE Sc_I_0.1 Sc_L_0.1
do 

input_bam=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/${sample}.bam
samtools fastq -1 /home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${sample}_r1.fq -2 /home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${sample}_r2.fq -0 /home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${sample}_r0.fq ${input_bam}

done
```

Now run the scripts:

```bash
A=DNA_13
B=DNA_11_PADRE
C=DNA_12_MADRE

D=Sc_I_0.1
E=Sc_L_0.1

source activate py27

declare -A triplets=()
# Set 1 - Bianco
triplets[A_B_B]=$(echo "DNA_13 DNA_11_PADRE DNA_11_PADRE")
triplets[A_C_C]=$(echo "DNA_13 DNA_12_MADRE DNA_12_MADRE")
triplets[B_A_A]=$(echo "DNA_11_PADRE DNA_13 DNA_13")
triplets[B_C_C]=$(echo "DNA_11_PADRE DNA_12_MADRE DNA_12_MADRE")
triplets[C_B_B]=$(echo "DNA_12_MADRE DNA_11_PADRE DNA_11_PADRE")
triplets[C_A_A]=$(echo "DNA_12_MADRE DNA_13 DNA_13")

#Set 2 - Limit of detection
triplets[B_B_A]=$(echo "DNA_11_PADRE DNA_11_PADRE DNA_13")
triplets[A_B_A]=$(echo "DNA_13 DNA_11_PADRE DNA_13")
triplets[A_B_D]=$(echo "DNA_13 DNA_11_PADRE Sc_I_0.1")
triplets[A_C_E]=$(echo "DNA_13 DNA_12_MADRE Sc_L_0.1")
m=6G

for strict in ON
do
    #for triplet in B_C_C C_B_B C_A_A A_B_D A_C_E
    #for triplet in A_B_B A_C_C
    #for triplet in A_B_A
    #for triplet in B_B_A
    for triplet in B_A_A B_C_C C_B_B C_A_A A_B_D A_C_E A_B_B A_C_C
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_140620190102/STRICT_${strict}/${triplet}
        mkdir -p ${base_out}
        s_time=$(date +"%d%m%Y%H%M%S")
        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[$i]}_r0.fq
            echo "~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA" |qsub -N ${s_name[$i]}_chim_s01_${s_time} -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.log -e ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.e -V -l h_vmem=${m} -q fast
        done

        #now the third element will be processed with the script 02
        sample_c=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[2]}_r0.fq
        ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
        
        ric_100_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
        ric_0_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
        not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

        echo "~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict}" | qsub -N ${s_name[2]}_chim_s02 -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[2]}_chim_s02.log -e ${base_out}/\$JOB_ID_${s_name[2]}_chim_s02.e -V -l h_vmem=${m} -q fast -hold_jid ${s_name[0]}_chim_s01_${s_time},${s_name[1]}_chim_s01_${s_time}
        sleep 2
    done
done
```

---
Check HLA HD options to use custom HLA reference

```bash
declare -A triplets=()
# Set 1 - Bianco
triplets[A_B_B]=$(echo "DNA_13 DNA_11_PADRE DNA_11_PADRE")
triplets[A_C_C]=$(echo "DNA_13 DNA_12_MADRE DNA_12_MADRE")
triplets[B_A_A]=$(echo "DNA_11_PADRE DNA_13 DNA_13")
triplets[B_C_C]=$(echo "DNA_11_PADRE DNA_12_MADRE DNA_12_MADRE")
triplets[C_B_B]=$(echo "DNA_12_MADRE DNA_11_PADRE DNA_11_PADRE")
triplets[C_A_A]=$(echo "DNA_12_MADRE DNA_13 DNA_13")

#Set 2 - Limit of detection
triplets[B_B_A]=$(echo "DNA_11_PADRE DNA_11_PADRE DNA_13")
triplets[A_B_A]=$(echo "DNA_13 DNA_11_PADRE DNA_13")
triplets[A_B_D]=$(echo "DNA_13 DNA_11_PADRE Sc_I_0.1")
triplets[A_C_E]=$(echo "DNA_13 DNA_12_MADRE Sc_L_0.1")
m=6G

for strict in NOCUSTOM
do
    #for triplet in B_C_C C_B_B C_A_A A_B_D A_C_E
    #for triplet in A_B_B A_C_C
    #for triplet in A_B_A
    #for triplet in B_B_A
    for triplet in B_A_A B_C_C C_B_B C_A_A A_B_D A_C_E A_B_B A_C_C
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        s_time=$(date +"%d%m%Y%H%M%S")
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_${s_time}/STRICT_${strict}/${triplet}
        mkdir -p ${base_out}
        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[$i]}_r0.fq
            echo "~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA" |qsub -N ${s_name[$i]}_chim_s01_${s_time} -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.log -e ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.e -V -l h_vmem=${m} -q fast
        done

        #now the third element will be processed with the script 02
        sample_c=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[2]}_r0.fq
        ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
        
        ric_100_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
        ric_0_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
        not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

        echo "~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict}" | qsub -N ${s_name[2]}_chim_s02 -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[2]}_chim_s02.log -e ${base_out}/\$JOB_ID_${s_name[2]}_chim_s02.e -V -l h_vmem=${m} -q fast -hold_jid ${s_name[0]}_chim_s01_${s_time},${s_name[1]}_chim_s01_${s_time}
        sleep 2
    done
done
```

---
#20/06/2019

Downloaded new data from Torrent site

/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019

rename files:

```bash
mv IonXpress_001_rawlib.basecaller.bam DNA_1_PADRE.bam
mv IonXpress_002_rawlib.basecaller.bam DNA_2_MADRE.bam
mv IonXpress_003_rawlib.basecaller.bam DNA_3.bam
mv IonXpress_004_rawlib.basecaller.bam Sc_A_0.1.bam
mv IonXpress_005_rawlib.basecaller.bam Sc_E_0.1.bam
```

To get:
* DNA_1_PADRE.bam
* DNA_2_MADRE.bam
* DNA_3.bam
* Sc_A_0.1.bam
* Sc_E_0.1.bam

Need to get back to fastq:

```bash
mkdir -p /home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ

for sample in DNA_1_PADRE DNA_2_MADRE DNA_3 Sc_A_0.1 Sc_E_0.1
do 

input_bam=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/${sample}.bam
samtools fastq -1 /home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${sample}_r1.fq -2 /home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${sample}_r2.fq -0 /home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${sample}_r0.fq ${input_bam}

done
```

Now run the scripts:

```bash
A=DNA_13
B=DNA_11_PADRE
C=DNA_12_MADRE

D=Sc_I_0.1
E=Sc_L_0.1

F=DNA_3
G=DNA_1_PADRE
H=DNA_2_MADRE
I=Sc_A_0.1
L=Sc_E_0.1

declare -A triplets=()
# Set 1 - Bianco
triplets[F_G_G]=$(echo "DNA_3 DNA_1_PADRE DNA_1_PADRE")
triplets[F_H_H]=$(echo "DNA_3 DNA_2_MADRE DNA_2_MADRE")
triplets[G_F_F]=$(echo "DNA_1_PADRE DNA_3 DNA_3")
triplets[G_H_H]=$(echo "DNA_1_PADRE DNA_2_MADRE DNA_2_MADRE")
triplets[H_G_G]=$(echo "DNA_2_MADRE DNA_1_PADRE DNA_1_PADRE")
triplets[H_F_F]=$(echo "DNA_2_MADRE DNA_3 DNA_3")

#Set 2 - Limit of detection
triplets[G_G_F]=$(echo "DNA_1_PADRE DNA_1_PADRE DNA_3")
triplets[H_H_F]=$(echo "DNA_2_MADRE DNA_2_MADRE DNA_3")
triplets[F_G_F]=$(echo "DNA_3 DNA_1_PADRE DNA_3")
triplets[F_H_F]=$(echo "DNA_3 DNA_2_MADRE DNA_3")
triplets[F_G_I]=$(echo "DNA_3 DNA_1_PADRE Sc_A_0.1")
triplets[F_H_L]=$(echo "DNA_3 DNA_2_MADRE Sc_E_0.1")
m=6G

for strict in NOCUSTOM
do
    for triplet in F_G_G F_H_H G_F_F G_H_H H_G_G H_F_F G_G_F H_H_F F_G_F F_H_F F_G_I F_H_L
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        s_time=$(date +"%d%m%Y%H%M%S")
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/test_${s_time}/STRICT_${strict}/${triplet}
        mkdir -p ${base_out}
        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[$i]}_r0.fq
            echo "~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA" |qsub -N ${s_name[$i]}_chim_s01_${s_time} -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.log -e ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.e -V -l h_vmem=${m} -q fast
        done

        #now the third element will be processed with the script 02
        sample_c=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[2]}_r0.fq
        ref_seq_2=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019.fasta
        
        ric_100_counts=${base_out}/${s_name[0]}_r0.fq/${s_name[0]}_r0.fq_align_most_freq.counts
        ric_0_counts=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_most_freq.counts
        not_typed_donor=${base_out}/${s_name[1]}_r0.fq/${s_name[1]}_r0.fq_align_not_typed.list

        echo "~/scripts/bash_scripts/02_chimeval_freq_extr_chim_align.sh ${base_folder} ${sample_c} ${ref_seq_2} ${base_out} ${ric_0_counts} ${ric_100_counts} ${not_typed_donor} ${strict}" | qsub -N ${s_name[2]}_chim_s02 -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[2]}_chim_s02.log -e ${base_out}/\$JOB_ID_${s_name[2]}_chim_s02.e -V -l h_vmem=${m} -q fast -hold_jid ${s_name[0]}_chim_s01_${s_time},${s_name[1]}_chim_s01_${s_time}
        sleep 2
    done
done


declare -A triplets=()
# Set 1 - Bianco
triplets[F_G_G]=$(echo "DNA_3 DNA_1_PADRE DNA_1_PADRE")
triplets[F_H_H]=$(echo "DNA_3 DNA_2_MADRE DNA_2_MADRE")
triplets[G_F_F]=$(echo "DNA_1_PADRE DNA_3 DNA_3")
triplets[G_H_H]=$(echo "DNA_1_PADRE DNA_2_MADRE DNA_2_MADRE")
triplets[H_G_G]=$(echo "DNA_2_MADRE DNA_1_PADRE DNA_1_PADRE")
triplets[H_F_F]=$(echo "DNA_2_MADRE DNA_3 DNA_3")

#Set 2 - Limit of detection
triplets[G_G_F]=$(echo "DNA_1_PADRE DNA_1_PADRE DNA_3")
triplets[H_H_F]=$(echo "DNA_2_MADRE DNA_2_MADRE DNA_3")
triplets[F_G_F]=$(echo "DNA_3 DNA_1_PADRE DNA_3")
triplets[F_H_F]=$(echo "DNA_3 DNA_2_MADRE DNA_3")
triplets[F_G_I]=$(echo "DNA_3 DNA_1_PADRE Sc_A_0.1")
triplets[F_H_L]=$(echo "DNA_3 DNA_2_MADRE Sc_E_0.1")
m=6G

    #for triplet in F_G_F F_H_F F_G_G F_H_H G_F_F G_H_H H_G_G H_F_F G_G_F H_H_F F_G_I F_H_L
for strict in NOCUSTOM
do
    for triplet in F_G_I F_H_L
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        s_time=$(date +"%d%m%Y%H%M%S")
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_${s_time}/STRICT_${strict}/${triplet}
        mkdir -p ${base_out}
        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[$i]}_r0.fq
            echo "~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA" |qsub -N ${s_name[$i]}_chim_s01_${s_time} -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.log -e ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.e -V -l h_vmem=${m} -q fast
        done

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
    done
done
```
---
#23/07/2019

Run the first steps with a more relaxed het threshold:

```bash
declare -A triplets=()
# Set 1 - Bianco
triplets[A_B_B]=$(echo "DNA_13 DNA_11_PADRE DNA_11_PADRE")
triplets[A_C_C]=$(echo "DNA_13 DNA_12_MADRE DNA_12_MADRE")
triplets[B_A_A]=$(echo "DNA_11_PADRE DNA_13 DNA_13")
triplets[B_C_C]=$(echo "DNA_11_PADRE DNA_12_MADRE DNA_12_MADRE")
triplets[C_B_B]=$(echo "DNA_12_MADRE DNA_11_PADRE DNA_11_PADRE")
triplets[C_A_A]=$(echo "DNA_12_MADRE DNA_13 DNA_13")

# Set 1 - Bianco
triplets[F_G_G]=$(echo "DNA_3 DNA_1_PADRE DNA_1_PADRE")
triplets[F_H_H]=$(echo "DNA_3 DNA_2_MADRE DNA_2_MADRE")
triplets[G_F_F]=$(echo "DNA_1_PADRE DNA_3 DNA_3")
triplets[G_H_H]=$(echo "DNA_1_PADRE DNA_2_MADRE DNA_2_MADRE")
triplets[H_G_G]=$(echo "DNA_2_MADRE DNA_1_PADRE DNA_1_PADRE")
triplets[H_F_F]=$(echo "DNA_2_MADRE DNA_3 DNA_3")

#Set 2 - Limit of detection
triplets[B_B_A]=$(echo "DNA_11_PADRE DNA_11_PADRE DNA_13")
triplets[C_C_A]=$(echo "DNA_12_MADRE DNA_12_MADRE DNA_13")
triplets[A_B_A]=$(echo "DNA_13 DNA_11_PADRE DNA_13")
triplets[A_C_A]=$(echo "DNA_13 DNA_12_MADRE DNA_13")
triplets[A_B_D]=$(echo "DNA_13 DNA_11_PADRE Sc_I_0.1")
triplets[A_C_E]=$(echo "DNA_13 DNA_12_MADRE Sc_L_0.1")
#Set 2 - Limit of detection
triplets[G_G_F]=$(echo "DNA_1_PADRE DNA_1_PADRE DNA_3")
triplets[H_H_F]=$(echo "DNA_2_MADRE DNA_2_MADRE DNA_3")
triplets[F_G_F]=$(echo "DNA_3 DNA_1_PADRE DNA_3")
triplets[F_H_F]=$(echo "DNA_3 DNA_2_MADRE DNA_3")
triplets[F_G_I]=$(echo "DNA_3 DNA_1_PADRE Sc_A_0.1")
triplets[F_H_L]=$(echo "DNA_3 DNA_2_MADRE Sc_E_0.1")
m=6G

for strict in NOCUSTOM
do
    for triplet in A_B_B A_C_C B_A_A B_C_C C_B_B C_A_A B_B_A C_C_A A_B_A A_C_A A_B_D A_C_E F_G_G F_H_H G_F_F G_H_H H_G_G H_F_F G_G_F H_H_F F_G_F F_H_F F_G_I F_H_L
    do
        echo ${triplet}
        s_name=(${triplets[${triplet}]})
        base_folder=~/analyses/michelangelo/chimerismo/NON_HLA
        ref_seq=/home/cocca/analyses/michelangelo/chimerismo/REF_FASTA/fasta_non-HLA_13022019
        s_time=$(date +"%d%m%Y%H%M%S")
        base_out=~/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_${s_time}/STRICT_${strict}/${triplet}
        mkdir -p ${base_out}
        #the first 2 elements of the triplet are processed with script 01, the third will be processed after the first two finished
        for i in 0 1
        do
            sample=/home/cocca/analyses/michelangelo/chimerismo/new_bams_5052019/FASTQ/${s_name[$i]}_r0.fq
            echo "~/scripts/bash_scripts/01_chimeval_align_reads_freq_extr.sh ${base_folder} ${sample} ${ref_seq} ${base_out} NON_HLA HET_strict" |qsub -N ${s_name[$i]}_chim_s01_${s_time} -m ea -M massimiliano.cocca@burlo.trieste.it -o ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.log -e ${base_out}/\$JOB_ID_${s_name[$i]}_chim_s01.e -V -l h_vmem=${m} -q all.q,fast
        done

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
    done
done
```
---
check tests:

```bash
ric_100_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_26062019173657/STRICT_NOCUSTOM/F_G_F/DNA_3_r0.fq/DNA_3_r0.fq_align_most_freq.counts
ric_0_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_26062019173657/STRICT_NOCUSTOM/F_G_F/DNA_1_PADRE_r0.fq/DNA_1_PADRE_r0.fq_align_most_freq.counts

ric_100_aln_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_26062019173657/STRICT_NOCUSTOM/F_G_F/DNA_3_r0.fq/DNA_3_r0.fq_align.counts
ric_0_aln_counts=/home/cocca/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_26062019173657/STRICT_NOCUSTOM/F_G_F/DNA_1_PADRE_r0.fq/DNA_1_PADRE_r0.fq_align.counts
```

####################

Aggiungere un giro su riallineamento in caso non si riesca a quantificare un aplicone:
    -tolgo gli alleli con reads < 10%
    -Riallineo agli alleli più frequenti

Estrarre contig che sono a 0 nella chimera per controllare se sono tutti non informativi

```bash
fgrep -f <(cut -f 1 ~/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_24072019181711/STRICT_NOCUSTOM/F_H_L/Sc_E_0.1_r0.fq_K/Sc_E_0.1_r0.fq_align_ric_100_res.counts) DNA_3_r0.fq_align_most_freq.counts <(fgrep -f <(cut -f 1 ~/analyses/michelangelo/chimerismo/NON_HLA/26062019/test_24072019181711/STRICT_NOCUSTOM/F_H_L/Sc_E_0.1_r0.fq_K/Sc_E_0.1_r0.fq_align_ric_100_res.counts) DNA_2_MADRE_r0.fq_align_most_freq.counts | cut -f 1 -d " ")
```

---
#9/11/2020

##Some example commands


First create some folders and populate them:

```bash
samples="/home/aloisio/analyses/chimeval/TGP_POPS"
regions="/home/aloisio/analyses/chimeval/REG_TAB"

mkdir -p ${samples} ${regions}

cp /genedata/cocca/analyses/michelangelo/TGP_POPS/*.* ${samples}
cp /genedata/cocca/analyses/michelangelo/REG_TAB/*.* ${regions}

```


```bash
for maf in 0.01
do

for pop in ALL
do

for chr in {1..22}
do

TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
subset_sample="/home/aloisio/analyses/chimeval/TGP_POPS/${pop}_samples_phase3.txt"
region_list="/home/aloisio/analyses/chimeval/REG_TAB/NON_HLA_ALLARGATO.tab"


prefix_filt=`date +"%d%m%Y%H"`
out_d="/home/aloisio/analyses/chimeval/${prefix_filt}_MAF_MICHILIST/${pop}"

mkdir -p ${out_d}
echo "/home/cocca/scripts/bash_scripts/hap_inform.sh ${TGP_input} ${out_d} LIST ${region_list} ${subset_sample} ${pop}"| qsub -m ea -M michelangelo.aloisio@burlo.trieste.it -N extract_LIST_${pop}_20bl_${chr} -o ${out_d}/\$JOB_ID_${pop}_${chr}.log -e ${out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G


#extract info haps
prefix=`date +"%d%m%Y%H"`
hap_TGP_input=${out_d}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.reg_list.${pop}.vcf.gz
hap_block_list="/home/aloisio/analyses/chimeval/REG_TAB/REG_TAB/NON_HLA_ALLARGATO.tab"
hap_out_d=/home/aloisio/analyses/chimeval/${prefix}_${maf}_MAF_MICHILIST_blocks/${pop}/
mkdir -p ${hap_out_d}

for n_snp in 5
do
echo "${pop} ${chr}"
echo "/home/cocca/scripts/bash_scripts/hap_inform_extract_reads.py -i ${hap_TGP_input} -n_snp ${n_snp} -chr ${chr} -pop ${pop} -o_path ${hap_out_d} -b_mode LIST -b_list ${hap_block_list} -amp_s 200 " | qsub -m ea -M michelangelo.aloisio@burlo.trieste.it -N hap_count_${pop}_${chr} -o ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log -e ${hap_out_d}/\$JOB_ID_${pop}_${chr}.log.e -V -l h_vmem=5G -hold_jid extract_LIST_${pop}_20bl_${chr}

done
done
done
done
```