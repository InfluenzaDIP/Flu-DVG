#!/bin/bash
cd Ref/InfluenzaRef/H3N2
Ref=(`ls *.fa`)
cd ../../../
pwd
for each in ${Ref[*]}
do
	echo $each
	ref=${each%.fa}
	bowtie2 -p 16 -x Ref/InfluenzaRef/H3N2/$ref -1 3.QC/MatePaired/H3N2_1_mt.fastq -2 3.QC/MatePaired/H3N2_2_mt.fastq -S 4.SAM/${ref}.sam
done 
