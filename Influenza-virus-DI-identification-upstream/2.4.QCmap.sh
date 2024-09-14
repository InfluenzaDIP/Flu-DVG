#!/bin/bash
for each in `ls *.sam`
do 
	samtools view -hS -F 4 -q 20 $each | grep -E "@|AS:i" | grep -v "XS:" | samtools view -bS | samtools sort -O bam -o BAM/${each%.sam}_HQ.bam
	samtools depth -a BAM/${each%.sam}_HQ.bam > depth/${each%.sam}_HQ.depth
done
