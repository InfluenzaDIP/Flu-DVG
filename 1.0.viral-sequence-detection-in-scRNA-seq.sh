makeblastdb -in ./ref/H3N2.fasta -dbtype nucl -parse_seqids -out ./blast_results/H3N2

for filename in `ls ./Rawdata/A`
do
nohup blastn -query ./blast_results/${filename%.*}.fasta -db ./blast_results/H3N2 -out ./blast_results/${filename%.*}.blast -num_threads 4 -outfmt "6 qseqid qlen sacc stitle pident qstart qend sstart send slen length evalue" > ./blast_results/${filename%.*}.log 2>&1 &
done

for filename in `ls ./Rawdata/B`
do
nohup blastn -query ./blast_results/${filename%.*}.fasta -db ./blast_results/H3N2 -out ./blast_results/${filename%.*}.blast -num_threads 4 -outfmt "6 qseqid qlen sacc stitle pident qstart qend sstart send slen length evalue" > ./blast_results/${filename%.*}.log 2>&1 &
done

for filename in `ls ./Rawdata/C`
do
nohup blastn -query ./blast_results/${filename%.*}.fasta -db ./blast_results/H3N2 -out ./blast_results/${filename%.*}.blast -num_threads 4 -outfmt "6 qseqid qlen sacc stitle pident qstart qend sstart send slen length evalue" > ./blast_results/${filename%.*}.log 2>&1 &
done

cd ./blast_results
awk -F'[ >]' '{if (NR%2==1) {print $2 > "B-1_S1_L002_R1_001.fastq.dat"}}' B-1_S1_L002_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' B-1_S1_L002_R1_001.fastq.fasta | paste -d' ' - B-1_S1_L002_R1_001.fastq.dat > temp.dat
mv temp.dat B-1_S1_L002_R1_001.fastq.dat
awk -F'[ >]' '{if (NR%2==1) {print $2 > "B-2_S1_L003_R1_001.fastq.dat"}}' B-2_S1_L003_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' B-2_S1_L003_R1_001.fastq.fasta | paste -d' ' - B-2_S1_L003_R1_001.fastq.dat > temp.dat
mv temp.dat B-2_S1_L003_R1_001.fastq.dat
awk -F'[ >]' '{if (NR%2==1) {print $2 > "B-3_S1_L003_R1_001.fastq.dat"}}' B-3_S1_L003_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' B-3_S1_L003_R1_001.fastq.fasta | paste -d' ' - B-3_S1_L003_R1_001.fastq.dat > temp.dat
mv temp.dat B-3_S1_L003_R1_001.fastq.dat

awk -F'[ >]' '{if (NR%2==1) {print $2 > "C-1_S1_L002_R1_001.fastq.dat"}}' C-1_S1_L002_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' C-1_S1_L002_R1_001.fastq.fasta | paste -d' ' - C-1_S1_L002_R1_001.fastq.dat > temp.dat
mv temp.dat C-1_S1_L002_R1_001.fastq.dat
awk -F'[ >]' '{if (NR%2==1) {print $2 > "C-2_S1_L004_R1_001.fastq.dat"}}' C-2_S1_L004_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' C-2_S1_L004_R1_001.fastq.fasta | paste -d' ' - C-2_S1_L004_R1_001.fastq.dat > temp.dat
mv temp.dat C-2_S1_L004_R1_001.fastq.dat
awk -F'[ >]' '{if (NR%2==1) {print $2 > "C-3_S1_L004_R1_001.fastq.dat"}}' C-3_S1_L004_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' C-3_S1_L004_R1_001.fastq.fasta | paste -d' ' - C-3_S1_L004_R1_001.fastq.dat > temp.dat
mv temp.dat C-3_S1_L004_R1_001.fastq.dat

awk 'NR==FNR{a[$2]=$1;next} {print $0 "\t" a[$1]}' B.dat B.blast > B.output
awk 'NR==FNR{a[$2]=$1;next} {print $0 "\t" a[$1]}' C.dat C.blast > C.output
