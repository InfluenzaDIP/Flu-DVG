makeblastdb -in /workdir/wangph/ref/H3N2.fasta -dbtype nucl -parse_seqids -out ./blast/H3N2

for filename in `ls /workdir/licd/Jida/X101SC22116113-Z01-J001/Rawdata/A`
do
# gunzip -c /workdir/licd/Jida/X101SC22116113-Z01-J001/Rawdata/A/$filename > ./blast/${filename%.*}
# sed -n '1~4s/^@/>/p;2~4p' ./blast/${filename%.*} > ./blast/${filename%.*}.fasta
nohup blastn -query ./blast/${filename%.*}.fasta -db ./blast/H3N2 -out ./blast/${filename%.*}.blast -num_threads 4 -outfmt "6 qseqid qlen sacc stitle pident qstart qend sstart send slen length evalue" > ./blast/${filename%.*}.log 2>&1 &
# rm ./blast/${filename%.*}
done

for filename in `ls /workdir/licd/Jida/X101SC22116113-Z01-J001/Rawdata/B`
do
# gunzip -c /workdir/licd/Jida/X101SC22116113-Z01-J001/Rawdata/B/$filename > ./blast/${filename%.*}
# sed -n '1~4s/^@/>/p;2~4p' ./blast/${filename%.*} > ./blast/${filename%.*}.fasta
nohup blastn -query ./blast/${filename%.*}.fasta -db ./blast/H3N2 -out ./blast/${filename%.*}.blast -num_threads 4 -outfmt "6 qseqid qlen sacc stitle pident qstart qend sstart send slen length evalue" > ./blast/${filename%.*}.log 2>&1 &
# rm ./blast/${filename%.*}
done

for filename in `ls /workdir/licd/Jida/X101SC22116113-Z01-J001/Rawdata/C`
do
# gunzip -c /workdir/licd/Jida/X101SC22116113-Z01-J001/Rawdata/C/$filename > ./blast/${filename%.*}
# sed -n '1~4s/^@/>/p;2~4p' ./blast/${filename%.*} > ./blast/${filename%.*}.fasta
nohup blastn -query ./blast/${filename%.*}.fasta -db ./blast/H3N2 -out ./blast/${filename%.*}.blast -num_threads 4 -outfmt "6 qseqid qlen sacc stitle pident qstart qend sstart send slen length evalue" > ./blast/${filename%.*}.log 2>&1 &
# rm ./blast/${filename%.*}
done

for filename in `ls /workdir/licd/Jida/X101SC22116113-Z01-J001/Rawdata/D`
do
# gunzip -c /workdir/licd/Jida/X101SC22116113-Z01-J001/Rawdata/D/$filename > ./blast/${filename%.*}
# sed -n '1~4s/^@/>/p;2~4p' ./blast/${filename%.*} > ./blast/${filename%.*}.fasta
nohup blastn -query ./blast/${filename%.*}.fasta -db ./blast/H3N2 -out ./blast/${filename%.*}.blast -num_threads 4 -outfmt "6 qseqid qlen sacc stitle pident qstart qend sstart send slen length evalue" > ./blast/${filename%.*}.log 2>&1 &
# rm ./blast/${filename%.*}
done

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

awk -F'[ >]' '{if (NR%2==1) {print $2 > "D-1_S1_L002_R1_001.fastq.dat"}}' D-1_S1_L002_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' D-1_S1_L002_R1_001.fastq.fasta | paste -d' ' - D-1_S1_L002_R1_001.fastq.dat > temp.dat
mv temp.dat D-1_S1_L002_R1_001.fastq.dat
awk -F'[ >]' '{if (NR%2==1) {print $2 > "D-2_S1_L004_R1_001.fastq.dat"}}' D-2_S1_L004_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' D-2_S1_L004_R1_001.fastq.fasta | paste -d' ' - D-2_S1_L004_R1_001.fastq.dat > temp.dat
mv temp.dat D-2_S1_L004_R1_001.fastq.dat
awk -F'[ >]' '{if (NR%2==1) {print $2 > "D-3_S1_L003_R1_001.fastq.dat"}}' D-3_S1_L003_R1_001.fastq.fasta
awk 'NR % 2 == 0 {print substr($0, 1, 16)}' D-3_S1_L003_R1_001.fastq.fasta | paste -d' ' - D-3_S1_L003_R1_001.fastq.dat > temp.dat
mv temp.dat D-3_S1_L003_R1_001.fastq.dat

awk 'NR==FNR{a[$2]=$1;next} {print $0 "\t" a[$1]}' B.dat B.blast > B.output
awk 'NR==FNR{a[$2]=$1;next} {print $0 "\t" a[$1]}' C.dat C.blast > C.output
awk 'NR==FNR{a[$2]=$1;next} {print $0 "\t" a[$1]}' D.dat D.blast > D.output