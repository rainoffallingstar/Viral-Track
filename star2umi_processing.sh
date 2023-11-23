# generate genomic files for virus

STAR --runThreadN 32 --runMode genomeGenerate 
--genomeDir virushg38/
--genomeFastaFiles genomes.fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa

# UMI-tools processing

# 文件夹权限

# step 1

cat /public3/home/scg5695/project/15_scCRC/multi/Fastq/mRNA/*_R1_001.fastq.gz > crcsinglecell_R1.fastq.gz

cat /public3/home/scg5695/project/15_scCRC/multi/Fastq/mRNA/*_R2_001.fastq.gz > crcsinglecell_R2.fastq.gz

# umi bc
# 哪种测序方式？一个样本大概多少细胞？--》 10X，可以自动估计细胞

# step 2
# about bc-pattern
zcat crcsinglecell_R1.fastq.gz | head -n2

umi_tools whitelist --stdin crcsinglecell_R1.fastq.gz --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --log2stderr > whitelist.txt
                     
# step 3
# by using .fastq as output to avoid O/I bottleneck

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin crcsinglecell_R1.fastq.gz \
                  --stdout crcsinglecell_R1_extracted.fastq \
                  --read2-in crcsinglecell_R2.fastq.gz \
                  --read2-out=crcsinglecell_R2_extracted.fastq \
                  --whitelist=whitelist.txt 





                     
                     
                     

