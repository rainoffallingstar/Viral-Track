#!/bin/bash

# 设置变量
input_dir="/public3/home/scg5695/project/15_scCRC/multi/Fastq/mRNA"
output_dir="/public3/home/scg9946/Viral-Track/umitools_output"
barcode_file="/public3/home/scg9946/Viral-Track/barcodes/barcodes.txt"
sample_list="/public3/home/scg9946/Viral-Track/barcodes/sample_list.txt"
exec 3< "$sample_list"
# 循环处理每个样本
while IFS= read -r sample_list <&3; do
        echo "正在处理样本 $sample_list 并压缩"
        # 运行 umi_tools 命令
        umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
            --stdin "${input_dir}/${sample_list}_R1_001.fastq.gz" \
            --stdout "${input_dir}/${sample_list}_R1_extracted.fastq" \
            --read2-in "${input_dir}/${sample_list}_R2_001.fastq.gz" \
            --read2-out="${output_dir}/${sample_list}_R2_extracted.fastq" \
            --log="${output_dir}/${sample_list}_extract_processed.log" \
            --whitelist="$barcode_file" 
        echo "cleaning"
        rm -rf "${input_dir}/${sample_list}_R1_extracted.fastq"
        echo "正在压缩样本"
        pigz -p 16 "${output_dir}/${sample_list}_R2_extracted.fastq"
        
done  < "$sample_list"
