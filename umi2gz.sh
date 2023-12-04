#!/bin/bash

# 设置变量
input_dir="/public3/home/scg5695/project/15_scCRC/multi/Fastq/mRNA"
output_dir="/public3/home/scg9946/Viral-Track/umitools_output"
barcode_file="/public3/home/scg9946/Viral-Track/barcodes.txt"
sample_list="/public3/home/scg9946/Viral-Track/sample_list.txt"

# 循环处理每个样本
while IFS= read -r sample_list <&3; do
        echo "正在处理样本 $R1sample_list 并压缩"
        # 运行 umi_tools 命令
        umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
            --stdin "${input_dir}/${sample_list}_R1_001.fastq.gz" \
            --read2-in "${input_dir}/${sample_list}_R2_001.fastq.gz" \
            --read2-stdout="${output_dir}/${sample_name}_R2_extracted.fastq" \
            --log="${output_dir}/${sample_name}_extract_processed.log" \
            --whitelist="$barcode_file" 
        echo "正在压缩样本"
        pigz -p 30 "${output_dir}/${sample_name}_R2_extracted.fastq"
        
done  < "$sample_list"