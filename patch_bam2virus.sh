#!/bin/bash

# 染色体名称和输出文件名的文本文件
chromosome_file="target_chr.txt"
output_file="target_chrname.txt"

# BAM文件路径
bam_file="/public3/home/scg9946/Viral-Track/virus_output/crcsinglecell_R2_extracted/crcsinglecell_R2_extracted_Aligned.sortedByCoord.out.bam"

# 逐行读取染色体名称和输出文件名
while IFS= read -r chromosome && IFS= read -r output <&3; do
    # 提取染色体序列并保存为BAM文件
    samtools view -b "$bam_file" "$chromosome" > "$output"
    echo "提取染色体 $chromosome 完成，保存为 $output"
done < "$chromosome_file" 3< "$output_file"