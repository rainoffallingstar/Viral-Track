#!/bin/bash

# 指定目录路径
directory="/public3/home/scg9946/Virus-Track/germoutput"

# 遍历目录下的文件夹
for folder in "$directory"/*; do
if [ -d "$folder" ]; then
# 获取文件夹名字作为样本名
sample_name=$(basename "$folder")

# 创建Viral_BAM_files目录
output_directory="$folder/Viral_BAM_files"
mkdir -p "$output_directory"

# 获取bam文件路径
bam_file="$folder/${sample_name}Aligned.sortedByCoord.out.bam"

# 获取chromosome文件和output文件的路径
chromosome_file="$folder/chr_target.txt"
output_file="$folder/chr_names.txt"

# 逐行读取chromosome文件和output文件的内容
while IFS= read -r chromosome && IFS= read -r output <&3; do
# 提取染色体序列并保存为BAM文件
samtools view -b "$bam_file" "$chromosome" > "$output_directory/$output"
echo "提取染色体 $chromosome 完成，保存为 $output_directory/$output"
done < "$chromosome_file" 3< "$output_file"
fi
done