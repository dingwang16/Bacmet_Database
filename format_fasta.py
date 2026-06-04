#!/usr/bin/env python3
"""
文件名: format_fasta.py
功能：标准化FASTA格式，将连续序列按固定行宽重新排版。
"""
import sys
import os

input_fasta = "/Bio/Project/GL_20250922_8894238869/clean_reads/bacmet_real_dna.fasta"
output_fasta = "/Bio/Project/GL_20250922_8894238869/clean_reads/bacmet_real_dna_formatted.fasta"
line_width = 60  # 标准FASTA每行字符数

print(f"正在格式化FASTA文件: {input_fasta}")
print(f"输出文件: {output_fasta}")
print(f"每行宽度: {line_width} 个字符")

try:
    with open(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
        header = ''
        seq_buffer = []
        seq_count = 0

        for line in fin:
            line = line.rstrip()
            if line.startswith('>'):
                # 如果之前有序列在缓存中，先写入它
                if seq_buffer:
                    # 将缓存的序列字符连接并分割为固定行宽
                    full_seq = ''.join(seq_buffer)
                    for i in range(0, len(full_seq), line_width):
                        fout.write(full_seq[i:i + line_width] + '\n')
                    seq_buffer = []

                # 写入新序列的头行
                fout.write(line + '\n')
                header = line
                seq_count += 1
            else:
                # 收集序列行，移除可能存在的内部空格或数字
                seq_buffer.append(line)

        # 处理最后一个序列
        if seq_buffer:
            full_seq = ''.join(seq_buffer)
            for i in range(0, len(full_seq), line_width):
                fout.write(full_seq[i:i + line_width] + '\n')

    print(f"格式化完成！共处理 {seq_count} 条序列。")
    print(f"格式化后的文件已保存为: {output_fasta}")

except FileNotFoundError:
    print(f"错误：找不到输入文件 {input_fasta}")
    sys.exit(1)
except Exception as e:
    print(f"处理过程中发生错误: {e}")
    sys.exit(1)