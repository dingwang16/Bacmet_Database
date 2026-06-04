import re

# 读取文件
with open('bacmet_protein_ids.txt', 'r') as f:
    content = f.read()

# 使用正则表达式匹配两种格式的ID
# 匹配 tr|ID| 格式 (如: tr|Q5FAM9|)
# 匹配 sp|ID| 格式 (如: sp|P0AE06|)
# (?:tr|sp) 匹配 tr 或 sp，但不捕获这个组
# \|([A-Z0-9_]+)\| 捕获两个竖线之间的ID
pattern = r'(?:tr|sp)\|([A-Z0-9_]+)\|'
extracted_ids = re.findall(pattern, content)

# 去重并保存
unique_ids = list(set(extracted_ids))

# 输出结果
print(f"找到 {len(extracted_ids)} 个ID（含重复）")
print(f"去重后得到 {len(unique_ids)} 个唯一ID")

# 保存到文件
with open('extracted_protein_ids.txt', 'w') as f:
    for pid in unique_ids:
        f.write(pid + '\n')

print(f"已提取到 {len(unique_ids)} 个UniProt/RefSeq格式的ID，并保存到 'extracted_protein_ids.txt'")