#!/usr/bin/env python3
"""

功能：从UniProt记录中提取核苷酸访问号，并下载对应DNA序列。
路径：UniProt_ID -> Nucleotide_Accession -> NCBI_FASTA
"""

import requests, time, sys, os, re, json
from typing import Optional, Tuple

# ==================== 配置 ====================
INPUT_ID_FILE = "/Bio/Project/GL_20250922_8894238869/clean_reads/extracted_protein_ids.txt"
OUTPUT_DIR = "/Bio/Project/GL_20250922_8894238869/clean_reads"
OUTPUT_FASTA = os.path.join(OUTPUT_DIR, "bacmet_real_dna.fasta")
DETAILED_LOG = os.path.join(OUTPUT_DIR, "dna_fetch_detailed.log")
FAILED_LIST = os.path.join(OUTPUT_DIR, "failed_processing.list")

# 目标核苷酸数据库类型 (在UniProt cross-references中的名称)
TARGET_NUC_DBS = ["EMBL", "GenBank", "DDBJ"]

# ==================== 核心提取函数 ====================
def get_uniprot_json(protein_id: str) -> Optional[dict]:
    """获取UniProt条目的JSON数据"""
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    try:
        resp = requests.get(url, timeout=45)
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.RequestException as e:
        print(f"  网络请求失败: {e}")
        return None

def extract_nucleotide_accession(uniprot_data: dict, pid: str) -> Optional[str]:
    """
    核心函数：从UniProt JSON中提取核苷酸访问号。
    策略：寻找交叉引用到EMBL/GenBank/DDBJ的记录，并提取其'id'。
    """
    if not uniprot_data or 'uniProtKBCrossReferences' not in uniprot_data:
        return None

    for xref in uniprot_data['uniProtKBCrossReferences']:
        db_type = xref.get('database', '')
        # 检查是否为目标核酸数据库
        if db_type in TARGET_NUC_DBS:
            nucleotide_id = xref.get('id')
            # 对ID进行清洗和基本验证
            if nucleotide_id and _looks_like_nucleotide_accession(nucleotide_id):
                return nucleotide_id
            # 如果id字段不符合，尝试从properties中寻找
            for prop in xref.get('properties', []):
                if prop.get('key') in ['SequenceId', 'Genomic_DNA', 'Accession']:
                    val = prop.get('value', '')
                    if _looks_like_nucleotide_accession(val):
                        return val
    return None

def _looks_like_nucleotide_accession(accession: str) -> bool:
    """
    启发式判断一个字符串是否为核苷酸访问号。
    规则：至少包含数字和字母，不以‘sp|’等蛋白质前缀开头，长度通常在6-20字符。
    """
    acc_clean = accession.split(';')[0].strip()  # 处理可能的分隔符
    # 排除明显的蛋白质访问号模式（包含多个‘|’或版本号在末尾）
    if re.match(r'^[A-Z]{3}\d{5}', acc_clean):  # 如‘ABC12345’格式
        return True
    if re.match(r'^[A-Z]{1,2}_?\d{5,9}(\.\d+)?$', acc_clean): # 如‘AB123456.1’或‘NC_123456’
        return True
    if re.match(r'^\w{2}\d{6}$', acc_clean): # 如‘HM362782’
        return True
    # 如果包含版本号（点），且不是蛋白质常见的模式，也接受
    if '.' in acc_clean and len(acc_clean) < 25:
        base = acc_clean.split('.')[0]
        if base.isalnum() and not base.startswith('sp|'):
            return True
    return False

def download_nucleotide_fasta(accession: str) -> Optional[str]:
    """使用核苷酸访问号从NCBI下载FASTA序列"""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        'db': 'nuccore',
        'id': accession,
        'rettype': 'fasta',
        'retmode': 'text'
    }
    try:
        resp = requests.get(url, params=params, timeout=45)
        if resp.status_code == 200:
            content = resp.text.strip()
            if content.startswith('>'):
                # 成功，返回FASTA内容
                return content
            else:
                # 可能返回了错误页面
                print(f"    NCBI返回非FASTA内容，访问号可能无效或已更新。")
                return None
    except Exception as e:
        print(f"    下载序列失败: {e}")
        return None

# ==================== 主流程 ====================
def main():
    print("="*70)
    print("BacMet蛋白质ID -> 核苷酸访问号 -> DNA序列 全自动获取")
    print("="*70)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 读取所有蛋白质ID
    try:
        with open(INPUT_ID_FILE, 'r') as f:
            all_ids = [line.strip() for line in f if line.strip()]
        total = len(all_ids)
        print(f"读取到 {total} 个待处理的蛋白质ID。")
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {INPUT_ID_FILE}")
        sys.exit(1)

    success_count = 0
    failed_details = []  # 记录失败原因

    with open(DETAILED_LOG, 'w') as log_f, open(OUTPUT_FASTA, 'w') as fasta_f:
        log_f.write(f"开始处理 {total} 个蛋白质ID\n")
        log_f.write("时间戳 | 蛋白质ID | 状态 | 核苷酸访问号 | 备注\n")
        log_f.write("-"*80 + "\n")

        for idx, pid in enumerate(all_ids, 1):
            print(f"[{idx:3d}/{total}] 处理 {pid}...")
            log_entry = {
                'time': time.strftime('%H:%M:%S'),
                'pid': pid,
                'status': '',
                'accession': '',
                'note': ''
            }

            # 步骤1：获取UniProt记录
            uniprot_data = get_uniprot_json(pid)
            if not uniprot_data:
                msg = "无法获取UniProt记录"
                print(f"  ❌ {msg}")
                log_entry['status'] = 'FAIL'
                log_entry['note'] = msg
                failed_details.append((pid, msg))
                _write_log(log_f, log_entry)
                continue

            # 步骤2：提取核苷酸访问号
            nuc_acc = extract_nucleotide_accession(uniprot_data, pid)
            if not nuc_acc:
                # 尝试输出一些调试信息到日志
                db_refs = uniprot_data.get('uniProtKBCrossReferences', [])
                db_list = [f"{x.get('database')}:{x.get('id')}" for x in db_refs[:3]]
                msg = f"未找到有效的核苷酸访问号。前几条交叉引用: {', '.join(db_list)}"
                print(f"  ❌ {msg[:60]}...")
                log_entry['status'] = 'FAIL'
                log_entry['note'] = msg
                failed_details.append((pid, "No nucleotide accession found"))
                _write_log(log_f, log_entry)
                continue

            log_entry['accession'] = nuc_acc
            print(f"  找到核苷酸访问号: {nuc_acc}")

            # 步骤3：从NCBI下载序列
            time.sleep(0.4)  # 礼貌延迟，遵守NCBI频率限制
            fasta_content = download_nucleotide_fasta(nuc_acc)
            if not fasta_content:
                msg = f"无法从NCBI下载序列 (访问号: {nuc_acc})"
                print(f"  ❌ {msg}")
                log_entry['status'] = 'FAIL'
                log_entry['note'] = msg
                failed_details.append((pid, f"Download failed for {nuc_acc}"))
                _write_log(log_f, log_entry)
                continue

            # 步骤4：保存成功的序列
            # 重写FASTA头，包含原始蛋白质ID和核苷酸访问号
            lines = fasta_content.split('\n')
            original_header = lines[0][1:]  # 去掉开头的'>'
            new_header = f">{pid} | nucleotide_accession:{nuc_acc} | {original_header}"
            fasta_f.write(new_header + "\n")
            fasta_f.write("\n".join(lines[1:]) + "\n")

            success_count += 1
            msg = f"成功，序列长度: {len(''.join(lines[1:]))} bp"
            print(f"  ✅ {msg}")
            log_entry['status'] = 'SUCCESS'
            log_entry['note'] = msg
            _write_log(log_f, log_entry)

        # 循环结束，写入最终统计
        log_f.write("\n" + "="*80 + "\n")
        log_f.write(f"处理完成于 {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_f.write(f"总计: {total}, 成功: {success_count}, 失败: {total-success_count}\n")

    # 保存失败的ID列表
    if failed_details:
        with open(FAILED_LIST, 'w') as f:
            f.write("ProteinID\tFailureReason\n")
            for pid, reason in failed_details:
                f.write(f"{pid}\t{reason}\n")

    # 打印最终报告
    print("\n" + "="*70)
    print("处理完成！")
    print(f"成功获取 {success_count} / {total} 条DNA序列")
    print(f"成功率: {success_count/total*100:.1f}%")
    print(f"\n输出文件:")
    print(f"  ✅ DNA序列 (FASTA): {OUTPUT_FASTA}")
    print(f"  📋 详细处理日志: {DETAILED_LOG}")
    if failed_details:
        print(f"  ⚠️  失败ID及原因列表: {FAILED_LIST}")
    print("\n下一步建议:")
    if success_count > 0:
        print(f"  1. 使用 {OUTPUT_FASTA} 作为Bowtie2的数据库进行丰度计算。")
    if success_count < total:
        print(f"  2. 对于失败的 {total-success_count} 个ID，可结合反向翻译的理论序列进行补充。")
    print("="*70)

def _write_log(log_file, entry):
    """格式化写入单条日志"""
    log_file.write(f"{entry['time']} | {entry['pid']} | {entry['status']:8} | "
                   f"{entry['accession']:15} | {entry['note']}\n")

if __name__ == "__main__":
    main()
