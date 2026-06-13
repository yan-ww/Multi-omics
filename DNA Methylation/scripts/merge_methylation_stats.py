#!/usr/bin/env python3
"""
合并所有样本的 Bismark 比对统计
"""

import pandas as pd
import glob
import re
from pathlib import Path

def parse_bismark_report(report_file):
    """从 Bismark 报告中提取关键统计"""
    stats = {}
    with open(report_file, 'r') as f:
        content = f.read()
        
        # 提取比对率
        mapping_rate = re.search(r'Mapping efficiency:\s+([\d.]+)%', content)
        if mapping_rate:
            stats['mapping_rate'] = float(mapping_rate.group(1))
        
        # 提取唯一比对率
        unique_rate = re.search(r'Unique alignment \(PE\):\s+([\d.]+)%', content)
        if unique_rate:
            stats['unique_rate'] = float(unique_rate.group(1))
        
        # 提取甲基化C数量
        mc_cpg = re.search(r'CpG methylated C\'s:\s+(\d+)', content)
        if mc_cpg:
            stats['methylated_CpG'] = int(mc_cpg.group(1))
        
        # 提取总C数量
        total_cpg = re.search(r'Total C\'s in CpG context:\s+(\d+)', content)
        if total_cpg:
            stats['total_CpG'] = int(total_cpg.group(1))
    
    return stats

def main():
    # 查找所有报告文件
    report_files = glob.glob("results/02_bismark_bams/*_bismark_report.txt")
    
    all_stats = []
    for report in report_files:
        sample = Path(report).stem.replace('_bismark_report', '')
        stats = parse_bismark_report(report)
        stats['sample'] = sample
        all_stats.append(stats)
    
    # 转换为 DataFrame
    df = pd.DataFrame(all_stats)
    
    # 计算甲基化水平
    if 'methylated_CpG' in df.columns and 'total_CpG' in df.columns:
        df['methylation_level'] = df['methylated_CpG'] / df['total_CpG'] * 100
    
    # 保存
    df.to_csv("results/05_reports/methylation_stats.csv", index=False)
    print(f"Saved methylation statistics for {len(df)} samples")

if __name__ == "__main__":
    main()