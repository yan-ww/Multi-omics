"""
定量模块
使用 featureCounts 进行基因水平计数
"""

rule featurecounts:
    """使用 featureCounts 进行基因计数"""
    input:
        bam="results/hisat2/{sample}.bam",
        annotation=config['reference']['gtf_file']
    output:
        counts="results/featurecounts/{sample}_counts.txt",
        summary="results/featurecounts/{sample}_counts.txt.summary"
    params:
        feature_type=config['featurecounts']['feature_type'],
        attribute_type=config['featurecounts']['attribute_type'],
        strand=config['featurecounts']['strand'],
        threads=config['threads']['featurecounts'],
        min_overlap=config.get('featurecounts', {}).get('min_overlap', 1)
    log:
        "logs/featurecounts/{sample}.log"
    conda:
        "../envs/featurecounts.yaml"
    run:
        cmd = ["featureCounts",
               f"-T {params.threads}",
               f"-t {params.feature_type}",
               f"-g {params.attribute_type}",
               f"-s {params.strand}",
               f"-a {input.annotation}",
               f"-o {output.counts}",
               f"--minOverlap {params.min_overlap}"]
        
        # 双端模式添加 -p 参数
        if config.get('paired_end'):
            cmd.append("-p")
        
        # 添加BAM文件
        cmd.append(input.bam)
        
        shell(" ".join(cmd) + f" > {log} 2>&1")

rule merge_counts:
    """合并所有样本的计数矩阵"""
    input:
        expand("results/featurecounts/{sample}_counts.txt", sample=SAMPLES)
    output:
        "results/featurecounts/gene_counts.txt"
    params:
        samples=SAMPLES
    log:
        "logs/featurecounts/merge_counts.log"
    conda:
        "../envs/featurecounts.yaml"
    run:
        import pandas as pd
        import numpy as np
        
        # 读取第一个文件
        first_file = input[0]
        counts_df = pd.read_csv(first_file, sep='\t', comment='#', 
                                header=0, index_col=0)
        
        # 提取计数列（第7列，索引6）
        counts_matrix = counts_df.iloc[:, [6]]
        counts_matrix.columns = [SAMPLES[0]]
        
        # 合并其他样本
        for i, file in enumerate(input[1:], 1):
            sample_counts = pd.read_csv(file, sep='\t', comment='#', 
                                       header=0, index_col=0)
            counts_matrix[SAMPLES[i]] = sample_counts.iloc[:, 6]
        
        # 添加基因信息
        gene_info = counts_df.iloc[:, :6]
        final_matrix = pd.concat([gene_info, counts_matrix], axis=1)
        
        # 保存合并的计数矩阵
        final_matrix.to_csv(output[0], sep='\t')
        
        print(f"合并完成: {len(SAMPLES)} 个样本, {len(final_matrix)} 个基因")