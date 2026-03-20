"""
比对模块
使用 HISAT2 进行比对，samtools 进行排序和索引
"""

rule hisat2_index:
    """构建 HISAT2 索引（如果不存在）"""
    input:
        fasta=config['reference']['genome_fasta'],
        gtf=config['reference']['gtf_file']
    output:
        idx_dir=directory(config['reference']['hisat2_index'])
    params:
        threads=config['threads']['hisat2_index'],
        splice_sites=config.get('hisat2', {}).get('splice_sites', '')
    log:
        "logs/hisat2/build_index.log"
    conda:
        "../envs/hisat2.yaml"
    run:
        # 检查索引是否已存在
        idx_files = list(Path(output.idx_dir).glob("*.ht2"))
        if idx_files:
            print(f"HISAT2索引已存在: {output.idx_dir}")
        else:
            shell("hisat2-build -p {params.threads} "
                  "{input.fasta} "
                  "{output.idx_dir}/index > {log} 2>&1")
            
            # 如果有GTF文件，提取剪接位点
            if params.splice_sites:
                shell("hisat2_extract_splice_sites.py {input.gtf} "
                      "> {output.idx_dir}/splice_sites.txt")

rule hisat2_align:
    """使用 HISAT2 进行比对"""
    input:
        index=config['reference']['hisat2_index'] + "/index",
        r1=get_trimmed_path("{sample}", 1),
        r2=get_trimmed_path("{sample}", 2) if config.get('paired_end') else None
    output:
        sam="results/hisat2/{sample}.sam",
        summary="results/hisat2/{sample}_summary.txt"
    params:
        threads=config['threads']['hisat2_align'],
        options=config.get('hisat2', {}).get('align_options', '--rna-strandness RF')
    log:
        "logs/hisat2/{sample}_align.log"
    conda:
        "../envs/hisat2.yaml"
    run:
        cmd = ["hisat2",
               f"-p {params.threads}",
               f"-x {input.index}",
               params.options,
               f"-S {output.sam}"]
        
        if config.get('paired_end'):
            cmd.extend([f"-1 {input.r1}",
                       f"-2 {input.r2}"])
        else:
            cmd.append(f"-U {input.r1}")
        
        # 添加摘要输出
        cmd.append(f"2> {output.summary}")
        
        shell(" ".join(cmd) + f" > {log} 2>&1")

rule samtools_sort:
    """排序并转换SAM到BAM"""
    input:
        "results/hisat2/{sample}.sam"
    output:
        bam="results/hisat2/{sample}.bam",
        bai="results/hisat2/{sample}.bam.bai"
    params:
        threads=config['threads']['samtools'],
        memory=config.get('samtools', {}).get('memory', '2G')
    log:
        "logs/hisat2/{sample}_sort.log"
    conda:
        "../envs/samtools.yaml"
    run:
        shell("samtools sort -@ {params.threads} "
              "-m {params.memory} "
              "-o {output.bam} {input} > {log} 2>&1")
        shell("samtools index -@ {params.threads} "
              "{output.bam} >> {log} 2>&1")

rule samtools_stats:
    """生成比对统计信息"""
    input:
        "results/hisat2/{sample}.bam"
    output:
        "results/hisat2/{sample}_stats.txt"
    params:
        threads=config['threads']['samtools']
    log:
        "logs/hisat2/{sample}_stats.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats -@ {params.threads} {input} > {output} 2> {log}"