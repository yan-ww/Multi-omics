"""
修剪模块
使用 fastp 进行质量修剪和适配器去除
"""

rule fastp:
    """使用 fastp 进行质量修剪"""
    input:
        r1=lambda wildcards: get_fastq_path(wildcards.sample, 1),
        r2=lambda wildcards: get_fastq_path(wildcards.sample, 2) if config.get('paired_end') else None
    output:
        r1=get_trimmed_path("{sample}", 1),
        r2=get_trimmed_path("{sample}", 2) if config.get('paired_end') else None,
        json="results/trimmed/{sample}_fastp.json",
        html="results/trimmed/{sample}_fastp.html"
    params:
        quality=config['trimming']['quality'],
        length=config['trimming']['min_length'],
        n_base_limit=config['trimming']['n_base_limit'],
        threads=config['threads']['fastp']
    log:
        "logs/fastp/{sample}.log"
    conda:
        "../envs/fastp.yaml"
    run:
        cmd = ["fastp",
               f"--thread {params.threads}",
               f"--qualified_quality_phred {params.quality}",
               f"--length_required {params.length}",
               f"--n_base_limit {params.n_base_limit}",
               f"--json {output.json}",
               f"--html {output.html}",
               f"--report_title {wildcards.sample}"]
        
        # 双端或单端模式
        if config.get('paired_end'):
            cmd.extend([f"--in1 {input.r1}",
                       f"--in2 {input.r2}",
                       f"--out1 {output.r1}",
                       f"--out2 {output.r2}"])
        else:
            cmd.extend([f"--in {input.r1}",
                       f"--out {output.r1}"])
        
        # 添加适配器选项
        if config.get('trimming', {}).get('adapter_trimming', True):
            cmd.append("--detect_adapter_for_pe")
        
        # 运行命令
        shell(" ".join(cmd) + f" > {log} 2>&1")