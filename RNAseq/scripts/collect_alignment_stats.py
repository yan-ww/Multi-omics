from pathlib import Path

with open(snakemake.output[0], "w", encoding="utf-8") as out:
    out.write("sample\tsummary_file\n")
    for path in snakemake.input:
        sample = Path(path).name.replace("_hisat2_summary.txt", "")
        out.write(f"{sample}\t{path}\n")
