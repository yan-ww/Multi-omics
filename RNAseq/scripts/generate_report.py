from pathlib import Path
from html import escape

alignment = Path(snakemake.input.alignment_stats).read_text(encoding="utf-8")
de_results = "\n".join(Path(path).name for path in snakemake.input.de_results)
multiqc = "qc/" + Path(snakemake.input.multiqc).name

html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>RNA-seq report</title></head>
<body><h1>Bulk RNA-seq report</h1>
<p>MultiQC report: <a href="{escape(multiqc)}">{escape(multiqc)}</a></p>
<h2>Alignment summaries</h2><pre>{escape(alignment)}</pre>
<h2>Differential-expression results</h2><pre>{escape(de_results)}</pre>
</body></html>
"""
Path(snakemake.output[0]).write_text(html, encoding="utf-8")
