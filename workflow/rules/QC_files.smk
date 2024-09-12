configfile : "config/config.yaml"
import os.path
rule use_fastqc:
	input:
		fasta_file = config["paths"]["data_folder"] + "/{fasta}.fastq.gz"
	output:
		fastqc_report = "{out_path}/qc_logs/fastqc/{fasta}_fastqc.html"
	threads: workflow.cores*0.3
	conda: "../envs/himer_align.yaml"
	shell:
		"fastqc -o  {wildcards.out_path}/qc_logs/fastqc -t {threads} {input.fasta_file}"

rule generate_bed_file:
	input:
		gtf_file = config["paths"]["genome_files"][config["org"]]["gtf"]
	output:
		bed_file = "{out_path}/generated.bed"
	conda: "../envs/himer_align.yaml"
	shell: "gxf2bed -i {input.gtf_file} -o {output.bed_file}"

rule rseqc_infer_experiment:
	input:
		bam_file = "{out_path}/sample{sample}_" + config["org"]  + "/star_aligned/Aligned.sortedByCoord.out.bam",
		bed_file = rules.generate_bed_file.output.bed_file
	output:
		infered = "{out_path}/qc_logs/rseqc/sample{sample}_"+  config["org"] +"/infer_experiment.txt"
	threads: workflow.cores*0.3
	conda: "../envs/himer_align.yaml"
	shell:
		"infer_experiment.py -r {input.bed_file} -i {input.bam_file} >{output.infered}"

rule index_bam:
	input:
		bam_file = "{out_path}/sample{sample}_" + config["org"]  + "/star_aligned/Aligned.sortedByCoord.out.bam"
	output:
		bai_file = "{out_path}/sample{sample}_" + config["org"]  + "/star_aligned/Aligned.sortedByCoord.out.bam.bai"
	conda: "../envs/himer_align.yaml"
	threads: workflow.cores*0.3
	shell:
		"samtools index {input.bam_file}"


rule geneBody_coverage:
	input:
		bam_file = "{out_path}/sample{sample}_" + config["org"]  + "/star_aligned/Aligned.sortedByCoord.out.bam",
		bai_file= "{out_path}/sample{sample}_" + config["org"]  + "/star_aligned/Aligned.sortedByCoord.out.bam.bai",
		bed_file = rules.generate_bed_file.output.bed_file
	output:
		gene_covr = "{out_path}/qc_logs/rseqc/sample{sample}_"+  config["org"] +"/" + "sample{sample}_"+  config["org"] + ".geneBodyCoverage.txt"
	params:
		coverage_dir = lambda wildcards, output: os.path.split(output.gene_covr)[0],
		file_prefix = lambda wildcards: wildcards.sample + "_" +  config["org"]

	conda: "../envs/himer_align.yaml"
	threads: workflow.cores*0.3
	shell:
		"geneBody_coverage.py -r {input.bed_file} -i {input.bam_file}  -o {params.coverage_dir}/{params.file_prefix}"

