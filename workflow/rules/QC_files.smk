configfile : "config/config.yaml"

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
		gtf_file = rules.dumb_merge_genomes.output.chimeric_gtf
	output:
		bed_file = "{out_path}/concatenate_ref/{first_org}_{second_org}_generated.bed"
	conda: "../envs/himer_align.yaml"
	shell: "gxf2bed -i {input.gtf_file} -o {output.bed_file}"

rule rseqc_infer_experiment:
	input:
		bam_file = "{out_path}/{sample}_{first_org}_{second_org}/{org}_star_aligned/Aligned.sortedByCoord.out.bam",
		bed_file = rules.generate_bed_file.output.bed_file
	output:
		infered = "{out_path}/qc_logs/rseqc/{sample}_{first_org}_{second_org}/{org}_infer_experiment.txt"
	threads: workflow.cores*0.3
	conda: "../envs/himer_align.yaml"
	shell:
		"infer_experiment.py -r {input.bed_file} -i {input.bam_file} >{output.infered}"


rule geneBody_coverage:
	input:
		bam_file = "{out_path}/{sample}_{first_org}_{second_org}/{org}_star_aligned/Aligned.sortedByCoord.out.bam",
		bed_file = rules.generate_bed_file.output.bed_file
	output:
		gene_covr = "{out_path}/qc_logs/rseqc/{sample}_{first_org}_{second_org}/{org}.geneBodyCoverage.txt"
	conda: "../envs/himer_align.yaml"
	threads: workflow.cores*0.3
	shell:
		"geneBody_coverage.py -r {input.bed_file} -i {input.bam_file}  -o {wildcards.out_path}/qc_logs/rseqc/{wildcards.sample}_{wildcards.first_org}_{wildcards.second_org}/{wildcards.org}"

