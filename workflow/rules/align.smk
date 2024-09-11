configfile : "config/config.yaml"

import os.path

rule get_star_index:
    input:
        first_ref = "{out_path}/{first_org}_ref_prepared/{first_org}_genome.fa",
        second_ref = "{out_path}/{second_org}_ref_prepared/{second_org}_genome.fa",
        first_gtf = "{out_path}/{first_org}_ref_prepared/{first_org}_genes.gtf",
        second_gtf = "{out_path}/{second_org}_ref_prepared/{second_org}_genes.gtf"
    output:
        star_chr_length = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/chrLength.txt",
        star_chr_name_length = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/chrNameLength.txt",
        star_chr_name = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/chrName.txt",
        star_chr_start = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/chrStart.txt",
        star_exon_ge_tr_info = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/exonGeTrInfo.tab",
        star_exon_info = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/exonInfo.tab",
        star_gene_info = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/geneInfo.tab",
        star_genome = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/Genome",
        star_genome_parameters = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/genomeParameters.txt",
        star_sa = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/SA",
        star_sa_index = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/SAindex",
        star_sjdb_info = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/sjdbInfo.txt",
        star_sjdb_list_from_gtf_out = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/sjdbList.fromGTF.out.tab",
        star_sjdb_list_out = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/sjdbList.out.tab",
        star_transcript_info = "{out_path}/concatenate_ref/star_{first_org}_{second_org}/transcriptInfo.tab",
    priority: 2
    params:
        out_dir=lambda wildcards, output: os.path.split(output.star_genome)[0],
	sjdb_overhang=lambda wildcards: int(config["mate_length"])-1
    threads: 4
    resources:
        mem_mb=64000
    conda: "../envs/himer_align.yaml"
    shell:"""
	STAR  --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.out_dir} --genomeFastaFiles {input.first_ref} {input.second_ref} --sjdbGTFfile {input.first_gtf} {input.second_gtf} --sjdbOverhang {params.sjdb_overhang} 
	"""
	

rule run_chimeric_star:
	input:
		first_org_read_1 = "{out_path}/{sample}_" +"_".join(config["orgs"]) + "/{sample}_{org}_R1.fq",
		first_org_read_2 = "{out_path}/{sample}_" +"_".join(config["orgs"]) + "/{sample}_{org}_R2.fq",
		chimeric_sjdb =  "/scratch1/fs1/martyomov/maksim/permanent/scn-pipeline_storage_2/resources/star/mm/sjdbList.fromGTF.out.tab"
	output:
		sorted_bam = "{out_path}/{sample}_" + "_".join(config["orgs"]) + "/{org}_star_aligned/Aligned.sortedByCoord.out.bam",
		sj_file = "{out_path}/{sample}_" + "_".join(config["orgs"]) + "/{org}_star_aligned/SJ.out.tab"
	conda: "../envs/himer_align.yaml"
	resources:
        	mem_mb=64000
	threads: 4
	params:
        	star_genome_dir=lambda wildcards, input: os.path.split(input.chimeric_sjdb)[0],
		sjdb_overhang=lambda wildcards: int(config["mate_length"])-1
	shell:
		"""
		STAR --genomeDir {params.star_genome_dir} \
		--runThreadN {threads} \
		--readFilesIn {input.first_org_read_1} {input.first_org_read_2}\
		""" + \
        "--outFileNamePrefix {wildcards.out_path}/{wildcards.sample}_" +"_".join(config["orgs"])  + " /{wildcards.org}_star_aligned/ " + \
        """
        \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard \
		--sjdbFileChrStartEnd {input.chimeric_sjdb}\
		--sjdbOverhang {params.sjdb_overhang} \
		"""
		
		
		
