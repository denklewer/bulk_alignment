configfile : "config/config.yaml"

rule align_both_bb:
        input:
                first_ref = "{out_path}/{first_org}_ref_prepared/{first_org}_genome.fa",
		second_ref = "{out_path}/{second_org}_ref_prepared/{second_org}_genome.fa",
                sample_read_1 = config["paths"]["data_folder"] + "/{sample}_R1_001.fastq.gz",
                sample_read_2 = config["paths"]["data_folder"] + "/{sample}_R2_001.fastq.gz"
        output:
                first_read_1 = "{out_path}/{sample}_{first_org}_{second_org}/{sample}_{first_org}_R1.fq",
		first_read_2 = "{out_path}/{sample}_{first_org}_{second_org}/{sample}_{first_org}_R2.fq",
                second_read_1 =  "{out_path}/{sample}_{first_org}_{second_org}/{sample}_{second_org}_R1.fq",
		second_read_2 =  "{out_path}/{sample}_{first_org}_{second_org}/{sample}_{second_org}_R2.fq"
        conda: "../envs/himer_align.yaml"
	threads: 4
        params:
                bbsplit = config["paths"]["bb_tools"] + "/bbsplit.sh -Xmx100g",
                basename = "{out_path}/{sample}_{first_org}_{second_org}/{sample}_%_R#.fq",
                path = "{out_path}/bbmap_index/{first_org}_{second_org}_%"
        shell: "{params.bbsplit} ref_{wildcards.first_org}={input.first_ref} ref_{wildcards.second_org}={input.second_ref} in={input.sample_read_1} in2={input.sample_read_2} path={params.path} basename={params.basename} ambiguous2=all"



def get_gtf_name(wildcards):
        gtf_name_param = wildcards.org + "_gtf"
        return {"orig_gtf": config["paths"][gtf_name_param]}
rule rename_chomosomes_gtf:
        input:
                unpack(get_gtf_name)
        output:
                renamed_gtf = "{out_path}/{org}_ref_prepared/{org}_genes.gtf"
        shell: """
                cat <(grep '^#' {input.orig_gtf}) <(grep -v '^#' {input.orig_gtf}  | sed 's/^/{wildcards.org}_/g') > {output.renamed_gtf}
        """


def get_ref_name(wildcards):
	ref_name_param = wildcards.org + "_ref"
	return {"orig_ref": config["paths"][ref_name_param]}

rule rename_chromosomes_fasta:
	input: 
		unpack(get_ref_name)
	output:
		renamed_fasta = "{out_path}/{org}_ref_prepared/{org}_genome.fa"
	conda: "../envs/himer_align.yaml"
        params:
                rename_reads = config["paths"]["bb_tools"] + "/rename.sh"
	shell:"{params.rename_reads} in={input.orig_ref} out={output.renamed_fasta} prefix={wildcards.org} addprefix=t addunderscore=t -Xmx20g"


def get_genome_files(wildcards):
    return expand(
                "{out_path}/{current_organism}_ref_prepared/{current_organism}_genome.fa",
                current_organism=config["orgs"], out_path= wildcards.out_path
            )
def get_gtf_files(wildcards):
    return expand(
                "{out_path}/{current_organism}_ref_prepared/{current_organism}_genes.gtf",
                current_organism=config["orgs"], out_path= wildcards.out_path
            )

rule dumb_merge_genomes:
        input:
            ref_files = get_genome_files,
            gtf_files = get_gtf_files
        output:
                chimeric_fa = "{out_path}/concatenate_ref/" + "_".join(config["orgs"]) +".fa",
                chimeric_gtf = "{out_path}/concatenate_ref/" + "_".join(config["orgs"]) + ".gtf"
        shell: """
                cat {ref_files} > {output.chimeric_fa}
                cat {gtf_files} > {output.chimeric_gtf}
                """

#not working
#rule splt_org_paired_fq:
#	input:
#		one_org_fq = rules.align_both_bb.output.first_reads   #{out_path}/{sample}_{first_org}_{second_org}/{sample}_{org}.fq
#	output:
#		fq_R1 = "{out_path}/{sample}_{first_org}_{second_org}/{org}_reads/{sample}_{org}_R1.fq",
#		fq_R2 = "{out_path}/{sample}_{first_org}_{second_org}/{org}_reads/{sample}_{org}_R2.fq"
#	conda: "../envs/himer_align.yaml"
#	params:
#                bbdemux = config["paths"]["bb_tools"] + "/reformat.sh",
#	shell: "{params.bbdemux} in={input.one_org_fq} out={wildcards.out_path}/{wildcards.sample}_{wildcards.first_org}_{wildcards.second_org}/{wildcards.org}_reads/{wildcards.sample}_{wildcards.org}_R#.fq"
