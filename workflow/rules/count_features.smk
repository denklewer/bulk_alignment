configfile : "config/config.yaml"

rule featureCounts:
    input:
        gtf_file = rules.dumb_merge_genomes.output.chimeric_gtf,
        sorted_bam = rules.run_chimeric_star.output.sorted_bam #"{out_path}/{sample}_ +  "_".join(config["orgs"]) +"/star_aligned/Aligned.sortedByCoord.out.bam"
    output:
        counts = "{out_path}/{sample}_" + "_".join(config["orgs"]) +"/counts/counts.txt",
        summary = "{out_path}/{sample}_"+ "_".join(config["orgs"]) +"/counts/counts.txt.summary"
    threads: workflow.cores
    conda: "../envs/himer_align.yaml"
    shell:
 #       "featureCounts -p --extraAttributes gene_name -t gene -g gene_id  -T {threads}  -a {input.gtf_file} -o {output.counts} {input.sorted_bam} --byReadGroup"
        "featureCounts -p --extraAttributes gene_name  -g gene_id  -T {threads}  -a {input.gtf_file} -o {output.counts} {input.sorted_bam}"
