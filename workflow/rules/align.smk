configfile : "config/config.yaml"

import os.path


rule run_star:
    input:
        first_org_read_1 = config["paths"]["data_folder"] + "/{sample}_R1_001.fastq.gz",
        first_org_read_2 =  config["paths"]["data_folder"] + "/{sample}_R2_001.fastq.gz",
        sjdb =  config["paths"]["star_preprocessed_files"][config["org"]]["sjdbList"]
    output:
        sorted_bam = "{out_path}/{sample}_" + config["org"]  + "/star_aligned/Aligned.sortedByCoord.out.bam",
        sj_file = "{out_path}/{sample}_" + config["org"] + "/star_aligned/SJ.out.tab"
    conda: "../envs/himer_align.yaml"
    resources:
            mem_mb=64000
    threads: 4
    params:
        star_genome_dir=lambda wildcards, input: os.path.split(input.sjdb)[0],
        aligned_dir= lambda  wildcards, output: os.path.split(output.sorted_bam)[0]
    shell:
        """
        STAR --genomeDir {params.star_genome_dir} \
        --runThreadN {threads} \
        --readFilesIn {input.first_org_read_1} {input.first_org_read_2}\
        --outFileNamePrefix {params.aligned_dir}/  \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --readFilesCommand zcat
        """

