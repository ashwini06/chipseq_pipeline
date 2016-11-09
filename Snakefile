# vim: set ft=python:
shell.prefix("set -euo pipefail;")

configfile: "config.yaml"
workdir: "/proj/b2012025/RAW_DATA/test"

localrules: all

dir_path="/proj/b2012025/RAW_DATA/test/Rawdata/"
all_samples = config["narrow_samples"] +  config["control_samples"]
all_treat = config["narrow_samples"]
#all_samples.update(config["control_samples"])
all_fastqc = expand("02fastqc/{sam}_fastqc.zip",sam=all_samples)
all_bam = expand("03aln/{sam}.sorted.bam 03aln/{sam}.sorted.bam.bai 03aln/{sam}.sorted.bam.flagstat".split(), sam = all_samples)
all_peaks = expand("04peak/{sam}_peaks.narrowPeak", sam = all_treat)
all_annos = expand("05annotation/{sam}_annotations.txt",sam=all_treat)
all_motifs = expand("06motifs/{sam}/",sam=all_treat)
rule all:
    input : all_fastqc + all_bam + all_peaks + all_annos + all_motifs


rule fastqc:
    input: "Rawdata/{sam}.fastq"
   # input : lambda wildcards: config["narrow_samples"][wildcards.sam]
    output: "02fastqc/{sam}_fastqc.zip", "02fastqc/{sam}_fastqc.html"
    log: "00log/{sam}_fastqc.log"
    priority: 50
    message: "Running fastqc :{input}"
    shell:
       """
       module load bioinfo-tools
       module load FastQC/0.11.5
       fastqc -o 02fastqc -f fastq --noextract {input} 2> {log}
       """  

rule bowtie_align:
    input: "Rawdata/{sam}.fastq"
   # input : lambda wildcards: config["narrow_samples"][wildcards.sam]
    output: temp("03aln/{sam}.unsorted.bam")
    threads: 16 
    priority: 40
    params:
        bowtie = "--sam --best --strata -m 1 -v 3 ",
        mem = "16G"
    message: "Bowtie aligning {input}: {threads} threads/{params.mem}"
    log : "00log/{sam}_bowtiealign"
    shell :
        """
        module load bioinfo-tools bowtie/1.1.2 samtools/1.3
        bowtie {params.bowtie} {config[phred_encod]} {config[index_bowtie]} {input} 2> {log} | samtools view -Sb -F4 - > {output}
        """
  
rule sort_bam:
    input:  "03aln/{sam}.unsorted.bam"
    output: "03aln/{sam}.sorted.bam"
    log:    "00log/{sam}.sort_bam"
    threads: 12
    params:
        mem  = "12G",
    message: "sort_bam {input}: {threads} threads / {params.mem}"
    shell: 
        """
        module load bioinfo-tools samtools/1.3
        samtools sort -m 1G -@ {threads} -O bam -T {output}.tmp {input} > {output}
        """

rule index_bam:
    input:  "03aln/{sam}.sorted.bam"
    output: "03aln/{sam}.sorted.bam.bai"
    log:    "00log/{sam}.index_bam"
    threads: 1
    params:
        mem   = "500M",
    message: "index_bam {input}: {threads} threads / {params.mem}"
    shell:
        """
        module load bioinfo-tools samtools/1.3
        samtools index {input}
        """

rule flagstat_bam:
    input:  "03aln/{sam}.sorted.bam"
    output: "03aln/{sam}.sorted.bam.flagstat"
    log:    "00log/{sam}.flagstat_bam"
    threads: 1
    params:
        mem   = "500M",
    message: "flagstat_bam {input}: {threads} threads / {params.mem}"
    shell:
        """
        module load bioinfo-tools samtools/1.3
        samtools flagstat {input} > {output}
        """


#if [[ ${{#inbam[@]}} -eq 1 ]]; then 
#           #ln -s $(cd $(dirname {input}) && pwd)/$(basename {input}) {output}


rule merge_controls:
    input:   expand("03aln/{sam}.sorted.bam", sam = config["control_samples"])
    output:  "03aln/control.sorted.bam"
    log:     "00log/merge_controls"
    threads: 2
    params:
        mem = "1G",
    message: "merge_controls {input}: {threads} threads / {params.mem}"
    shell:
        """
        inbam=( {input} )
        if [[ ${{#inbam[@]}} -eq 1 ]]; then 
            cp {input} {output}
        else
            module load bioinfo-tools samtools/1.3
            samtools merge -r -@{threads} {output} {input}
        fi
        """


rule find_narrow_peaks:
    input:  "03aln/{sam}.sorted.bam", "03aln/control.sorted.bam"
    output: "04peak/{sam}_model.r", "04peak/{sam}_peaks.narrowPeak",
            "04peak/{sam}_peaks.xls", "04peak/{sam}_summits.bed",
            "04peak/{sam}_model.pdf"
    log:    "00log/{sam}.find_narrow_peaks"
    threads: 2
    params:
        mem  = "2G",
    message: "find_narrow_peaks {input}: {threads} threads / {params.mem}"
    shell:
        """
        module load bioinfo-tools MACS/2.1.0 R/3.2.3
        macs2 callpeak -t {input[0]} -c {input[1]} -f BAM -g {config[macs_g]} --outdir 04peak -n {wildcards.sam} -q 0.001 &> {log}
        cd 04peak && Rscript {wildcards.sam}_model.r
        """

rule find_broad_peaks:
    input:  "03aln/{sam}.sorted.bam", "03aln/control.sorted.bam"
    output: "04peak/{sam}_peaks.xls", "04peak/{sam}_peaks.broadPeak"
    log:    "00log/{sam}.find_broad_peaks"
    threads: 2
    params:
        mem  = "2G",
    message: "find_broad_peaks {input}: {threads} threads / {params.mem}"
    shell:
        """
        module load macs/2.1.0
        macs2 callpeak -t {input[0]} -c {input[1]} -f BAM -g {config[macs_g]} --broad --broad-cutoff 0.1 --nomodel --extsize 150 --outdir 04peak -n {wildcards.sam} -q 0.001 &> {log}
        """

rule annotate_peaks:
    input: "04peak/{sam}_peaks.narrowPeak"
    output: temp("05annotation/{sam}_annotemp.bed"), "05annotation/{sam}_annotations.txt"
    log: "00log/{sam}.annotate_peaks"
    threads: 1
    message: "annotatePeaks {input} : {threads} threads"
    shell:
        """
        echo {input}
        sed 's/^chr//g;s/^\(.\)/chr\\1/g' {input} > {output[0]}
        annotatePeaks.pl {output[0]} {config[homer_g]} > {output[1]} 2> {log}
        """

rule find_motifs:
    input: temp("05annotation/{sam}_annotemp.bed")
    output: "06motifs/{sam}/"
    log: "00log/{sam}.find_motifs"
    threads: 1
    message : "find motifs for {input}: {threads} threads"
    shell : 
        """
        findMotifsGenome.pl {input} {config[homer_g]} {output} -mask 2> {log}
        """
  
