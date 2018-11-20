#vim:set ft=python
shell.prefix("set -euo pipefail;")

configfile:"/proj/uppstore2017048/P10161/01-Results/scripts/config.yaml"
workdir:"/proj/uppstore2017048/P10161/01-Results"

localrules: all

#CASES=config["broad_peaks"]["D3.5_test"]
#CONTROLS=config["control_peaks"]["D3.5_test"]

## D3.5 = WT:series 3,4, Foxa2:series=1,2
# D5.5 = WT:series 1,2 , Foxa2:series=1,2

CASES=config["broad_peaks"]["mark"]
CONTROLS=config["control_peaks"]["mark"]
DAY=config["day"]
COND=config["cond"]
SERIES=config["series"]
useRep=config["deepRep"]
deep_flnm=config["deeptools"]["filename"]
markers=config["markers"]
enhancers=config["enhancers"][config['day'][0]]

ALL_SAMPLES = CASES + CONTROLS 

single_nm = expand("{day}_{cond}_{all_sam}-{series}",all_sam=ALL_SAMPLES,cond=COND,day=DAY,series=SERIES)

mixed_nm = expand("{day}_{cond}_{ip}-{series}_vs_{day}_{cond}_{input}-{series}",day=DAY,cond=COND,ip=CASES,input=CONTROLS,series=SERIES)

marker_nm = expand("{day}_{cond}_{mark}",day=DAY,cond=COND,mark=CASES)


ALL_bamCompare = expand("deepTools_ashwini/bamCompare/{sample}_bamCompare.bw", sample=mixed_nm)
ALL_bamCoverage = expand("deepTools_ashwini/bamCoverage/{sample}_bamCoverage.bw", sample=single_nm)

ALL_broadPEAKS = expand("macs/{nam}_peaks.xls".split(), nam=mixed_nm)
ALL_sicerPEAKS = expand("SICER/{nam}_epic.csv SICER/{nam}_epic.bed".split(), nam=mixed_nm)

ALL_peakAnno = expand("peakanno/{pk_nm}_filtered.txt peakanno/{pk_nm}_filtered.bed peakanno/{pk_nm}_filtered_annotated.txt peakanno/SICER/{pk_nm}_filtered_epic.txt peakanno/SICER/{pk_nm}_filtered_epic.bed peakanno/SICER/{pk_nm}_filtered_annotated_epic.txt".split(),pk_nm=marker_nm)

ALL_bigBED = expand("bigBED/{pk_nm}_filtered.bb bigBED/SICER/{pk_nm}_filtered_epic.bb".split(),pk_nm=marker_nm)


ALL_deepPROFILES = expand("deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_plotProfile.png deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_clusters.bed deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_heatmap1.png deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_heatmap2.png".split(),day=DAY,cond=COND,rep=useRep,nam=deep_flnm)

ALL_clusterANNOTATE = expand("deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_clusters_annotated.bed",day=DAY,cond=COND,rep=useRep,nam=deep_flnm)


ALL_genePROFILES = expand("deepTools_ashwini/profiles/{day}_{cond}-{rep}/markers/{day}_{cond}-{rep}_{genename}_TSS_plotProfile.png",day=DAY,cond=COND,rep=useRep,genename=markers)

ALL_enhancerPROFILES = expand("deepTools_ashwini/profiles/{day}_{cond}-{rep}/enhancers/{day}_{cond}-{rep}_{genename}_FoxA2peaks_plotProfile.png",day=DAY,cond=COND,rep=useRep,genename=enhancers)

rule all:
#     input: ALL_bamCompare + ALL_bamCoverage + ALL_broadPEAKS + ALL_peakAnno +  ALL_sicerPEAKS + ALL_peakAnno + ALL_bigBED + ALL_deepPROFILES + ALL_clusterANNOTATE
    input: ALL_enhancerPROFILES 


rule bamCompare:
    input:
        case = "/proj/uppstore2017048/P10161/01-Results/picard/{case_id}.dedup.sorted.bam", control="/proj/uppstore2017048/P10161/01-Results/picard/{control_id}.dedup.sorted.bam"
    output:
        bamCom = "deepTools_ashwini/bamCompare/{case_id}_vs_{control_id}_bamCompare.bw"
    log : "00log/{case_id}_vs_{control_id}_bamCompare.log"
    threads: 18
    priority: 1
    message: "generating bigwig files - bamCompare for {input}"
    shell:
        """
        module load bioinfo-tools deepTools/3.1.0
        #bamCompare -b1 {input.case} -b2 {input.control} --binSize 25 --scaleFactorsMethod SES --operation subtract --sampleLength 6000 --numberOfSamples 500000 -o {output.bamCom} &> {log}
        #--extendReads --effectiveGenomeSize 2652783500 --ignoreDuplicates 
        bamCompare -b1 {input.case} -b2 {input.control} --binSize 25 --scaleFactorsMethod "None" --operation subtract  --normalizeUsing RPKM -o {output.bamCom} -p "max" --ignoreDuplicates &> {log}
        """


rule bamCoverage:
    input: "picard/{sample}.dedup.sorted.bam"
    output : "deepTools_ashwini/bamCoverage/{sample}_bamCoverage.bw"
    threads: 12
    priority: 20
    log : "00log/{sample}_bamCoverage.log"
    message : "generating bigwig files - bamcoverage for {input}"
    shell :
        """
        module load bioinfo-tools deepTools/3.1.0
        bamCoverage -b {input} --normalizeUsing RPKM --binSize 25 --extendReads 200 -o {output} &> {log}
        """


rule macs2_broadPeaks:
    input: control = "picard/{control_id}.dedup.sorted.bam", case="picard/{case_id}.dedup.sorted.bam"
    #output: "macs2/{case_id}_vs_{control_id}_peaks.broadPeak", "macs2/{case_id}_vs_{control_id}_peaks.xls", "macs2/{case_id}_vs_{control_id}_model.r"
    output: "macs/{case_id}_vs_{control_id}_peaks.xls"
    log: "00log/{case_id}_vs_{control_id}_call_broadpeaks.log"
    threads: 12
    priority : 15
    params:
        broad_nm = "{case_id}_vs_{control_id}"
    message: "call_peaks macs2 broad {input}: {threads} threads"
    shell:
        """
        module load bioinfo-tools MACS/2.1.0
        macs2 callpeak -t {input.case} -c {input.control} -f BAM -g mm --outdir macs -n {params.broad_nm} --broad -q 0.01 &> {log}
        """


rule sicer_epicPeaks:
    input: control = "picard/{control_id}.dedup.sorted.bed", case="picard/{case_id}.dedup.sorted.bed"
    output: csv = "SICER/{case_id}_vs_{control_id}_epic.csv",
            bed = "SICER/{case_id}_vs_{control_id}_epic.bed",
    params:
        chrom_sizes="~/mm10.chrom.sizes",
        genome_size=0.8
    threads: 24
    priority:14
    shell:
        """
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}}; {{print "chr"$1,$2,$3,"U0","0",$6}}' {input.case} > SICER/case_temp.bed
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}}; {{print "chr"$1,$2,$3,"U0","0",$6}}' {input.control} > SICER/control_temp.bed
        epic --number-cores {threads} --genome mm10 --gaps-allowed 3 --window-size 200 --fragment-size 200 --bed {output.bed} --treatment SICER/case_temp.bed --control SICER/control_temp.bed --effective-genome-fraction {params.genome_size} --chromsizes {params.chrom_sizes}  > {output.csv}
        rm SICER/case_temp.bed SICER/control_temp.bed
       """


import glob
def macs2_peakfls(wildcards):
    return glob.glob ( "/proj/uppstore2017048/P10161/01-Results/" + "macs/" + wildcards.pk_nm + "*" + "_vs_" + "*" + "_peaks.xls" )  


rule macs2_peakOverlaps_annotate:
    """
    find the replicate files from macs2 for each mark
    """
    input: macs2_peakfls
    output:
        txt="peakanno/{pk_nm}_filtered.txt",
        bed="peakanno/{pk_nm}_filtered.bed",
        annofl="peakanno/{pk_nm}_filtered_annotated.txt"
    params:
        rlocate = "/pica/h1/ashwini/anaconda2/lib/R/library",
        REF_macs = "mm",
        filtering = "~/mm10.blacklist.bed",
       # gtf = "/sw/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
        gtf = "/proj/uppstore2017048/P9253/02-Results/GRCm38_ensembl.gtf"
    threads : 18
    priority: 10
    message : "Annotating macs peaks for {input}"
    shell :
        """
        Rscript ~/scripts/chipseq_pipeline/post_peak_calling_processing.v3.r {params.rlocate} {params.REF_macs} {params.filtering} {params.gtf} {output.txt} {output.bed} {output.annofl} {input}
        """

import glob
def sicer_peakfls(wildcards):
    return glob.glob ( "/proj/uppstore2017048/P10161/01-Results/" + "SICER/" + wildcards.pk_nm + "*" + "_vs_" + "*" + "_epic.bed" )


rule sicer_peakOverlaps_annotate:
    """
    find the replicate files from sicer for each mark
    """
    input: sicer_peakfls
    output:
        txt="peakanno/SICER/{pk_nm}_filtered_epic.txt",
        bed="peakanno/SICER/{pk_nm}_filtered_epic.bed",
        annofl="peakanno/SICER/{pk_nm}_filtered_annotated_epic.txt"
    params:
        rlocate = "/pica/h1/ashwini/anaconda2/lib/R/library",
        REF_macs = "mm",
        filtering = "~/mm10.blacklist.bed",
        gtf = "/proj/uppstore2017048/P9253/02-Results/GRCm38_ensembl.gtf"
    threads : 12
    priority: 10
    message : "Finding SICER peak replicate overlaps and annotating for {input}"
    shell :
        """
        Rscript ~/scripts/chipseq_pipeline/SICER_post_peak_calling_processing.r {params.rlocate} {params.REF_macs} {params.filtering} {params.gtf} {output.txt} {output.bed} {output.annofl} {input}
        """


### check once input bed file, whether to add prefix "chr" or not and change awk line accordingly
rule bed2bigBED:
    input: "peakanno/SICER/{pk_nm}_filtered_epic.bed"
    output: "bigBED/SICER/{pk_nm}_filtered_epic.bb"
    threads: 10
    log:"00log/{pk_nm}_bed2bB.log"
    message: "converting bed files to bigBed"
    shell:
        """
        module load bioinfo-tools ucsc-utilities/v345
        awk '{{if ($2 != '-1') {{print $1"\\t"$2"\\t"$3}}}}' {input} | sort -k1,1 -k2,2n > bigBED/temp_file
        bedToBigBed bigBED/temp_file ~/mm10.chrom.sizes {output} &> {log}
        rm bigBED/temp_file
       """



import glob
import fnmatch
import os 
import re

def fls4deep(wildcards):
    return glob.glob("/proj/uppstore2017048/P10161/01-Results/" + "deepTools_ashwini/bamCompare/" + wildcards.day + "_" + wildcards.cond + "_" +  "*" + "-" + wildcards.rep + "*" + ".bw" )
    


def samplenams(wildcards):
    names=[]
    for i in fls4deep(wildcards):
        sample_nm=os.path.basename(i)
        nam=re.sub('_vs_.*bw','',sample_nm)
        if re.search(r'FoxA2mut',nam):
            nam=re.sub('FoxA2mut','Fmut',nam)
        names.append(nam)
    return names



rule deepTools_heatmap:
    input: fls4deep
    output:
        matrix="deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_matrix.gz",
        profile="deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_plotProfile.png",
        clusters="deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_clusters.bed",
        hmap1="deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_heatmap1.png",
        hmap2="deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_heatmap2.png"
    threads:18
    log: "00log/{day}_{cond}-{rep}_HM_deepTools.log"
    message: "deepTools computeMatrix and plot heatmap for {wildcards.day}_{wildcards.cond}-{wildcards.rep}"
    params:
        bedfl = lambda wildcards: config["peak_files"][wildcards.day],
        name="plot profile for {day}_{cond}-{rep}",
        snames=samplenams,
        bedfl_nam=lambda wildcards:"Foxa2peaks" + "_" + wildcards.day,
        xlab=config["deeptools"]["xname"],
        refPoint=config["deeptools"]["refpoint"],
        suffix=config["deeptools"]["regionlabel"],
        colors=config["deeptools"]["colors"],
        nclust=config["deeptools"]["kmeans_ncluster"] 
    shell:
        """
        module load bioinfo-tools deepTools/3.1.0
        regionlabel=\"{params.suffix}{params.bedfl_nam}\"
        computeMatrix reference-point --referencePoint TSS -S {input} -R {params.bedfl} --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --skipZeros -o {output.matrix} --samplesLabel {params.snames} --numberOfProcessors {threads}
       # computeMatrix scale-regions -S {input} -R {params.bedfl} --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --regionBodyLength 10000 --skipZeros -o {output.matrix} --samplesLabel {params.snames} --numberOfProcessors {threads}
        plotProfile -m {output.matrix} -out {output.profile} --yAxisLabel "log2ratios(ChIP/Input)" --plotType="fill" --averageType="mean" --perGroup --colors {params.colors}  --refPointLabel \"{params.refPoint}\" --regionsLabel $regionlabel
       plotHeatmap -m {output.matrix} -out {output.hmap2}  --heatmapHeight 15 --heatmapWidth 5  --refPointLabel \"{params.refPoint}\" --kmeans {params.nclust}  --outFileSortedRegions {output.clusters} --xAxisLabel \"{params.xlab}\" --whatToShow "plot, heatmap and colorbar" --samplesLabel {params.snames}  
       plotHeatmap -m {output.matrix} -out {output.hmap1}  --heatmapHeight 15 --heatmapWidth 5 --refPointLabel \"{params.refPoint}\" --regionsLabel $regionlabel --xAxisLabel \"{params.xlab}\"  --whatToShow "plot, heatmap and colorbar" --samplesLabel {params.snames}
        """

rule clusters_annotate:
    input: "deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_clusters.bed"
    output: "deepTools_ashwini/profiles/{day}_{cond}-{rep}/{day}_{cond}-{rep}_Foxa2TFPeaks_{nam}_clusters_annotated.bed"
    threads:8
    params:
        rlocate = "/pica/h1/ashwini/anaconda2/lib/R/library",
        bedfl = lambda wildcards: config["peak_files"][wildcards.day]
    shell:
        """
        awk -F\\\\t '{{if(FNR==1){{fl=fl+1}}; if(fl==1){{k=$1":"$2"-"$3;ck[k]=$4}} else{{k=$4;if (k in ck){{print $0"\\t"ck[k]}}}}}}' {params.bedfl} {input} > {output}
        #Rscript ~/scripts/chipseq_pipeline/clusters_annotate.r {params.rlocate} {input} {output}
        """

rule gene_profiles:
    input: fls4deep
    output:
        matrix="deepTools_ashwini/profiles/{day}_{cond}-{rep}/markers/{day}_{cond}-{rep}_{genename}_TSS_matrix.gz",
        profile="deepTools_ashwini/profiles/{day}_{cond}-{rep}/markers/{day}_{cond}-{rep}_{genename}_TSS_plotProfile.png"
        
    params:
        snames=samplenams,
        bedfl=config["UCSC_bed"],
        genename=lambda wildcards: {wildcards.genename},
        colors=config["deeptools"]["colors"]
    message: "Plotting the profiles for gene {wildcards.genename}"
    threads:18
    shell:
        """
        module load bioinfo-tools deepTools/3.1.0
        computeMatrix scale-regions -S {input} -R temp.bed --beforeRegionStartLength 10000 --regionBodyLength 5000 --afterRegionStartLength 10000 --skipZeros -o {output.matrix} --samplesLabel {params.snames} --numberOfProcessors {threads}
        plotProfile -m {output.matrix} -out {output.profile} --yAxisLabel "mean Coverage" --yMin -30 --yMax 500 --plotType=fill --perGroup --colors {params.colors} --regionsLabel {params.genename}       
       """ 
   

rule enhancer_profiles:
    input: fls4deep
    output:
        matrix="deepTools_ashwini/profiles/{day}_{cond}-{rep}/enhancers/{day}_{cond}-{rep}_{genename}_FoxA2peaks_matrix.gz",
        profile="deepTools_ashwini/profiles/{day}_{cond}-{rep}/enhancers/{day}_{cond}-{rep}_{genename}_FoxA2peaks_plotProfile.png"
    params:
        snames=samplenams,
        bedfl= lambda wildcards: config["peak_files"][wildcards.day],
        genename=lambda wildcards:  {wildcards.genename},
        colors=config["deeptools"]["colors"][config["day"][0]][config["cond"][0]],
        refPoint=config["deeptools"]["refpoint"]
    message: "Plotting the profiles for gene {wildcards.genename}"
    threads:12
    shell:
        """
        module load bioinfo-tools deepTools/3.1.0
        grep -h {params.genename} {params.bedfl} | uniq > temp.bed
        region_name=$(grep -h {params.genename} {params.bedfl} | awk '{{print $4":"$1"-"$2":"$3}}' )
        computeMatrix reference-point --referencePoint TSS -S {input} -R temp.bed --beforeRegionStartLength 10000 --afterRegionStartLength 10000 --skipZeros -o {output.matrix} --samplesLabel {params.snames} --numberOfProcessors {threads}
        plotProfile -m {output.matrix} -out {output.profile} --yAxisLabel "mean Coverage" --plotType=fill --perGroup --colors {params.colors} --regionsLabel $region_name  --refPointLabel \"{params.refPoint}\" --yMin -10 --yMax 70
        """ 
