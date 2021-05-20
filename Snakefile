# snakemake workflow for detecting viral sequences

import os

# preparing files
configfile: "config.yml"

INPUT_DIR = config["input"]
OUTPUT_DIR = config["output"]

IDS, = glob_wildcards(INPUT_DIR+"/{id}.fa")

os.system("chmod -R +x tools")

for sample in IDS:
    if not os.path.exists(OUTPUT_DIR+"/"+sample+"/logs"):
        os.makedirs(OUTPUT_DIR+"/"+sample+"/logs")

# rule that specifies the final expected output files
rule all:
    input:  
        expand(OUTPUT_DIR+"/{id}/final_predictions.fa", id=IDS),
        expand(OUTPUT_DIR+"/{id}/final_predictions_tax.tsv", id=IDS)

# virsorter 2
rule vs2_detect:
    input:  
        INPUT_DIR+"/{id}.fa"
    output: 
        OUTPUT_DIR+"/{id}/virsorter2/final-viral-combined.fa"
    params:
        outdir = OUTPUT_DIR+"/{id}/virsorter2",
        database = config['data']+"/vs2"
    conda:
        "envs/virsorter2.yml"
    shell:  
        "virsorter run -j 8 -i {input} --db-dir {params.database} -w {params.outdir} --min-length 10000 all || touch {output}"

rule vs2_rename:
    input:
        OUTPUT_DIR+"/{id}/virsorter2/final-viral-combined.fa"
    output:
        OUTPUT_DIR+"/{id}/virsorter2/viral_sequences.fa"
    conda:
        "envs/virsorter2.yml"
    shell:
        "tools/rename_multifasta_prefix.py -f {input} -p VS2 > {output}"

# deepvirfinder
rule dvf_detect:
    input:
        INPUT_DIR+"/{id}.fa"
    output:
        OUTPUT_DIR+"/{id}/deepvirfinder/{id}.fa_gt10000bp_dvfpred.txt"
    params:
        outdir = OUTPUT_DIR+"/{id}/deepvirfinder",
        database = config['data']+"/models"
    conda:
        "envs/deepvirfinder.yml"
    shell:
        "tools/DeepVirFinder/dvf.py -i {input} -m {params.database} -o {params.outdir} -l 10000 -c 4 || touch {output}"

rule dvf_filter:
    input:
        dvf = OUTPUT_DIR+"/{id}/deepvirfinder/{id}.fa_gt10000bp_dvfpred.txt",
        fa = INPUT_DIR+"/{id}.fa"
    output:
        OUTPUT_DIR+"/{id}/deepvirfinder/viral_sequences.fa"
    params:
        contigs = OUTPUT_DIR+"/{id}/deepvirfinder/viral_contigs.txt",
        fa_ori = OUTPUT_DIR+"/{id}/deepvirfinder/viral_contigs.fa",
    conda:
        "envs/deepvirfinder.yml"
    shell:
        """
        awk '{{if($3 > 0.9 && $4 < 0.01)print$1}}' {input.dvf} > {params.contigs}
        tools/select_seqs_by_IDs.py -i {input.fa} -d {params.contigs} -o {params.fa_ori}
        tools/rename_multifasta_prefix.py -f {params.fa_ori} -p DVF > {output}
        rm {params.fa_ori} {params.contigs}
        """

# vibrant
rule vibrant_detect:
    input:
        INPUT_DIR+"/{id}.fa"
    output:
        OUTPUT_DIR+"/{id}/vibrant/VIBRANT_{id}/VIBRANT_phages_{id}/{id}.phages_combined.fna"
    params:
        outdir = OUTPUT_DIR+"/{id}/vibrant",
        database = config['data']+"/databases"
    conda:
        "envs/vibrant.yml"
    shell:
        "VIBRANT_run.py -t 4 -d {params.database} -i {input} -folder {params.outdir} -l 10000 -no_plot || mkdir $(dirname {output}) && touch {output}"

rule vibrant_rename:
    input:
        OUTPUT_DIR+"/{id}/vibrant/VIBRANT_{id}/VIBRANT_phages_{id}/{id}.phages_combined.fna"
    output:
        OUTPUT_DIR+"/{id}/vibrant/viral_sequences.fa"
    conda:
        "envs/vibrant.yml"
    shell:
        "tools/rename_multifasta_prefix.py -f {input} -p VIBRANT > {output}"

# checkv
rule checkv_input:
    input:
        vs2 = OUTPUT_DIR+"/{id}/virsorter2/viral_sequences.fa",
        dvf = OUTPUT_DIR+"/{id}/deepvirfinder/viral_sequences.fa",
        vibrant = OUTPUT_DIR+"/{id}/vibrant/viral_sequences.fa"
    output:
        OUTPUT_DIR+"/{id}/checkv/viral_sequences.fa"
    shell:
        "cat {input.vs2} {input.dvf} {input.vibrant} > {output}"

rule checkv_analysis:
    input:
        OUTPUT_DIR+"/{id}/checkv/viral_sequences.fa"
    output:
        OUTPUT_DIR+"/{id}/checkv/checkv_output/quality_summary.tsv",
        OUTPUT_DIR+"/{id}/checkv/checkv_output/proviruses.fna",
        OUTPUT_DIR+"/{id}/checkv/checkv_output/viruses.fna" 
    params:
        outdir = OUTPUT_DIR+"/{id}/checkv/checkv_output/",
        database = config['data']+"/checkv/checkv-db-v1.0"
    conda:
        "envs/checkv.yml"
    shell:
        "checkv end_to_end -t 8 -d {params.database} {input} {params.outdir}"

rule checkv_select:
    input:
        OUTPUT_DIR+"/{id}/checkv/checkv_output/quality_summary.tsv"
    output:
        OUTPUT_DIR+"/{id}/checkv/checkv_output/quality_summary_filtered.tsv"
    shell:
        "awk -F '\t' '{{if($6 > $7 && $2 >= 10000 && $13 <= 1.0)print$0}}' {input} | grep -v 'contig >1.5x' > {output}"

rule checkv_filter:
    input:
        pred = OUTPUT_DIR+"/{id}/checkv/checkv_output/quality_summary_filtered.tsv",
        prov = OUTPUT_DIR+"/{id}/checkv/checkv_output/proviruses.fna",
        vir = OUTPUT_DIR+"/{id}/checkv/checkv_output/viruses.fna"
    output:
        all = OUTPUT_DIR+"/{id}/checkv/checkv_output/trimmed_predictions.fa",
        final = OUTPUT_DIR+"/{id}/checkv/checkv_output/filtered_predictions.fa"
    conda:
        "envs/checkv.yml"
    shell:
        """
        cat {input.prov} {input.vir} > {output.all}
        tools/filter_checkv.py {output.all} {input.pred} {output.final}
        """

# cluster sequences
rule dereplicate:
    input:
        OUTPUT_DIR+"/{id}/checkv/checkv_output/filtered_predictions.fa"
    output:
        pred = OUTPUT_DIR+"/{id}/final_predictions.fa",
        clst_dir = directory(OUTPUT_DIR+"/{id}/cdhit"),
        clst_file = OUTPUT_DIR+"/{id}/cdhit/filtered_predictions_nr.clstr"
    shell:
        """
        tools/cdhit/cd-hit-est -c 0.99 -i {input} -o {output.clst_dir}/filtered_predictions_nr -T 0 -M 0 -d 0
        mv {output.clst_dir}/filtered_predictions_nr {output.pred}
        """

# taxonomic classification
rule tax_class:
    input:
        OUTPUT_DIR+"/{id}/final_predictions.fa"
    output:
        tax = OUTPUT_DIR+"/{id}/final_predictions_tax.tsv",
        contig_ids = OUTPUT_DIR+"/{id}/demovir/trembl_ublast.viral.u.contigID.txt"
    params:
        demovir_dir = directory(OUTPUT_DIR+"/{id}/demovir"),
        database = config['data']+"/demovir"
    conda:
        "envs/checkv.yml"
    shell:
        """
        prodigal -a {params.demovir_dir}/proteins.faa -i {input} -p meta &> /dev/null
        tools/usearch -ublast {params.demovir_dir}/proteins.faa -db {params.database}/uniprot_trembl.viral.udb -evalue 1e-5 -trunclabels -blast6out {params.demovir_dir}/trembl_ublast.viral.txt -threads 4 &> /dev/null
        sort -u -k1,1 {params.demovir_dir}/trembl_ublast.viral.txt > {params.demovir_dir}/trembl_ublast.viral.u.txt
        cut -f 1,2 {params.demovir_dir}/trembl_ublast.viral.u.txt | sed 's/_[0-9]\+\t/\t/' | cut -f 1 | paste {params.demovir_dir}/trembl_ublast.viral.u.txt - > {output.contig_ids}
        tools/demovir.R {params.demovir_dir} {params.database}
        mv {params.demovir_dir}/DemoVir_assignments.txt {output.tax}
        """
