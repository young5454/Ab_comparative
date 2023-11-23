import os
import glob
configfile: "groups_original.yml"


GROUP, STRAIN = glob_wildcards("/jupyterdem/assembly/{group}_{strain}/genome")
REF, = glob_wildcards("/jupyterdem/assembly/ref/genome/{ref}.fasta")

ADDGROUP, ADDSTRAIN = glob_wildcards("/jupyterdem/added_assembly/{group}_{strain}/genome")


# Rule to run Polypolish for all strains within a group
rule all:
    input:
        ### Polypolish I/O files
        # Long-read assembly
        expand("/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}.fasta", zip, group=GROUP, strain=STRAIN),
        
        # Illumina paired-end short reads
        expand("/jupyterdem/assembly/{group}_{strain}/reads/{group}_{strain}_1.fastq.gz", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/assembly/{group}_{strain}/reads/{group}_{strain}_2.fastq.gz", zip, group=GROUP, strain=STRAIN),

        # Raw alignments
        expand("/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_aligned_1.sam", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_aligned_2.sam", zip, group=GROUP, strain=STRAIN),

        # Filtered alignments
        expand("/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_filtered_1.sam", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_filtered_2.sam", zip, group=GROUP, strain=STRAIN),

        # Polished assembly
        expand("/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta", zip, group=GROUP, strain=STRAIN),
       
        # Added assembly
        expand("/jupyterdem/added_assembly/{group}_{strain}/genome/{group}_{strain}.fasta", zip, group=ADDGROUP, strain=ADDSTRAIN),
        
        ### Busco directory
        expand("/jupyterdem/busco/{group}_{strain}_busco/", zip, group=GROUP, strain=STRAIN),

        ### QUAST directory
        expand("/jupyterdem/quast/{group}_{strain}_quast/", zip, group=GROUP, strain=STRAIN),

        ### Prokka I/O files
        # Reference genome & genbank
        expand("/jupyterdem/assembly/ref/genome/{ref}.fasta", ref=REF),

        # Prokka outputs - REFERENCE 
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.err", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.ffn", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.fsa", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.gff", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.sqn", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.tsv", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.faa", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.fna", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.gbk", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.log", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.tbl", ref=REF),
        expand("/jupyterdem/annotation/ref/{ref}/{ref}.txt", ref=REF),
        
        # Prokka outputs - STRAINS 
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.err", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.ffn", zip, group=GROUP, strain=STRAIN),  
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.fsa", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gff", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.sqn", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.tsv", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.faa", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.fna", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gbk", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.log", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.tbl", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.txt", zip, group=GROUP, strain=STRAIN),
        
        # # Prokka outputs - STRAINS - ADDED
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.err", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.ffn", zip, group=ADDGROUP, strain=ADDSTRAIN),  
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.fsa", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.gff", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.sqn", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.tsv", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.faa", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.fna", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.gbk", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.log", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.tbl", zip, group=ADDGROUP, strain=ADDSTRAIN),
        # expand("/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.txt", zip, group=ADDGROUP, strain=ADDSTRAIN),
       
        ### Roary I/O files
        # Roary outputs - STRAIN-to-REF pairwise
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory.header.embl", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory.tab", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory_binary_genes.fa", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory_graph.dot", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/blast_identity_frequency.Rtab", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/clustered_proteins", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory.header.embl", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory.tab", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory_graph.dot", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_alignment_header.embl", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_gene_alignment.aln", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/gene_presence_absence.Rtab", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/gene_presence_absence.csv", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_conserved_genes.Rtab", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_genes_in_pan_genome.Rtab", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_new_genes.Rtab", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_unique_genes.Rtab", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/pan_genome_reference.fa", zip, group=GROUP, strain=STRAIN),
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/summary_statistics.txt", zip, group=GROUP, strain=STRAIN),

        # Roary output - Within-Group
        expand("/jupyterdem/pangenome2/{group}", group=GROUP),

        # Roary directory - Roary tmp for holding strain GFFs
        expand("/jupyterdem/roary_tmp/{group}/", group=GROUP),

        # Group directory
        expand("/jupyterdem/all_groups/{group}/gene_lists/", group=GROUP),
        expand("/jupyterdem/all_groups/{group}/faas/", group=GROUP),
        expand("/jupyterdem/all_groups/{group}/FASTA/", group=GROUP),

        # Roary plots directory#####
        expand("/jupyterdem/all_groups/{group}/roary_plots/", group=GROUP),

        ### EggNog-Mapper COG analysis output
        expand("/jupyterdem/all_groups/{group}/emapper-core/", group=GROUP),
        expand("/jupyterdem/all_groups/{group}/emapper-shells/", group=GROUP),
        expand("/jupyterdem/all_groups/{group}/cog_plots/", group=GROUP),


        # expand("/jupyterdem/all_groups/{group}/emapper-core/plot/", group=GROUP),
        # expand("/jupyterdem/all_groups/{group}/emapper-shells/plot/", group=GROUP),

        # Roary input - Within-GROUP Roary without reference
        # expand("/jupyterdem/roary_tmp/{group}/{group}_{strain}.gff", zip, group=GROUP, strain=STRAIN),

       
        # Roary output - Within-GROUP Roary without reference
        # expand("/jupyterdem/roary_tmp/{group}/accessory.header.embl", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/accessory.tab", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/accessory_binary_genes.fa", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/accessory_graph.dot", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/blast_identity_frequency.Rtab", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/clustered_proteins", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/core_accessory.header.embl", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/core_accessory.tab", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/core_accessory_graph.dot", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/core_alignment_header.embl", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/core_gene_alignment.aln", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/gene_presence_absence.Rtab", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/gene_presence_absence.csv", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/number_of_conserved_genes.Rtab", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/number_of_genes_in_pan_genome.Rtab", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/number_of_new_genes.Rtab", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/number_of_unique_genes.Rtab", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/pan_genome_reference.fa", group=GROUP),
        # expand("/jupyterdem/roary_tmp/{group}/summary_statistics.txt", group=GROUP)

        # expand("/jupyterdem/pangenome2/{group}_within_group/accessory.header.embl", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/accessory.tab", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/accessory_binary_genes.fa", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/accessory_graph.dot", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/blast_identity_frequency.Rtab", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/clustered_proteins", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/core_accessory.header.embl", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/core_accessory.tab", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/core_accessory_graph.dot", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/core_alignment_header.embl", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/core_gene_alignment.aln", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/gene_presence_absence.Rtab", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/gene_presence_absence.csv", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/number_of_conserved_genes.Rtab", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/number_of_genes_in_pan_genome.Rtab", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/number_of_new_genes.Rtab", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/number_of_unique_genes.Rtab", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/pan_genome_reference.fa", group=GROUP),
        # expand("/jupyterdem/pangenome2/{group}_within_group/summary_statistics.txt", group=GROUP)


# Rule to run Polypolish for polishing long-read assemblies with Illumina paired-end short reads
rule polypolish:
    input:
        genome_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}.fasta",
        reads1="/jupyterdem/assembly/{group}_{strain}/reads/{group}_{strain}_1.fastq.gz",
        reads2="/jupyterdem/assembly/{group}_{strain}/reads/{group}_{strain}_2.fastq.gz"
    output:
        aligned1="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_aligned_1.sam",
        aligned2="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_aligned_2.sam",
        filtered1="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_filtered_1.sam",
        filtered2="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_filtered_2.sam",
        polished_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta"
    params:
        threads=8, # Number of threads to use for multi-threaded processes
        path="/jupyterdem/assembly/{group}_{strain}/genome/"
    conda:
        "env_polypolish"
    shell:
        """
        bwa index {input.genome_fasta}
        bwa mem -t {params.threads} -a {input.genome_fasta} {input.reads1} > {output.aligned1}
        bwa mem -t {params.threads} -a {input.genome_fasta} {input.reads2} > {output.aligned2}
        polypolish_insert_filter.py --in1 {output.aligned1} --in2 {output.aligned2} --out1 {output.filtered1} --out2 {output.filtered2}
        polypolish {input.genome_fasta} {output.filtered1} {output.filtered2} > {output.polished_fasta}
        """


# Rule to run BUSCO assessments for polished assemblies and plot assessment plots
rule busco:
    input:
        # Input FASTAs are polished FASTAs from Polypolish
        strain_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta"
    output:
        out_dir=directory("/jupyterdem/busco/{group}_{strain}_busco/")
    params:
        lineage_path="/jupyterdem/busco_downloads/lineages/bacteria_odb10",
        output="{group}_{strain}_busco",
        out_path="/jupyterdem/busco/",
        sum_dir="/jupyterdem/busco/summaries/"
    conda:
        "env_busco"
    shell:
        """
        busco -m genome -i {input.strain_fasta} -o {params.output} --out_path {params.out_path} -l {params.lineage_path}
        mkdir -p {params.sum_dir}
        cd {output.out_dir}
        cp *.txt {params.sum_dir}
        
        generate_plot.py -wd {params.sum_dir}
        """


# Rule to run BUSCO assessments for polished assemblies and plot assessment plots
rule quast:
    input:
        # Input FASTAs are polished FASTAs from Polypolish
        strain_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta"
    output:
        out_dir=directory("/jupyterdem/quast/{group}_{strain}_quast/")
    params:
        threads=4
    conda:
        "env_quast"
    shell:
        """
        quast -o {output.out_dir} --threads {params.threads} {input.strain_fasta}
        """


# Rule to run Prokka for annotating REFERENCE FASTAs
rule prokka_ref:
    input:
        ref_fasta="/jupyterdem/assembly/ref/genome/{ref}.fasta",
        ref_genbank="/jupyterdem/assembly/ref/{ref}.gb"
    output:
        ref_err="/jupyterdem/annotation/ref/{ref}/{ref}.err",
        ref_ffn="/jupyterdem/annotation/ref/{ref}/{ref}.ffn",
        ref_fsa="/jupyterdem/annotation/ref/{ref}/{ref}.fsa",
        ref_gff="/jupyterdem/annotation/ref/{ref}/{ref}.gff",
        ref_sqn="/jupyterdem/annotation/ref/{ref}/{ref}.sqn",
        ref_tsv="/jupyterdem/annotation/ref/{ref}/{ref}.tsv",
        ref_faa="/jupyterdem/annotation/ref/{ref}/{ref}.faa",
        ref_fna="/jupyterdem/annotation/ref/{ref}/{ref}.fna",
        ref_gbk="/jupyterdem/annotation/ref/{ref}/{ref}.gbk",
        ref_log="/jupyterdem/annotation/ref/{ref}/{ref}.log",
        ref_tbl="/jupyterdem/annotation/ref/{ref}/{ref}.tbl",
        ref_txt="/jupyterdem/annotation/ref/{ref}/{ref}.txt"
    params:
        threads=8,
        out_dir="/jupyterdem/annotation/ref/{ref}",
        prefix="{ref}",
        locustag="{ref}",
        kingdom="Bacteria"
    conda:
        "env_prokka"
    shell:
        """
        prokka --outdir {params.out_dir} --prefix {params.prefix} --locustag {params.locustag} --cpus {params.threads} --kingdom {params.kingdom} --addgenes --force --proteins {input.ref_genbank} {input.ref_fasta}
        """


# Rule to run Prokka for annotating polished STRAIN FASTAs
rule prokka_strain:
    input:
        # Input FASTAs are polished FASTAs from Polypolish
        strain_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta",
        strain_ref_genbank=expand("/jupyterdem/assembly/ref/{ref}.gb", ref=REF)
    output:
        strain_err="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.err",
        strain_ffn="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.ffn",
        strain_fsa="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.fsa",
        strain_gff="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gff",
        strain_sqn="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.sqn",
        strain_tsv="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.tsv",
        strain_faa="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.faa",
        strain_fna="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.fna",
        strain_gbk="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gbk",
        strain_log="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.log",
        strain_tbl="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.tbl",
        strain_txt="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.txt"
    params:
        threads=8,
        out_dir="/jupyterdem/annotation/{group}_{strain}/",
        prefix="{group}_{strain}",
        locustag="{group}_{strain}",
        kingdom="Bacteria"
    conda:
        "env_prokka"
    shell:
        """
        prokka --outdir {params.out_dir} --prefix {params.prefix} --locustag {params.locustag} --cpus {params.threads} --kingdom {params.kingdom} --addgenes --force --proteins {input.strain_ref_genbank} {input.strain_fasta}
        """


# # Rule to run Prokka for annotating polished STRAIN FASTAs - For ADDED STRAINS
# rule prokka_strain_added:
#     input:
#         added_strain_fasta="/jupyterdem/added_assembly/{group}_{strain}/genome/{group}_{strain}.fasta",
#         strain_ref_genbank=expand("/jupyterdem/assembly/ref/{ref}.gb", ref=REF)
#     output:
#         added_strain_err="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.err",
#         added_strain_ffn="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.ffn",
#         added_strain_fsa="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.fsa",
#         added_strain_gff="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.gff",
#         added_strain_sqn="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.sqn",
#         added_strain_tsv="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.tsv",
#         added_strain_faa="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.faa",
#         added_strain_fna="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.fna",
#         added_strain_gbk="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.gbk",
#         added_strain_log="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.log",
#         added_strain_tbl="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.tbl",
#         added_strain_txt="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.txt"
#     params:
#         threads=8,
#         out_dir="/jupyterdem/added_annotation/{group}_{strain}/",
#         prefix="{group}_{strain}",
#         locustag="{group}_{strain}",
#         kingdom="Bacteria"
#     shell:
#         """
#         prokka --outdir {params.out_dir} --prefix {params.prefix} --locustag {params.locustag} --cpus {params.threads} --kingdom {params.kingdom} --addgenes --force --proteins {input.strain_ref_genbank} {input.added_strain_fasta}
#         """


# Rule to run Roary for STRAIN-to-REF 1:1 pairwise Pangenome analysis
rule roary_strain_ref_pairwise:
    input:
        # Input GFFs are annotated GFFs from Prokka
        ref_gff=expand("/jupyterdem/annotation/ref/{ref}/{ref}.gff", ref=REF),
        strain_gff="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gff"
    output:
        accessory_header="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory.header.embl",
        accessory_tab="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory.tab",
        accessory_binary_genes="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory_binary_genes.fa",
        accessory_graph="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory_graph.dot",
        blast_identity_frequency="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/blast_identity_frequency.Rtab",
        clustered_proteins="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/clustered_proteins",
        core_accessory_header="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory.header.embl",
        core_accessory="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory.tab",
        core_accessory_graph="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory_graph.dot",
        core_alignment_header="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_alignment_header.embl",
        core_gene_alignment="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_gene_alignment.aln",
        gene_presence_absence_rtab="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/gene_presence_absence.Rtab",
        gene_presence_absence_csv="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/gene_presence_absence.csv",
        number_of_conserved_genes="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_conserved_genes.Rtab",
        number_of_genes_in_pan_genome="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_genes_in_pan_genome.Rtab",
        number_of_new_genes="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_new_genes.Rtab",
        number_of_unique_genes="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_unique_genes.Rtab",
        pan_genome_reference="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/pan_genome_reference.fa",
        summary_statistics="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/summary_statistics.txt"
    params:
        threads=8,
        out_dir="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/",
        pident=90
    conda:
        "env_roary"
    shell:
        """
        cd {params.out_dir}
        roary -e -p {params.threads} -i {params.pident} -v {input.ref_gff} {input.strain_gff}
        """


# Rule to make roary_tmp group folder
rule move_gff_files:
    ###
    input:
        strain_gff=expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gff", zip, group=GROUP, strain=STRAIN)
    output:
        tmp_dir=directory("/jupyterdem/roary_tmp/{group}/")
    params:
        python="python3.6",
        script="move_gffs.sh",
        group_info="groups_original.yml",
        tmp_dir=directory("/jupyterdem/roary_tmp/"),
    conda:
        "env_roary"
    shell:
        """
        python shellmake.py --group_yml {params.group_info} --save_path {params.tmp_dir} --script {params.script}
        bash {params.script}
        """


# Rule to run Roary for STRAINS within GROUPS
rule roary_within_group:
    input:
        tmp_dir="/jupyterdem/roary_tmp/{group}/"
    output:
        out_dir=directory("/jupyterdem/pangenome2/{group}/")
    params:
        threads=8,
        pident=90,
    conda:
        "env_roary"
    shell:
        """
        mkdir {output.out_dir}
        cd {input.tmp_dir}
        roary -e -p {params.threads} -i {params.pident} -v *.gff
        mv *.gff {output.out_dir}
        """


#############################################################################################################################
# Rule to make six gene lists from ROARY gene_presence_absence.csv
# core_all, core_hypo, core_nonhypo, shells_all, shells_hypo, shells_nonhypo 
rule gene_list_maker:
    input:
        dummy=rules.roary_within_group.output.out_dir
    output:
        out_dir=directory("/jupyterdem/all_groups/{group}/gene_lists/")
    params:
        tmp_dir=directory("/jupyterdem/roary_tmp/{group}/")
    conda:
        "env_emapper"
    shell:
        """
        mkdir -p {output.out_dir}
        python gene_list_maker.py --csv {params.tmp_dir}/gene_presence_absence.csv --save_path {output.out_dir}/
        """


# Rule to copy faa files from PROKKA output to GROUP directory
rule move_faa_files:
    input:
        strain_faa=expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.faa", zip, group=GROUP, strain=STRAIN)
    output:
        faa_dir=directory("/jupyterdem/all_groups/{group}/faas/")
    params:
        script="move_faa.sh",
        group_info="groups_original.yml",
        group_dir="/jupyterdem/all_groups/"
    conda:
        "env_emapper"
    shell:
        """
        python shellmake2.py --group_yml {params.group_info} --save_path {params.group_dir} --script {params.script}
        bash {params.script}
        """


# Rule to curate six FASTA files based gene lists created
rule fasta_curation:
    input:
        faa_dir=rules.move_faa_files.output.faa_dir,
        text_path=rules.gene_list_maker.output.out_dir
    output:
        fasta_dir=directory("/jupyterdem/all_groups/{group}/FASTA/")
    params:
        tmp_dir=directory("/jupyterdem/roary_tmp/{group}/"),
        statistics=directory("/jupyterdem/all_groups/{group}/statistics.txt")
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.fasta_dir}
        python core_shell_statistics.py --text_path {input.text_path}/ --save_path {output.fasta_dir}/ --faa_path {input.faa_dir}/ \
        --gpa {params.tmp_dir}/gene_presence_absence.csv --summary {params.tmp_dir}/summary_statistics.txt > {params.statistics}
        """


### 11/20/23 Change
# Rule to conduct COG analysis by running Egg-NOG-mapper for each group's core nonhypo FASTAs
rule cog_analysis_core:
    input:
        # Inputs are curated FASTAs from rule fasta_curation
        fasta_dir=rules.fasta_curation.output.fasta_dir
    output:
        emapper_dir=directory("/jupyterdem/all_groups/{group}/emapper-core/")
    params:
        # Parameters for Egg-NOG-mapper
        pident=30,
        evalue=0.001,
        score=60,
        query_cover=70,
        subject_cover=70,
        output="{group}-core",
        cpu=8
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.emapper_dir}
        emapper.py -i {input.fasta_dir}/core_nonhypo.fasta --pident {params.pident} --evalue {params.evalue} \
        --score {params.score} --query_cover {params.query_cover} --subject_cover {params.subject_cover} \
        --output {params.output} --output_dir {output.emapper_dir} --cpu {params.cpu}
        """


# Rule to conduct COG analysis by running Egg-NOG-mapper for each group's shell nonhypo FASTAs
rule cog_analysis_shells:
    input:
        # Inputs are curated FASTAs from rule fasta_curation
        fasta_dir=rules.fasta_curation.output.fasta_dir
    output:
        emapper_dir=directory("/jupyterdem/all_groups/{group}/emapper-shells/")
    params:
        # Parameters for Egg-NOG-mapper
        pident=30,
        evalue=0.001,
        score=60,
        query_cover=70,
        subject_cover=70,
        output="{group}-shells",
        cpu=8
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.emapper_dir}
        emapper.py -i {input.fasta_dir}/shells_nonhypo.fasta --pident {params.pident} --evalue {params.evalue} \
        --score {params.score} --query_cover {params.query_cover} --subject_cover {params.subject_cover} \
        --output {params.output} --output_dir {output.emapper_dir} --cpu {params.cpu}
        """


# Rule to draw nested pie chart for each group's core/shell COG analysis results
rule cog_visualization:
    input:
        # Inputs are annotated tsv files generated from Egg-NOG-mapper
        core_emapper_dir=rules.cog_analysis_core.output.emapper_dir,
        shells_emapper_dir=rules.cog_analysis_shells.output.emapper_dir
    output:
        plot_dir=directory("/jupyterdem/all_groups/{group}/cog_plots/")
    params:
        core_hypo_path="/jupyterdem/all_groups/{group}/FASTA/core_hypo.fasta",
        shells_hypo_path="/jupyterdem/all_groups/{group}/FASTA/shells_hypo.fasta",

        group_name="{group}",
        types_core="Core",
        types_shells="Shells",
        # save_path="/jupyterdem/all_groups/{group}/emapper-core/",
        core_statistics="/jupyterdem/all_groups/{group}/emapper-core/statistics.txt",
        shells_statistics="/jupyterdem/all_groups/{group}/emapper-shells/statistics.txt"
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.plot_dir}

        python cog_analysis_nested.py --tsv_file {input.core_emapper_dir}/*.annotations \
        --hypo_path {params.core_hypo_path} --group_name {params.group_name} \
        --types {params.types_core} --save_path {output.plot_dir}/ > {params.core_statistics}

        python cog_analysis_nested.py --tsv_file {input.shells_emapper_dir}/*.annotations \
        --hypo_path {params.shells_hypo_path} --group_name {params.group_name} \
        --types {params.types_shells} --save_path {output.plot_dir}/ > {params.shells_statistics}
        """


# Rule to draw ROARY plots - slightly modified version of roary_plots.py by Marco Galardini
# https://github.com/sanger-pathogens/Roary/blob/master/contrib/roary_plots/README.md
rule roary_visualization:
    input:
        dummy=rules.roary_within_group.output.out_dir
    output:
        plot_dir=directory("/jupyterdem/all_groups/{group}/roary_plots/")
    params:
        tmp_dir=directory("/jupyterdem/roary_tmp/{group}/"),
        group="{group}"
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.plot_dir}
        FastTree -nt -gtr {params.tmp_dir}/core_gene_alignment.aln > {params.group}_tree.newick

        python roary_plots.py --labels {params.group}_tree.newick {params.tmp_dir}gene_presence_absence.csv \
        --save_path {output.plot_dir}

        mv {params.group}_tree.newick {output.plot_dir}
        """







# # Rule to draw nested pie chart for each group's shells COG analysis results
# rule cog_visualization_shells:
#     input:
#         # Inputs are annotated tsv files generated from Egg-NOG-mapper
#         emapper_dir=rules.cog_analysis_shells.output.emapper_dir
#     params:
#         hypo_path="/jupyterdem/all_groups/{group}/FASTA/shells_hypo.fasta",
#         group_name="{group}",
#         types="Shells",
#         # save_path="/jupyterdem/all_groups/{group}/emapper-shells/",
#         statistics="/jupyterdem/all_groups/{group}/emapper-shells/statistics.txt"
#     conda:
#         "env_emapper"
#     shell:
#         """
#         python cog_analysis_nested.py --tsv_file {input.emapper_dir}/*.annotations \
#         --hypo_path {params.hypo_path} --group_name {params.group_name} \
#         --types {params.types} --save_path {output.plot_dir} > {params.statistics}
#         """



        









# rule gene_list_maker:
#     input:
#         gpa=rules.roary_within_group.output.gpa
#     output:
#         out_dir=directory("/jupyterdem/all_groups/{group}/gene_lists/")
#     conda:
#         "env_emapper"
#     shell:
#         """
#         mkdir -p {output.out_dir}
#         python gene_list_maker.py --csv {input.gpa} --save_path {output.out_dir}
#         """



# # Rule to run Roary for STRAIN-to-REF 1:1 pairwise Pangenome analysis - For ADDED STRAINS
# rule roary_strain_ref_pairwise_added:
#     input:
#         # Input GFFs are annotated GFFs from Prokka
#         ref_gff=expand("/jupyterdem/annotation/ref/{ref}/{ref}.gff", ref=REF),
#         added_strain_gff="/jupyterdem/added_annotation/{group}_{strain}/{group}_{strain}.gff"
#     output:
#         accessory_header="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory.header.embl",
#         accessory_tab="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory.tab",
#         accessory_binary_genes="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory_binary_genes.fa",
#         accessory_graph="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/accessory_graph.dot",
#         blast_identity_frequency="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/blast_identity_frequency.Rtab",
#         clustered_proteins="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/clustered_proteins",
#         core_accessory_header="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory.header.embl",
#         core_accessory="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory.tab",
#         core_accessory_graph="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_accessory_graph.dot",
#         core_alignment_header="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_alignment_header.embl",
#         core_gene_alignment="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/core_gene_alignment.aln",
#         gene_presence_absence_rtab="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/gene_presence_absence.Rtab",
#         gene_presence_absence_csv="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/gene_presence_absence.csv",
#         number_of_conserved_genes="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_conserved_genes.Rtab",
#         number_of_genes_in_pan_genome="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_genes_in_pan_genome.Rtab",
#         number_of_new_genes="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_new_genes.Rtab",
#         number_of_unique_genes="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/number_of_unique_genes.Rtab",
#         pan_genome_reference="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/pan_genome_reference.fa",
#         summary_statistics="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/summary_statistics.txt"
#     params:
#         threads=8,
#         out_dir="/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/",
#         pident=90
#     shell:
#         """
#         cd {params.out_dir}
#         roary -e -p {params.threads} -i {params.pident} -v {input.ref_gff} {input.strain_gff}
#         """

# ## Rule to make roary_tmp group folder

# # Rule to run Roary for STRAINS within GROUPS
# rule roary_within_group_added:
#     output:
#         out_dir=directory("/jupyterdem/pangenome2/{group}")
#     params:
#         threads=8,
#         tmp_dir=directory("/jupyterdem/roary_tmp/{group}/"),
#         pident=90
#     shell:
#         """
#         mkdir {output.out_dir}
#         cd {params.tmp_dir}
#         roary -e -p {params.threads} -i {params.pident} -v *.gff
#         mv *.gff {output.out_dir}
#         """


# # Rule to copy and save STRAIN GFFs into a tmp folder for within-group Roary
# rule copy_strain_gffs:
#     input:
#         strain_gff=expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gff", zip, group=GROUP, strain=STRAIN)
#     output:
#         out_dir=directory("/jupyterdem/roary_tmp/{group}/")
#     run:
#         input_string_gff_list = str(input.strain_gff).split(' ')
#         print(input_string_gff_list)
#         # Create GROUP directories
#         shell("mkdir -p {output.out_dir}")
#         # Check if the 'GROUP' part of the 'strain_gff' path matches the {group} in 'out_dir'
#         for strain in input_string_gff_list:
#             group_id = strain.split("/")[-1].split('_')[0]
#             if group_id == output.out_dir.split("/")[-1]:
#                 # Part 2: Execute if the if statement is true
#                 shell("cp {strain} {output.out_dir}")


# def list_files_in_directory(directory_path):
#     try:
#         # Get a list of all files in the specified directory
#         files = os.listdir(directory_path)
        
#         # Create a string with each file on a new line
#         file_list_string = " ".join(files)

#         return file_list_string

#     except Exception as e:
#         return str(e)

# # Example usage:
# directory_path = "/jupyterdem/roary_tmp/hvCRAB/"
# file_list = list_files_in_directory(directory_path)

# print(file_list)

# rule roary_within_group:
#     input:
#         expand("/jupyterdem/roary_tmp/{group}/{group}_{strain}.gff", zip, group=GROUP, strain=STRAIN)
#     output:
#         "roary.done"
#     params:
#         threads=8,
#         pident=90,
#         work_dir="/jupyterdem/roary_tmp/{group}/",
#         out_dir="/jupyterdem/pangenome2/{group}/"
#     shell:
#         """
#         cd {params.work_dir}
#         roary -e -p {params.threads} -f {params.out_dir} -i {params.pident} -v {input} && touch {output}
#         """


# rule copy_strain_gffs:
#     input:
#         strain_gff=expand("/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gff", zip, group=GROUP, strain=STRAIN)
#     output:
#         out_dir=directory("/jupyterdem/roary_tmp/{group}/")
#     shell:
#         """
#         mkdir {output.out_dir}
#         cp {input.strain_gff} {output.out_dir}
#         """



    # run:
    #     shell("cd {params.tmp_dir}")
    #     shell("roary -e -p {params.threads} -i {params.pident} -v *.gff")
    #     shell("mv *.gff {output.out_dir}")
    # output:
    #     accessory_header="/jupyterdem/roary_tmp/{group}/accessory.header.embl",
    #     accessory_tab="/jupyterdem/roary_tmp/{group}/accessory.tab",
    #     accessory_binary_genes="/jupyterdem/roary_tmp/{group}/accessory_binary_genes.fa",
    #     accessory_graph="/jupyterdem/roary_tmp/{group}/accessory_graph.dot",
    #     blast_identity_frequency="/jupyterdem/roary_tmp/{group}/blast_identity_frequency.Rtab",
    #     clustered_proteins="/jupyterdem/roary_tmp/{group}/clustered_proteins",
    #     core_accessory_header="/jupyterdem/roary_tmp/{group}/core_accessory.header.embl",
    #     core_accessory="/jupyterdem/roary_tmp/{group}/core_accessory.tab",
    #     core_accessory_graph="/jupyterdem/roary_tmp/{group}/core_accessory_graph.dot",
    #     core_alignment_header="/jupyterdem/roary_tmp/{group}/core_alignment_header.embl",
    #     core_gene_alignment="/jupyterdem/roary_tmp/{group}/core_gene_alignment.aln",
    #     gene_presence_absence_rtab="/jupyterdem/roary_tmp/{group}/gene_presence_absence.Rtab",
    #     gene_presence_absence_csv="/jupyterdem/roary_tmp/{group}/gene_presence_absence.csv",
    #     number_of_conserved_genes="/jupyterdem/roary_tmp/{group}/number_of_conserved_genes.Rtab",
    #     number_of_genes_in_pan_genome="/jupyterdem/roary_tmp/{group}/number_of_genes_in_pan_genome.Rtab",
    #     number_of_new_genes="/jupyterdem/pangenome2/roary_tmp/{group}/number_of_new_genes.Rtab",
    #     number_of_unique_genes="/jupyterdem/pangenome2/roary_tmp/{group}/number_of_unique_genes.Rtab",
    #     pan_genome_reference="/jupyterdem/pangenome2/roary_tmp/{group}/pan_genome_reference.fa",
    #     summary_statistics="/jupyterdem/pangenome2/roary_tmp/{group}/summary_statistics.txt"
    # params:
    #     threads=8,
    #     out_dir="/jupyterdem/roary_tmp/{group}/",
    #     pident=90
    # shell:
    #     """
    #     cd {params.out_dir}
    #     roary -e -p {params.threads} -i {params.pident} -v *.gff
    #     """


# rule roary_within_group:
#     input:
#         # Input GFFs are annotated GFFs from Prokka
#         strain_gff=expand("/jupyterdem/roary_tmp/{group}/{group}_{strain}.gff", zip, group=GROUP, strain=STRAIN)
#     output:
#         accessory_header="/jupyterdem/pangenome2/{group}_within_group/accessory.header.embl",
#         accessory_tab="/jupyterdem/pangenome2/{group}_within_group/accessory.tab",
#         accessory_binary_genes="/jupyterdem/pangenome2/{group}_within_group/accessory_binary_genes.fa",
#         accessory_graph="/jupyterdem/pangenome2/{group}_within_group/accessory_graph.dot",
#         blast_identity_frequency="/jupyterdem/pangenome2/{group}_within_group/blast_identity_frequency.Rtab",
#         clustered_proteins="/jupyterdem/pangenome2/{group}_within_group/clustered_proteins",
#         core_accessory_header="/jupyterdem/pangenome2/{group}_within_group/core_accessory.header.embl",
#         core_accessory="/jupyterdem/pangenome2/{group}_within_group/core_accessory.tab",
#         core_accessory_graph="/jupyterdem/pangenome2/{group}_within_group/core_accessory_graph.dot",
#         core_alignment_header="/jupyterdem/pangenome2/{group}_within_group/core_alignment_header.embl",
#         core_gene_alignment="/jupyterdem/pangenome2/{group}_within_group/core_gene_alignment.aln",
#         gene_presence_absence_rtab="/jupyterdem/pangenome2/{group}_within_group/gene_presence_absence.Rtab",
#         gene_presence_absence_csv="/jupyterdem/pangenome2/{group}_within_group/gene_presence_absence.csv",
#         number_of_conserved_genes="/jupyterdem/pangenome2/{group}_within_group/number_of_conserved_genes.Rtab",
#         number_of_genes_in_pan_genome="/jupyterdem/pangenome2/{group}_within_group/number_of_genes_in_pan_genome.Rtab",
#         number_of_new_genes="/jupyterdem/pangenome2/{group}_within_group/number_of_new_genes.Rtab",
#         number_of_unique_genes="/jupyterdem/pangenome2/{group}_within_group/number_of_unique_genes.Rtab",
#         pan_genome_reference="/jupyterdem/pangenome2/{group}_within_group/pan_genome_reference.fa",
#         summary_statistics="/jupyterdem/pangenome2/{group}_within_group/summary_statistics.txt"
#     params:
#         threads=8,
#         work_dir="/jupyterdem/roary_tmp/{group}/",
#         out_dir="/jupyterdem/pangenome2/{group}_within_group/",
#         pident=90
#     shell:
#         """
#         cd {params.work_dir}
#         roary -f {params.out_dir} -e -p {params.threads} -i {params.pident} -v *.gff
#         """

# # Rule to run ABRicate for annotating virulence factors within REFERENCE FASTAs
# rule abricate_ref:
#     input:
#         # Input FASTAs are polished FASTAs using Polypolish
#         ref_fasta="/jupyterdem/assembly/ref/genome/{ref}.fasta"
#     output:
#         ref_abricate_tab="/jupyterdem/annotation/ref/{ref}/{ref}_abricate.tsv"
#     params:
#         threads=8,
#         db="vfdb"
#     shell:
#         """
#         abricate-get_db --db {params.db} --force
#         abricate --db {params.db} {input.ref_fasta}
#         """


# # Rule to run ABRicate for annotating virulence factors within STRAIN FASTAs
# rule abricate_strain:
#     input:
#         # Input FASTAs are polished FASTAs using Polypolish
#         strain_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta"
#     output:
#         strain_abricate_tsv="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}_abricate.tsv"
#     params:
#         threads=8,
#         db="vfdb"
#     shell:
#         """
#         abricate-get_db --db {params.db} --force
#         abricate --db {params.db} {input.strain_fasta}
#         """