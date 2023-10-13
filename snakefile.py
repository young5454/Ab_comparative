
import os
import glob
configfile: "config.yaml"

# GROUP, STRAIN = config['groups']
# print(GROUP)
# print(STRAIN)
# wildcard_set = set()
# Define your dictionary of valid group-strain combinations
# INFO = {
#     "hvCRAB": ["A0006", "A0197", "A0220"],
#     "hvCRAB-MDR": ["A0062", "A0277", "A0220"]
# }

GROUP, STRAIN = glob_wildcards("/jupyterdem/assembly/{group}_{strain}/genome")
REF, = glob_wildcards("/jupyterdem/assembly/ref/genome/{ref}.fasta")

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
       
        ### Prokka I/O files
        # Reference genome & genbank
        expand("/jupyterdem/assembly/ref/genome/{ref}.fasta", ref=REF),
        expand("/jupyterdem/assembly/ref/{ref}.gb", ref=REF),

        # Prokka outputs - STRAINS
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
        
        # Prokka outputs - REFERENCE FASTA
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
        expand("/jupyterdem/pangenome/{group}_{strain}_ref_pairwise/summary_statistics.txt", zip, group=GROUP, strain=STRAIN)


# Rule to run Polypolish for polishing long-read assemblies with Illumina short reads
rule polypolish:
    input:
        genome_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}.fasta",
        reads1="/jupyterdem/assembly/{group}_{strain}/reads/{group}_{strain}_1.fastq.gz",
        reads2="/jupyterdem/assembly/{group}_{strain}/reads/{group}_{strain}_2.fastq.gz"
    output:
        aligned1 = "/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_aligned_1.sam",
        aligned2 = "/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_aligned_2.sam",
        filtered1 = "/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_filtered_1.sam",
        filtered2 = "/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_filtered_2.sam",
        polished_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta"
    params:
        threads=8, # Number of threads to use for multi-threaded processes
        path="/jupyterdem/assembly/{group}_{strain}/genome/"
    shell:
        """
        bwa index {input.genome_fasta}
        bwa mem -t {params.threads} -a {input.genome_fasta} {input.reads1} > {output.aligned1}
        bwa mem -t {params.threads} -a {input.genome_fasta} {input.reads2} > {output.aligned2}
        polypolish_insert_filter.py --in1 {output.aligned1} --in2 {output.aligned2} --out1 {output.filtered1} --out2 {output.filtered2}
        polypolish {input.genome_fasta} {output.filtered1} {output.filtered2} > {output.polished_fasta}
        """

### Need to add QUAST and BUSCO for QC

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
    shell:
        """
        prokka --outdir {params.out_dir} --prefix {params.prefix} --locustag {params.locustag} --cpus {params.threads} --kingdom {params.kingdom} --addgenes --force --proteins {input.ref_genbank} {input.ref_fasta}
        """


# Rule to run Prokka for annotating polished STRAIN FASTAs
rule prokka_strain:
    input:
        # Input FASTAs are polished FASTAs from Polypolish
        strain_fasta="/jupyterdem/assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta",
        strain_ref_genbank = expand("/jupyterdem/assembly/ref/{ref}.gb", ref=REF)
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
    shell:
        """
        prokka --outdir {params.out_dir} --prefix {params.prefix} --locustag {params.locustag} --cpus {params.threads} --kingdom {params.kingdom} --addgenes --force --proteins {input.strain_ref_genbank} {input.strain_fasta}
        """


# Rule to run Roary for STRAIN-to-REF 1:1 pairwise Pangenome analysis
rule roary_strain_ref_pairwise:
    input:
        # Input GFFs are annotated GFFs from Prokka
        ref_gff= expand("/jupyterdem/annotation/ref/{ref}/{ref}.gff", ref=REF),
        strain_gff="/jupyterdem/annotation/{group}_{strain}/{group}_{strain}.gff"
    output:
        # "./accessory.header.embl",
        # "./accessory.tab",
        # "./accessory_binary_genes.fa",
        # "./accessory_graph.dot",
        # "./blast_identity_frequency.Rtab",
        # "./clustered_proteins",
        # "./core_accessory.header.embl",
        # "./core_accessory.tab",
        # "./core_accessory_graph.dot",
        # "./core_alignment_header.embl",
        # "./core_gene_alignment.aln",
        # "./gene_presence_absence.Rtab",
        # "./gene_presence_absence.csv",
        # "./number_of_conserved_genes.Rtab",
        # "./number_of_genes_in_pan_genome.Rtab",
        # "./number_of_new_genes.Rtab",
        # "./number_of_unique_genes.Rtab",
        # "./pan_genome_reference.fa",
        # "./summary_statistics.txt"
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
    shell:
        """
        cd {params.out_dir}
        roary -e -p {params.threads} -v {input.ref_gff} {input.strain_gff}
        """




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