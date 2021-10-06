
def checkpoint_output_split_pfam_fasta_by_identifier(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.split_pfam_fasta_by_identifier.get(**wildcards).output[0]    
    file_names = expand("outputs/pfam_sigs/{pfam}.fa",
                        pfam = glob_wildcards(os.path.join(checkpoint_output, "{pfam}")).pfam)
    return file_names

rule all:
    input:
        checkpoint_output_split_pfam_fasta_by_identifier

rule download_pfam_A:
    output: "inputs/Pfam-A.fasta.gz"
    resources: mem_mb = 500
    threads: 1
    shell:'''
    wget -O {output} ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz
    '''

rule decompress_pfam_A:
    input: "inputs/Pfam-A.fasta.gz"
    output: "inputs/Pfam-A.fasta"
    resources: mem_mb = 500
    threads: 1
    shell:'''
    gunzip -c {input} > {output}
    '''

checkpoint split_pfam_fasta_by_identifier:
    input: "outputs/pfam/Pfam-A.fasta"
    output: directory("outputs/pfam_fastas")
    resources: mem_mb = 4000
    threads: 1
    conda: "envs/sourmash.yml"
    script: "scripts/split_fasta_by_pfam.py"

rule sourmash_signatures:
    input: "outputs/pfam_fastas/{pfam}.fa"
    output: "outputs/pfam_sigs/{pfam}_k10_scaled1.sig"
    conda: "envs/sourmash.yml"
    resources: mem_mb = 1000
    threads: 1
    shell:'''
    sourmash sketch protein -p k=10,scaled=1 -o {output} --name {wildcards.pfam} {input}
    '''

