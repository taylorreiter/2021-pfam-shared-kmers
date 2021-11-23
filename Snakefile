minscaled=1
scaled= [1,10,100]
alphabet_info = {'protein': [7,10], 'dayhoff': ['16']}

p_str = []
alpha_ksizes=[]

for alpha, ksizes in alphabet_info.items():
    ak=expand(f"{alpha}-k{{k}}", k=ksizes)
    alpha_ksizes += ak
    # sourmash param strings
    kinf = expand("k={k}", k=ksizes)
    kstr = ",".join(kinf)
    this_pstr = f"-p {alpha},{kstr},scaled={minscaled},abund"
    p_str.append(this_pstr)

all_p_str = " ".join(p_str)
print(all_p_str)
print(alpha_ksizes)

def checkpoint_output_split_pfam_fasta_by_identifier(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.split_pfam_fasta_by_identifier.get(**wildcards).output[0]    
    file_names = expand("outputs/pfam_sigs/{pfam}.sig",
                        pfam = glob_wildcards(os.path.join(checkpoint_output, "{pfam}.fa")).pfam)
    return file_names


rule all:
    input:
        #"outputs/pfam_compare/pfam_compare.csv",
        expand("outputs/pfam_index/pfam_clans.{ak}-scaled{scaled}.sbt.zip", ak=alpha_ksizes, scaled=scaled),

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
    input: "inputs/Pfam-A.fasta"
    output: directory("outputs/pfam_fastas")
    resources: mem_mb = 4000
    threads: 1
    conda: "envs/sourmash.yml"
    script: "scripts/split_fasta_by_pfam.py"

rule sourmash_sketch:
    input: "outputs/pfam_fastas/{pfam}.fa"
    output: temp("outputs/pfam_sigs/{pfam}.sig")
    params:
        param_string= all_p_str
    conda: "envs/sourmash.yml"
    resources: mem_mb=lambda wildcards, attempt: attempt * 3000 
    threads: 1
    shell:'''
    sourmash sketch protein -p {params.param_string} -o {output} --name {wildcards.pfam} {input}
    '''

rule write_sketchlist:
    input: checkpoint_output_split_pfam_fasta_by_identifier
    output: temp("outputs/pfam_clans.siglist.txt")
    resources: mem_mb=lambda wildcards, attempt: attempt * 3000 
    threads: 1
    run:
        with open(str(output), 'w') as out:
            for inF in input:
                sketch_path = str(inF)
                if os.path.exists(sketch_path):
                    out.write(sketch_path + "\n")
                else:
                    print(f"Missing sketchfile! This should not happen {sketch_path}\n")


rule zip_sketches:
    input: "outputs/pfam_clans.siglist.txt"
    output: f"outputs/pfam_index/pfam_clans.{{alphabet}}-k{{ksize}}-scaled{minscaled}.zip"
    log: f"outputs/logs/sourmash/zip_sketches/pfam_clans.{{alphabet}}-k{{ksize}}-scaled{minscaled}.zip.log"
    conda: "envs/sourmash.yml"
    resources: mem_mb=lambda wildcards, attempt: attempt * 5000 
    threads: 1
    shell:
        '''
        sourmash sig cat --from-file {input} -o {output} 2> {log}
        '''

rule sbt_index:
    input: f"outputs/pfam_index/pfam_clans.{{alphabet}}-k{{ksize}}-scaled{minscaled}.zip"
    output: "outputs/pfam_index/pfam_clans.{alphabet}-k{ksize}-scaled{scaled}.sbt.zip"
    params:
        alpha_cmd= lambda w: f"--{w.alphabet}",
    log: "outputs/logs/sourmash/sbt-index/pfam_clans.{alphabet}-k{ksize}-scaled{scaled}.sbt.log"
    conda: "envs/sourmash.yml"
    resources: mem_mb=lambda wildcards, attempt: attempt * 50000 
    threads: 1
    shell:
        '''
        sourmash index {output} {input} {params.alpha_cmd} --ksize {wildcards.ksize} --scaled {wildcards.scaled} 2> {log}
        '''


rule sourmash_compare:
    #input: checkpoint_output_split_pfam_fasta_by_identifier
    input: f"outputs/pfam_index/pfam_clans.{{alphabet}}-k{{ksize}}-scaled{minscaled}.zip"
    output: 
        csv = f"outputs/pfam_compare/pfam_compare.{{alphabet}}-k{{ksize}}-scaled{minscaled}.csv",
        comp = f"outputs/pfam_compare/pfam_compare.{{alphabet}}-k{{ksize}}-scaled{minscaled}.comp"
    params:
        alpha_cmd= lambda w: f"--{w.alphabet}",
    log: f"outputs/logs/sourmash/pfam_compare/pfam_clans.{{alphabet}}-k{{ksize}}-scaled{minscaled}.compare.log"
    conda: "envs/sourmash.yml"
    resources:  mem_mb=lambda wildcards, attempt: attempt * 200000
    threads: 1
    shell:
        '''
        sourmash compare {input} {params.alpha_cmd} --ksize {wildcards.ksize} -o {output.comp} --csv {output.csv} 2> {log}
        '''
