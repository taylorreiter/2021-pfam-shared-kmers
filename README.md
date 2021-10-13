## Measuring shared sequence similarity between PFAM clans

## Introduction

Biological sequence annotation underlies inference into the functional capacity of organisms.
Annotation can refer to many things, including prediction of open reading frame, gene name, ortholog family, or function. 
The strategies used to perform annotation vary depending on the final goal, but in general many strategies have been developed for each type of annotation. 
This document/repository focuses on ortholog annotation of protein sequences and is a precursor analysis supporting the development of an annotation approach that is robust for annotation fragmented and divergent sequences, particularly as are found in a metagenome assembly graph.

In this repository, we measure the jaccard similarity between protein sequences in approximately 19,000 PFAM families. 
The goal of this analysis was to determine whether PFAM ortholog annotation could be accurately acheived by protein k-mer overlap between a query and PFAM families.
In particular, we used jaccard similarity to inform whether a `spacegraphcats multifasta` approach might work for compact de Bruijn Graph (cDBG) and dominating set annotation.
Previously, we have used spacegraphcats `multifasta query` to peform cDBG annotation.
`multifasta query` transfers annotations from a nucleotide query sequence to any and all nodes in a cDBG that have at least one overlapping k-mer (*k* = 31). 
Given the specificity of long nucleotide k-mers, we have generally had success with this approach (low false positive, ok recall). 
However, the multifasta file that has to be used to annotate the graph needs to be highly specific, including sequences from the exact species that is present in the graph. 
We thought that working in amino acid space, as opposed to nucleotide space, might improve the recall of the `multifasta query` annotation approach. 
However, in its current implementation, `multifasta query` only requires a single k-mer to overlap with a node in the graph for the query annotation to transfer to that node.
Given this, there would have to be relatively little overlap between protein sequences in different ortholog families. 
Thus, this repository asseses the jaccard similarity between protein sequences in PFAM families. 

### Methods

We downloaded Pfam\_A from the PFAM database release 34.0. 
Pfam\_A contains 42,420,656 protein sequences classified as 19,179 PFAM families. 
The sequences are 90% non-redundant.

We generated a sourmash sketch (protein, k = 10, scaled = 1) for each PFAM family. 
Then we used sourmash compare to pairwise compare the jaccard similarity between each PFAM family.
We post-processed this matrix to remove comparisons to self.

### Results

The mean jaccard similarity among all comparisons was low. 

+ 593 PFAM families had at least 0.1 jaccard similarity with at least one other PFAM family.
+ 3,897 PFAM families had at least 0.01 jaccard similarity with at least one other PFAM family.
+ 13,816 PFAM families had at least 0.001 jaccard similarity with at least one other PFAM family.

The maximum jaccard similarity between two PFAM families was 0.921, however it was shared between PF08057 (Erythromycin resistance leader 2) and PF08051 (Erythromycin resistance leader 1).
Other PFAM families with high similarity display similar patterns, where they encode similar functions. 
  
### Conclusions

The fraction of shared k-mers is high enough that `multifasta query` approach may lend itself to high false positives, however many of these false positives may be functionally insignificant in a biological context.
A gather-like approach, where annotation is assigned based on maximum shared k-mers, may be more appropriate if it can be efficiently implemented to work on all cDBG nodes or all dominating sets.  

## Getting started

```
conda env create --name pfam --file environment.yml
conda activate pfam

snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=200000 --cluster "sbatch -t 10080 -J pfam -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```
