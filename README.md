## Getting started
```
conda env create --name pfam --file environment.yml
conda activate pfam

snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=200000 --cluster "sbatch -t 10080 -J pfam -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```
