#! /usr/bin/env python

import os
import sys
import screed

fasta = str(snakemake.input)
log  = str("tmp_log.txt")
output_dirname = "outputs/pfam_fastas"

if not os.path.exists(output_dirname):
    try:
        os.makedirs(output_dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

with open(log, "w") as out_log:
    out_log.write(f"Splitting {fasta} by pfam. Writing files to {output_dirname} \n")
    pfam_records=[]
    pfam_name=""
    filenum = 0
    # looks like they're all in order, so we can loop through, write to file when encounter new pfam
    for n, record in enumerate(screed.open(fasta)):
        if n > 0 and n % 100000 == 0:
            out_log.write(f"working on {str(n)}th contig\n")
        this_pfam = record.name.rsplit(" ")[-1]
        # initiate pfam_name
        if not pfam_name:
            pfam_name = this_pfam
        if this_pfam != pfam_name:
            pfam_filename = pfam_name.split(";")[0]
            outfile = os.path.join(output_dirname, f"{pfam_filename}.fa")
            filenum+=1
            with open(outfile, "w") as out:
                for record in pfam_records:
                    out.write(f">{record.name}\n{record.sequence}\n")
            pfam_records=[]
            pfam_name = this_pfam
        # add record to pfam_records
        pfam_records+=[record]

    # catch last pfam
    pfam_filename = pfam_name.split(";")[0]
    outfile = os.path.join(output_dirname, f"{pfam_filename}.fa")
    filenum+=1
    with open(outfile, "w") as out:
        for record in pfam_records:
            out.write(f">{record.name}\n{record.sequence}\n")


    out_log.write(f"{str(n)} contigs written to {str(filenum)} individual pfam fasta files\n")

