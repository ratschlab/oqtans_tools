#!/usr/bin/env python 
"""
Program to prepare the genome of an organism in an internal format called GIO.
Usage:
python GenomePrep.py in.fasta out_dir 

Required:
	BioPython   
"""
import sys, os
from Bio import SeqIO

def __main__():
    try:
        fasta_in=sys.argv[1]
        out_dir=sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
    os.system("mkdir -p "+out_dir+"/genome")
    fasfh=open(fasta_in, 'rU')
    fasrec=[]
    for rec in SeqIO.parse(fasfh, "fasta"):
        fasrec.append(rec.id)
        out_genome=open(out_dir+"/genome/"+rec.id+'.flat', "w")
        gen_seq=rec.seq.tomutable()
        gen_seq=gen_seq.tostring().lower()
        out_genome.write(gen_seq)
        out_genome.close()
    fasfh.close()
    out_genome_config=open(out_dir+"/genome.config", "w")
    out_genome_config.write('BASEDIR '+out_dir+'\n\n')
    out_genome_config.write('CONTIGS '+str(len(fasrec))+'\n')
    for cid in fasrec:
        out_genome_config.write(cid+'\tgenome/'+cid+'.flat\tgenome/'+cid+'.dna\n')
    out_genome_config.write('\nALPHABET acgt\n\n')
    out_genome_config.write('ESTFILES 0\n\n')
    out_genome_config.write('CDNAFILES 0\n\n')
    out_genome_config.write('ANNOTATIONFILES 0')
    out_genome_config.close()

if __name__=="__main__":
    __main__()
