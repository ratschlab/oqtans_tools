#!/usr/bin/env python
"""
Extract genome annotation from a GFF3 (a tab delimited format
for storing sequence features and annotations:
http://www.sequenceontology.org/gff3.shtml) file.

Usage: ParseGFF.py in.gff3 out.mat 

Requirements: 
    Scipy :- http://scipy.org/ 
    Numpy :- http://numpy.org/ 

Copyright (C) 2010-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany 
"""
import re, sys, os
import scipy.io as sio
import numpy as np

def createExon(strand_p, five_p_utr, cds_cod, three_p_utr):
    """Create exon cordinates from UTR's and CDS region
    """
    exon_pos = []
    if strand_p == '+':        
        utr5_start, utr5_end = 0, 0
        if five_p_utr != []:
            utr5_start, utr5_end = five_p_utr[-1][0], five_p_utr[-1][1] 
        cds_5start, cds_5end = cds_cod[0][0], cds_cod[0][1]
        jun_exon = []
        if cds_5start-utr5_end == 0 or cds_5start-utr5_end == 1:
            jun_exon = [utr5_start, cds_5end]    
        if len(cds_cod) == 1:
            five_prime_flag = 0
            if jun_exon != []:
                five_p_utr = five_p_utr[:-1]
                five_prime_flag = 1
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
            jun_exon = []
            utr3_start, utr3_end = 0, 0
            if three_p_utr != []: 
                utr3_start = three_p_utr[0][0]
                utr3_end = three_p_utr[0][1]
            if utr3_start-cds_5end == 0 or utr3_start-cds_5end == 1:
                jun_exon = [cds_5start, utr3_end]
            three_prime_flag = 0
            if jun_exon != []: 
                cds_cod = cds_cod[:-1]
                three_p_utr = three_p_utr[1:]
                three_prime_flag = 1
            if five_prime_flag == 1 and three_prime_flag == 1:
                exon_pos.append([utr5_start, utr3_end])
            if five_prime_flag == 1 and three_prime_flag == 0:
                exon_pos.append([utr5_start, cds_5end])
                cds_cod = cds_cod[:-1]
            if five_prime_flag == 0 and three_prime_flag == 1:
                exon_pos.append([cds_5start, utr3_end])
            for cds in cds_cod:
                exon_pos.append(cds)
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
        else:    
            if jun_exon != []:
                five_p_utr = five_p_utr[:-1]
                cds_cod = cds_cod[1:]
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
            exon_pos.append(jun_exon) if jun_exon != [] else ''
            jun_exon = []
            utr3_start, utr3_end = 0, 0
            if three_p_utr != []:
                utr3_start = three_p_utr[0][0]
                utr3_end = three_p_utr[0][1]
            cds_3start = cds_cod[-1][0]
            cds_3end = cds_cod[-1][1]
            if utr3_start-cds_3end == 0 or utr3_start-cds_3end == 1:       
                jun_exon = [cds_3start, utr3_end]
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                three_p_utr = three_p_utr[1:]
            for cds in cds_cod:
                exon_pos.append(cds)
            exon_pos.append(jun_exon) if jun_exon != [] else ''
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
    elif strand_p == '-':
        utr3_start, utr3_end = 0, 0        
        if three_p_utr != []:
            utr3_start = three_p_utr[-1][0]
            utr3_end = three_p_utr[-1][1]
        cds_3start = cds_cod[0][0]
        cds_3end = cds_cod[0][1]
        jun_exon = []
        if cds_3start-utr3_end == 0 or cds_3start-utr3_end == 1:
            jun_exon = [utr3_start, cds_3end]  
        if len(cds_cod) == 1:    
            three_prime_flag = 0
            if jun_exon != []:
                three_p_utr = three_p_utr[:-1]
                three_prime_flag = 1
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
            jun_exon = []
            (utr5_start, utr5_end) = (0, 0)
            if five_p_utr != []:
                utr5_start = five_p_utr[0][0]
                utr5_end = five_p_utr[0][1]
            if utr5_start-cds_3end == 0 or utr5_start-cds_3end == 1:
                jun_exon = [cds_3start, utr5_end]
            five_prime_flag = 0
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                five_p_utr = five_p_utr[1:]
                five_prime_flag = 1
            if three_prime_flag == 1 and five_prime_flag == 1:
                exon_pos.append([utr3_start, utr5_end])
            if three_prime_flag == 1 and five_prime_flag == 0:
                exon_pos.append([utr3_start, cds_3end])
                cds_cod = cds_cod[:-1]
            if three_prime_flag == 0 and five_prime_flag == 1:
                exon_pos.append([cds_3start, utr5_end])        
            for cds in cds_cod:
                exon_pos.append(cds)
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
        else:
            if jun_exon != []:
                three_p_utr = three_p_utr[:-1]
                cds_cod = cds_cod[1:]
            for utr3 in three_p_utr:
                exon_pos.append(utr3)   
            if jun_exon != []:
                exon_pos.append(jun_exon)
            jun_exon = []
            (utr5_start, utr5_end) = (0, 0)
            if five_p_utr != []:
                utr5_start = five_p_utr[0][0]
                utr5_end = five_p_utr[0][1]    
            cds_5start = cds_cod[-1][0]
            cds_5end = cds_cod[-1][1]
            if utr5_start-cds_5end == 0 or utr5_start-cds_5end == 1:
                jun_exon = [cds_5start, utr5_end]
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                five_p_utr = five_p_utr[1:]
            for cds in cds_cod:
                exon_pos.append(cds)
            if jun_exon != []:
                exon_pos.append(jun_exon)    
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
    return exon_pos

def init_gene():
    """Initializing the gene structure
    """
    gene_details=dict(chr='', 
                    exons=[], 
                    gene_info={}, 
                    id='', 
                    is_alt_spliced=0, 
                    name='', 
                    source='', 
                    start='', 
                    stop='', 
                    strand='', 
                    transcripts=[])
    return gene_details

def FeatureValueFormat(singlegene):
    """Make feature value compactable to write in a .mat format
    """
    comp_exon = np.zeros((len(singlegene['exons']),), dtype=np.object)
    for i in range(len(singlegene['exons'])):
        comp_exon[i]= np.array(singlegene['exons'][i])
    singlegene['exons'] = comp_exon
    comp_transcripts = np.zeros((len(singlegene['transcripts']),), dtype=np.object)
    for i in range(len(singlegene['transcripts'])):
        comp_transcripts[i] = np.array(singlegene['transcripts'][i])
    singlegene['transcripts'] = comp_transcripts
    return singlegene 

def CreateGeneModels(genes_cmpt, transcripts_cmpt, exons_cmpt, utr3_cmpt, utr5_cmpt, cds_cmpt):
    """Creating Coding/Non-coding gene models from parsed GFF objects.
    """
    gene_counter, gene_models=1, []
    for gene_entry in genes_cmpt: ## Figure out the genes and transcripts associated feature 
        if gene_entry in transcripts_cmpt:
            gene=init_gene() 
            gene['id']=gene_counter
            gene['name']=gene_entry[1]
            gene['chr']=genes_cmpt[gene_entry]['chr']
            gene['source']=genes_cmpt[gene_entry]['source']
            gene['start']=genes_cmpt[gene_entry]['start']
            gene['stop']=genes_cmpt[gene_entry]['stop']
            gene['strand']=genes_cmpt[gene_entry]['strand']
            if not gene['strand'] in ['+', '-']:
                gene['strand']='.' # Strand info not known replaced with a dot symbol instead of None, ?, . etc.
            gene['gene_info']=dict(ID=gene_entry[1])
            if len(transcripts_cmpt[gene_entry])>1:
                gene['is_alt_spliced'] = 1
            for tids in transcripts_cmpt[gene_entry]: ## transcript section related tags 
                gene['transcripts'].append(tids['ID'])
                if len(exons_cmpt) != 0: 
                    if (gene['chr'], tids['ID']) in exons_cmpt:
                        exon_cod=[[feat_exon['start'], feat_exon['stop']] for feat_exon in exons_cmpt[(gene['chr'], tids['ID'])]]
                else: ## build exon coordinates from UTR3, UTR5 and CDS
                    utr5_pos, cds_pos, utr3_pos = [], [], []
                    if (gene['chr'], tids['ID']) in utr5_cmpt:
                        utr5_pos=[[feat_utr5['start'], feat_utr5['stop']] for feat_utr5 in utr5_cmpt[(gene['chr'], tids['ID'])]]
                    if (gene['chr'], tids['ID']) in cds_cmpt:
                        cds_pos=[[feat_cds['start'], feat_cds['stop']] for feat_cds in cds_cmpt[(gene['chr'], tids['ID'])]]
                    if (gene['chr'], tids['ID']) in utr3_cmpt:
                        utr3_pos=[[feat_utr3['start'], feat_utr3['stop']] for feat_utr3 in utr3_cmpt[(gene['chr'], tids['ID'])]]
                    exon_cod=createExon(gene['strand'], utr5_pos, cds_pos, utr3_pos) 
                if gene['strand']=='-':
                    if len(exon_cod) >1:
                        if exon_cod[0][0] > exon_cod[-1][0]:
                            exon_cod.reverse()
                if exon_cod: 
                    gene['exons'].append(exon_cod)
            gene=FeatureValueFormat(gene) # get prepare for MAT writing 
            gene_counter+=1
            gene_models.append(gene)
    return gene_models    

def GFFParse(gff_file):
    """Parsing GFF file based on feature relationship.
    """
    genes, utr5, exons=dict(), dict(), dict()
    transcripts, utr3, cds=dict(), dict(), dict()
    # TODO Include growing key words of different non-coding/coding transcripts 
    features=['mRNA', 'transcript', 'ncRNA', 'miRNA', 'pseudogenic_transcript', 'rRNA', 'snoRNA', 'snRNA', 'tRNA', 'scRNA']
    gff_handle=open(gff_file, "rU")
    for gff_line in gff_handle:
        gff_line=gff_line.strip('\n\r').split('\t')
        if re.match(r'#|>', gff_line[0]): # skip commented line and fasta identifier line 
            continue
        if len(gff_line)==1: # skip fasta sequence/empty line if present 
            continue 
        assert len(gff_line)==9, '\t'.join(gff_line) # not found 9 tab-delimited fields in this line     
        if '' in gff_line: # skip this line if there any field with an empty value
            print 'Skipping..', '\t'.join(gff_line)
            continue
        if gff_line[-1][-1]==';': # trim the last ';' character 
            gff_line[-1]=gff_line[-1].strip(';')
        if gff_line[2] in ['gene', 'pseudogene']:
            gid, gene_info=None, dict()
            gene_info['start']=int(gff_line[3])
            gene_info['stop']=int(gff_line[4])
            gene_info['chr']=gff_line[0]
            gene_info['source']=gff_line[1]
            gene_info['strand']=gff_line[6]
            for attb in gff_line[-1].split(';'):
                attb=attb.split('=') # gff attributes are separated by key=value pair 
                if attb[0]=='ID':
                    gid=attb[1]
                    break
            genes[(gff_line[0], gid)]=gene_info # store gene information based on the chromosome and gene symbol.
        elif gff_line[2] in features: 
            gid, mrna_info=None, dict() 
            mrna_info['start']=int(gff_line[3])
            mrna_info['stop']=int(gff_line[4])
            mrna_info['chr']=gff_line[0]
            mrna_info['strand']=gff_line[6]
            for attb in gff_line[-1].split(';'):
                attb=attb.split('=')
                if attb[0]=='Parent':
                    gid=attb[1]
                elif attb[0]=='ID':
                    mrna_info[attb[0]]=attb[1]
            for fid in gid.split(','): # child may be mapped to multiple parents ex: Parent=AT01,AT01-1-Protein 
                if (gff_line[0], fid) in transcripts:
                    transcripts[(gff_line[0], fid)].append(mrna_info)
                else:
                    transcripts[(gff_line[0], fid)]=[mrna_info]
        elif gff_line[2] in ['exon']:
            tids, exon_info=None, dict()
            exon_info['start']=int(gff_line[3])
            exon_info['stop']=int(gff_line[4])
            exon_info['chr']=gff_line[0]
            exon_info['strand']=gff_line[6]
            for attb in gff_line[-1].split(';'):
                attb=attb.split('=')
                if attb[0]=='Parent':
                    tids=attb[1]
                    break
            for tid in tids.split(','):
                if (gff_line[0], tid) in exons:
                    exons[(gff_line[0], tid)].append(exon_info)
                else:
                    exons[(gff_line[0], tid)]=[exon_info]
        elif gff_line[2] in ['five_prime_UTR']:
            utr5_info, tids=dict(), None
            utr5_info['start']=int(gff_line[3])
            utr5_info['stop']=int(gff_line[4])
            utr5_info['chr']=gff_line[0]
            utr5_info['strand']=gff_line[6]
            for attb in gff_line[-1].split(';'):
                attb=attb.split('=')
                if attb[0]=='Parent':
                    tids=attb[1]
                    break
            for tid in tids.split(','):
                if (gff_line[0], tid) in utr5:
                    utr5[(gff_line[0], tid)].append(utr5_info)
                else:
                    utr5[(gff_line[0], tid)]=[utr5_info]
        elif gff_line[2] in ['CDS']:
            cds_info, tids=dict(), None
            cds_info['start']=int(gff_line[3])
            cds_info['stop']=int(gff_line[4])
            cds_info['chr']=gff_line[0]
            cds_info['strand']=gff_line[6]
            for attb in gff_line[-1].split(';'):
                attb=attb.split('=')
                if attb[0]=='Parent':
                    tids=attb[1]
                    break
            for tid in tids.split(','):
                if (gff_line[0], tid) in cds:
                    cds[(gff_line[0], tid)].append(cds_info)
                else:
                    cds[(gff_line[0], tid)]=[cds_info]
        elif gff_line[2] in ['three_prime_UTR']:
            utr3_info, tids=dict(), None
            utr3_info['start']=int(gff_line[3])
            utr3_info['stop']=int(gff_line[4])
            utr3_info['chr']=gff_line[0]
            utr3_info['strand']=gff_line[6]
            for attb in gff_line[-1].split(';'):
                attb=attb.split('=')
                if attb[0]=='Parent':
                    tids=attb[1]
                    break
            for tid in tids.split(','):
                if (gff_line[0], tid) in utr3:
                    utr3[(gff_line[0], tid)].append(utr3_info)
                else:
                    utr3[(gff_line[0], tid)]=[utr3_info]
    gff_handle.close()
    return genes, transcripts, exons, utr3, utr5, cds

def __main__():
    """This function provides a best way to extract genome feature 
       information from a GFF3 file for the downstream processing.
    """
    try:
        gff_file = sys.argv[1]
        mat_file = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
    genes, transcripts, exons, utr3, utr5, cds=GFFParse(gff_file) 
    gene_models=CreateGeneModels(genes, transcripts, exons, utr3, utr5, cds)
    # TODO Write to matlab/octave struct instead of cell arrays.
    sio.savemat(mat_file, 
                    mdict=dict(genes=gene_models), 
                    format='5', 
                    oned_as='row')

if __name__=='__main__':
    __main__()
