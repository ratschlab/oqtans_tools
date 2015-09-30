#!/usr/bin/env python
"""
Common utility functions
"""

import os 
import re
import sys 
import gzip 
import bz2
import numpy 

def init_gene_DE():
    """
    Initializing the gene structure for DE
    """
    gene_det = [('id', 'f8'), 
                    ('chr', 'S15'), 
                    ('exons', numpy.dtype),
                    ('gene_info', numpy.dtype),
                    ('is_alt_spliced', 'f8'), 
                    ('name', 'S25'),
                    ('source', 'S25'),
                    ('start', 'f8'),
                    ('stop', 'f8'), 
                    ('strand', 'S2'), 
                    ('transcripts', numpy.dtype)]

    return gene_det

def _open_file(fname):
    """
    Open the file (supports .gz .bz2) and returns the handler
    """
    try:
        if os.path.splitext(fname)[1] == ".gz":
            FH = gzip.open(fname, 'rb')
        elif os.path.splitext(fname)[1] == ".bz2":
            FH = bz2.BZ2File(fname, 'rb')
        else:
            FH = open(fname, 'rU')
    except Exception as error:
        sys.exit(error)
    return FH

def make_Exon_cod(strand_p, five_p_utr, cds_cod, three_p_utr):
    """
    Create exon cordinates from UTR's and CDS region
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
