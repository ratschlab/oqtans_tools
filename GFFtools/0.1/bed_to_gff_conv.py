#!/usr/bin/env python
"""
Convert data in 12 column Browser Extensible Data (BED) format file 
to Generic Feature Format Version 3 (GFF3).

Usage: python bed_to_gff_conv.py in.bed > out.gff 

Copyright (C) 
    2009-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany.
    2012-2013 Memorial Sloan-Kettering Cancer Center New York City, USA.
"""
import re, sys
from common_util import _open_file

def bed_parse(qfile, source_name):
    """
    Process BED file 
    """

    BEDfh = _open_file(qfile)
    print '##gff-version 3'
    
    for rec in BEDfh:
        rec = rec.strip( '\n\r' )

        if not rec or rec[0] in  ['#']:
            continue
        if not re.search('\t', rec):
            continue

        line = rec.split('\t')
        assert len(line) >= 12, rec
        # checking the consistency b/w start of exon and number of exons
        if len(line[-1].split(',')) != len(line[-2].split(',')):
            continue 

        rstart = line[-1].split(',')
        if rstart[-1] == '': 
            rstart.pop()
        exon_len = line[-2].split(',')
        if exon_len[-1] == '': 
            exon_len.pop()

        if line[5] != '+' and line[5] != '-':
            line[5] = '.' #replace the unknown strand with '.' 

        pline = [str(line[0]),
                    source_name,
                    'gene',
                    str(int(line[1]) + 1),
                    line[2],
                    line[4],
                    line[5],
                    '.',
                    'ID=Gene:' + line[3] + ';Name=Gene:' + line[3]]
        print '\t'.join(pline) 

        pline = [str(line[0]),
                    source_name,
                    'transcript',
                    str(int(line[1]) + 1),
                    line[2],
                    line[4],
                    line[5],
                    '.',
                    'ID=' + line[3] + ';Name=' + line[3] + ';Parent=Gene:' + line[3]]
        print '\t'.join(pline)

        st = int(line[1])
        for ex_cnt in range(int(line[-3])):
            start = st + int(rstart[ex_cnt]) + 1
            stop = start + int(exon_len[ex_cnt]) - 1
            
            if ex_cnt > 0:
                pline = [str(line[0]),
                        source_name,
                        'intron',
                        str(intron_start),
                        str(start-1),
                        line[4],
                        line[5],
                        '.',
                        'Parent=' + line[3]]
                print '\t'.join(pline) 
                
            pline = [str(line[0]),
                        source_name,
                        'exon',
                        str(start),
                        str(stop),
                        line[4],
                        line[5],
                        '.',
                        'Parent=' + line[3]]
            print '\t'.join(pline) 
            intron_start = stop+1

    BEDfh.close()

def __main__():
   
    source_name = 'bed2gff'

    try:
        query_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    bed_parse(query_file, source_name)

if __name__ == "__main__": 
    __main__()
