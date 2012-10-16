#!/usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2010 Vipin T Sreedharan
# Copyright (C) 2010 Max Planck Society
#

# Description : Provide available source, feature types from a GFF file 

import re, sys
import time 
import collections

def available_limits(gff_handle):
    """Figure out the available feature types from the given GFF file"""

    filter_info = dict(gff_id = [0], gff_source_type = [1, 2],
                gff_source = [1], gff_type = [2])
    cur_limits = dict()
    for filter_key in filter_info.keys():
        cur_limits[filter_key] = collections.defaultdict(int)
    for line in gff_handle:
        if line.strip('\n\r')[0] != "#":
            parts = [p.strip() for p in line.split('\t')]
            if len(parts) == 1 and re.search(r'\w+', parts[0]):continue ## GFF files with FASTA sequence together 
            assert len(parts) == 9, line
            for filter_key, cur_indexes in filter_info.items():
                cur_id = tuple([parts[i] for i in cur_indexes])
                cur_limits[filter_key][cur_id] += 1
    # get rid of the default dicts
    final_dict = dict()
    for key, value_dict in cur_limits.items():
        if len(key) == 1:
            key = key[0]
        final_dict[key] = dict(value_dict)
   
    return final_dict

if __name__=='__main__':
    
    stime = time.asctime( time.localtime(time.time()) )
    print '-------------------------------------------------------'
    print 'FeatureScan started on ' + stime
    print '-------------------------------------------------------'

    try:
        gff_handle = open(sys.argv[1], 'rU')
    except:
        sys.stderr.write("Can't open the GFF3 file, terminating...\n")
        sys.stderr.write("USAGE: gff_available_limits.py <gff file>\n")
        sys.exit(-1)
    final_dict = available_limits(gff_handle)
    gff_handle.close()
    print 
    print "==Overview of available source(s) and feature type(s) from GFF file=="
    print 
    print "Chromosome identifier(s) and corresponding count:"
    for contig, cnt in sorted(final_dict['gff_id'].items()):
        print '\t' + str(contig[0]) + '\t' + str(cnt)
    print
    print "Source(s) of feature and corresponding count:"
    for source, cnt in sorted(final_dict['gff_source'].items()):
        print '\t' + str(source[0]) + '\t' + str(cnt)
    print
    print "Feature type(s) and corresponding count:"
    for ftype, cnt in sorted(final_dict['gff_type'].items()):
        print '\t' + str(cnt) + '\t' + str(ftype[0]) 
    print
    print "Unique combination of Feature type(s), Source(s) and corresponding count:"
    for sftype, cnt in sorted(final_dict['gff_source_type'].items()):
        print '\t' + str(cnt) + '\t' + str(sftype[0]) + ', '+ str(sftype[1])
    print
    stime = time.asctime( time.localtime(time.time()) )
    print '-------------------------------------------------------'
    print 'FeatureScan finished at ' + stime
    print '-------------------------------------------------------'
