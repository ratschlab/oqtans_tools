"""

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2009-2010 Andre Kahles
  Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

  This program extracts intron features from a given annotation in GFF3 format
  
  For detailed usage information type:

    python gen_intronlist_from_annotation.py

"""

import GFFParser
import cPickle
import sys
from Bio.SeqFeature import SeqFeature
import pdb

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--annotation', dest='anno', metavar='FILE', help='annotation file in gff3 format', default='-')
    required.add_option('-o', '--output', dest='outfile', metavar='FILE', help='annotation intron list', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-S', '--show_sources', dest='show_sources', action='store_true', help='only show available sources of gff file', default=False)
    optional.add_option('-s', '--sources', dest='sources', metavar='SOURCELIST', help='list of comma-separated sources to use from annotation', default='')
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 3 or '-' in [options.anno, options.outfile] :
        parser.print_help()
        sys.exit(2)

    return (options, args)


def fix_structure(trans_dict):
    """Adapts the structure of the trans_dict to different gff-file structures"""

    for idx in range(len(trans_dict.features)):
        old_trans_dict = trans_dict.features[idx]
        if trans_dict.features[idx].type == 'exon':
            exon = trans_dict.features[idx]
            trans_dict.features[idx] = SeqFeature(exon.location, type = 'gene', strand = exon.strand, id = exon.id)
            trans_dict.features[idx].sub_features = [SeqFeature(exon.location, type = 'Transcript', strand = exon.strand, id = exon.id)] 
            trans_dict.features[idx].sub_features[0].sub_features = [exon]
        elif len(trans_dict.features[idx].sub_features) > 0 and trans_dict.features[idx].sub_features[0].type == 'exon':
            exon = trans_dict.features[idx]
            trans_dict.features[idx] = SeqFeature(exon.location, type = 'gene', strand = exon.strand, id = exon.id)
            trans_dict.features[idx].sub_features = [exon]

def main():
    """Main function handling the program flow and organizing the GFF parser."""

    (options, args) = parse_options(sys.argv)

    iterator = GFFParser.GFFAddingIterator()    
    examiner = GFFParser.GFFExaminer()

    exon_map = dict()

    id_dict = examiner.available_limits(options.anno)['gff_id']
    intron_lists = dict()

    ### collect all available sources from gff-file
    source_dict = examiner.available_limits(options.anno)['gff_source_type']
    taken_sources = set()
    #types = ['gene', 'mRNA', 'exon', 'CDS']
    types = ['exon']

    ### parse only for exons and let the GFFparser 
    ### infer the respective parents (otherwise doubled entries occured)
    ### we sanitize the structure later on anyways
    for key in [source[0] for source in source_dict.keys() if source[1] in types]:
        taken_sources.add(key)

    ### try different type, if sources are empty    
    if len(taken_sources) == 0:
        types = ['CDS']
        for key in [source[0] for source in source_dict.keys() if source[1] in types]:
            taken_sources.add(key)

    ### print taken_sources
    if len(taken_sources) == 0:
        print >> sys.stderr, 'No suitable sources found!'
        sys.exit(-1)

    ### only show available sources - if neccessary
    if options.show_sources:
        print 'Parsed file %s\n' % options.anno
        print 'Following sources are available:\n'
        for source in taken_sources:
            print source    
        print '\nUse option -s to specify a comma-separated list of sources (-s source1,source2,source3), otherwise all sources are taken'
        sys.exit(0)

    if options.sources != '':
        user_sources = set(options.sources.split(','))
        taken_sources = taken_sources.intersection(user_sources)
        if len(taken_sources) == 0:
            print >> sys.stderr, 'The specified sources do not match any of the available sources - Please use option -S to get a list of available sources'
            sys.exit(-1)

    if options.verbose:
        print "take sources %s" % str(list(taken_sources))

    ### build up gff-parsing filter
    gff_sources = []
    for source in taken_sources:
        gff_sources.extend(zip([source] * len(types), types))

    ### parse gff-file
    for idx in id_dict.keys():
        print 'parsing chromosome %s' % idx
        if len(gff_sources) > 0:
            trans_dict = iterator.get_all_features(options.anno, {'gff_source_type':gff_sources, 'gff_id':idx})
        else:
            trans_dict = iterator.get_all_features(options.anno, {'gff_id':idx})
        ### since we parse only one chromosome, this loop is evaluated only once
        for chrm in trans_dict.keys():
            ### verify/sanitize the created dictionairy
            fix_structure(trans_dict[chrm])
            intron_lists[chrm] = dict()
            for gene in trans_dict[chrm].features:
                for trans in gene.sub_features:
                    if trans.type == 'exon':
                        print "WARNING: Exon on transcript level:"
                        print trans
                        print 'will continue\n'
                        continue
                    elif len(trans.sub_features) > 1: ### at least two exons for one intron ...
                        strand = trans.sub_features[0].strand
                        contig_list = [(trans.sub_features[i].location.nofuzzy_start, trans.sub_features[i].location.nofuzzy_end) for i in range(len(trans.sub_features))]
                        contig_list.sort(lambda u, v:u[0]-v[0])
                        for exon in range(len(contig_list) - 1):
                            ### update intron lists
                            if contig_list[exon][1] - contig_list[exon + 1][0] == 0:
                                continue
                            try:
                                assert(contig_list[exon][1] < contig_list[exon + 1][0])
                            except AssertionError:
                                print >> sys.stderr, 'exon_1 %i, exon_2 %i' % (contig_list[exon][1], contig_list[exon + 1][0]) 
                                print >> sys.stderr, contig_list[exon]
                                print >> sys.stderr, contig_list[exon+1]
                                print >> sys.stderr, exon
                                sys.exit(-1)
                            ### for now strand information is only dummy
                            intron_lists[chrm][(0, contig_list[exon][1], contig_list[exon + 1][0])] = strand
                     
                        ### update exon map
                        for exon in range(len(contig_list)):
                            if not exon_map.has_key(chrm):
                                exon_map[chrm] = dict()

                            if not exon_map[chrm].has_key(trans.id):
                                exon_map[chrm][trans.id] = dict()
                            ### we assume, that an exon cannot occurr twice in the same transcript!
                            ### the value in the dict is a binary encoding, if the left/right end is intronic 10 = 2 means, 5' end is intronic
                            if len(contig_list) == 1:
                                exon_map[chrm][trans.id][contig_list[exon]] = 0 ### 00 -> should never occurr
                            elif exon == 0:
                                exon_map[chrm][trans.id][contig_list[exon]] = 2 ### 10
                            elif exon == len(contig_list) - 1:
                                exon_map[chrm][trans.id][contig_list[exon]] = 1 ### 01
                            else:
                                exon_map[chrm][trans.id][contig_list[exon]] = 3 ### 11 

    outfile = open(options.outfile, 'w')
    cPickle.dump(intron_lists, outfile)
    outfile.close()
    
    outfile = open(options.outfile + '.' + 'cov', 'w')
    cPickle.dump(exon_map, outfile)
    outfile.close()

if __name__ == '__main__':
    main()

