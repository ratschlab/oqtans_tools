"""

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2095-2010 Andre Kahles
  Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

  This script compares an alignment to the annotation.
  
  For detailed usage information type:

    python compare_to_features.py 

"""

import cPickle
import sys

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED INPUTS')
    required.add_option('-i', '--annotation_introns', dest='anno_int', metavar='FILE', help='annotation intron list', default='-')
    required.add_option('-f', '--features', dest='features', metavar='FILE', help='alignment intron features', default='-')
    filtering = OptionGroup(parser, 'FILTER SETTINGS')
    filtering.add_option('-S', '--splice_consensus', dest='splice_consensus', action='store_true', help='enforce introns to have splice consensus', default=False)
    filtering.add_option('-I', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [100000000]', default='100000000')
    filtering.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    filtering.add_option('-M', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismacthes [1000000]', default=1000000)
    filtering.add_option('-c', '--min_coverage', dest='min_coverage', metavar='INT', type='int', help='minimal coverage to count an intron [0]', default='0')
    filtering.add_option('-C', '--min_coverage_2', dest='min_coverage_2', metavar='INT', type='int', help='minimal coverage after applying all other filters [1]', default='1')
    filtering.add_option('-x', '--exclude_chrms', dest='exclude_chrm', metavar='STRINGLIST', help='list of comma separated chromosomes to exclude from evaluation', default='-')
    filtering.add_option('-E', '--exclude_introns', dest='exclude_introns', metavar='STRINGLIST', help='list of comma separated intron files to exclude from submitted features', default='-')
    iosettings = OptionGroup(parser, 'I/O SETTINGS')
    iosettings.add_option('-s', '--strand_specific', dest='strand_specific', action='store_true', help='data is strand specific', default=False)
    iosettings.add_option('-p', '--performance_log', dest='performance_log', action='store_true', help='store the intron recovery performance in extra log [off]', default=False)
    iosettings.add_option('-o', '--outfile_base', dest='outfile_base', metavar='PATH', help='basedir for performance log', default='-')
    optional = OptionGroup(parser, 'OTHERS')
    optional.add_option('-X', '--max_feat_mismatches', dest='max_feat_mismatches', metavar='INT', type='int', help='max number of mismatches for feat generation [80] (only change this, if you are absolutely sure!)', default=80)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(filtering)
    parser.add_option_group(iosettings)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 3 or '-' in [options.features, options.anno_int] :
        parser.print_help()
        sys.exit(2)

    return options

def build_intron_list(options):
    """Builds up an intron list from the given alignment features."""

    counter = 0
    filter_counter = 0
    intron_lists = dict()

    for line in open(options.features, 'r'):
        
        if counter % 10000 == 0 and options.verbose:
            print 'lines read: [ %s (taken: %s / filtered: %s)]' % (counter, counter - filter_counter, filter_counter)
        counter += 1

        sl = line.strip().split('\t')
        chrm = sl[0]

        ### find valid intron
        if int(sl[3]) < options.min_coverage:
            filter_counter += 1
            continue 

        if options.splice_consensus:
            ### load genome if necessary
            if not options.genome_dict.has_key(chrm):
                if chrm in options.gio.contig_names:
                    if options.verbose:
                        print 'Reading chromosome %s' % chrm
                    fgen = open('%s/genome/%s.flat' % (options.gio.basedir, chrm))
                    options.genome_dict[chrm] = fgen.readline()
                    fgen.close()
                elif options.ignore_missing_chrm:
                    continue
                else:       
                    print >> sys.stderr, 'Chromosome Names do not match genome information file. Chromosome %s can not be found!' % chrm
                    exit(2)



        for _sl in sl[5:]:
            (val, count) = _sl.split(':')
            ex = int(val) / (options.max_feat_mismatches + 1)
            mm = int(val) % (options.max_feat_mismatches + 1)
            if ex >= options.min_exon_len and mm <= options.max_mismatches:
                ### intron = (strand [-1, 1], start, end)  ... end in the pythonic way
                intron = (0, int(sl[1]), int(sl[2]))
                try:
                    intron_lists[chrm][intron] += int(count)
                except KeyError:
                    try:
                        intron_lists[chrm][intron] = int(count)
                    except KeyError:
                        intron_lists[chrm] = {intron:int(count)}
            else:
                filter_counter += 1
            

    print '\n parsed file %s' % options.features
    print 'read %s lines' % counter
    if options.max_mismatches < 1000000 or options.min_exon_len > 0 or options.min_coverage > 1:
        print 'filter criteria:'
        print '    max mismatches: %s' % options.max_mismatches
        print '    min exon length: %s' % options.min_exon_len
        print '    min coverage: %s\n' % options.min_coverage
        print 'filtered %s lines' % filter_counter

    return intron_lists


def main():
    """Main function for comparing two intron lists """    
    
    options = parse_options(sys.argv)

    if options.outfile_base != '-':
        outfile_base = (options.outfile_base + '/' + options.features.split('/')[-1])
    else:
        outfile_base = options.features

    ### set filter tags for filenames
    if options.min_exon_len > 0:
        outfile_base += '_me%s' % options.min_exon_len
    if options.max_mismatches < 1000000:
        outfile_base += '_mm%s' % options.max_mismatches
    if options.max_intron_len < 100000000:
        outfile_base += '_mi%s' % options.max_intron_len
    if options.min_coverage > 0:
        outfile_base += '_mc%s' % options.min_coverage

    ### load genome, if necessary
    if options.splice_consensus:
        import genome_utils
        options.gio = genome_utils.GenomeInfo(options.genome)
        options.gio.contig_names.sort()
        options.genome_dict = dict()

    ### parse feature file for creating intron_list and coverage_map
    alignment_list = build_intron_list(options)
    
    ### get list of annotated introns
    annotation_list = cPickle.load(open(options.anno_int, 'r'))

    ### filter annotated and predicted introns for excluded chromosomes
    if options.exclude_chrm != '-':
        _ex_chrm = options.exclude_chrm.strip().split(',')
        ### handle leading or trailing commas
        if _ex_chrm[0] == '':
            _ex_chrm = _ex_chrm[1:]
        if _ex_chrm[-1] == '':
            _ex_chrm = _ex_chrm[:-1]
        for chrm in _ex_chrm:
            if annotation_list.has_key(chrm):
                del annotation_list[chrm]

            if alignment_list.has_key(chrm):
                del alignment_list[chrm]

    ### filter predicted introns for excluded introns
    if options.exclude_introns != '-':
        _ex_introns = options.exclude_introns.strip().split(',')
        ### handle leading or trailing commas
        if _ex_introns[0] == '':
            _ex_introns = _ex_introns[1:]
        if _ex_introns[-1] == '':
            _ex_introns = _ex_introns[:-1]
        for _infile in _ex_introns:
            _ex_intron = cPickle.load(open(_infile, 'r'))
            for chrm in _ex_intron.keys():
                if alignment_list.has_key(chrm):
                    for _intron in _ex_intron[chrm].keys():
                        try:
                            del alignment_list[chrm][_intron]
                        except KeyError:
                            continue

    ### filter intron lists for max intron length
    print '\nFiltering intron list for max intron len'
    print '-----------------------------------------'
    skipped = 0
    for chrm in annotation_list.keys():
        skiplist = set()
        for intron in annotation_list[chrm].keys():
            if (intron[2] - intron[1]) > options.max_intron_len:
                skiplist.add(intron)
        for intron in skiplist:
            del annotation_list[chrm][intron]
        skipped += len(skiplist)
    print '%s introns removed from annotation' % skipped

    skipped = 0
    for chrm in alignment_list.keys():
        skiplist = set()
        for intron in alignment_list[chrm].keys():
            if (intron[2] - intron[1]) > options.max_intron_len:
                skiplist.add(intron)
        for intron in skiplist:
            del alignment_list[chrm][intron]
        skipped += len(skiplist)
    print '%s introns removed from alignment' % skipped
    del skiplist

    ### filter intron lists for min coverage
    if options.min_coverage_2 > 1:
        print '\nFiltering intron list for min support after filtering'
        print '-----------------------------------------------------'
        skipped = 0
        for chrm in alignment_list.keys():
            skiplist = set()
            for intron in alignment_list[chrm].keys():
                if alignment_list[chrm][intron] < options.min_coverage_2:
                    skiplist.add(intron)
            for intron in skiplist:
                del alignment_list[chrm][intron]
                if annotation_list[chrm].has_key(intron):
                    del annotation_list[chrm][intron]
            skipped += len(skiplist)
        print '%s introns removed from alignment\n' % skipped
        del skiplist

    ### match intron lists
    non_matched = dict()
    total_precision = float(0)
    total_recall = float(0)
    key_count = 0
    for chrm in annotation_list.keys():
        if alignment_list.has_key(chrm):
            matches = len(set(annotation_list[chrm].keys()).intersection(set(alignment_list[chrm].keys())))
            _curr_recall = float(matches) / float(max(1, len(annotation_list[chrm].keys())))
            _curr_precision = float(matches) / float(max(1, len(alignment_list[chrm].keys())))
            _curr_fscore = (2 * _curr_precision * _curr_recall) / float(max(1, _curr_precision + _curr_recall))
            print '-----------------------------'
            print ' in Chromosome %s ' % chrm 
            print '-----------------------------'
            print 'recall: %s' % _curr_recall
            print 'precision: %s' % _curr_precision
            print 'F-Score: %s' % _curr_fscore
            total_precision += _curr_precision
            total_recall += _curr_recall
            non_matched[chrm] = set(alignment_list[chrm].keys()).difference(set(annotation_list[chrm].keys()))
            ### do not include chromosomes with zero values into average
            if matches > 0:
                key_count += 1
   
    total_precision /= max(1.0, float(key_count))
    total_recall /= max(1.0, float(key_count))
    total_fscore = (2 * total_precision * total_recall) / float(max(1, total_precision + total_recall))
    print '-----------------------------'
    print ' average over all chromosomes '
    print '-----------------------------'
    print 'recall: %s' % total_recall
    print 'precision: %s' % total_precision
    print 'F-Score: %s' % total_fscore

    if options.performance_log:
        outf = open(outfile_base + '_performance.log', 'w')
        non_matched = dict()
        total_precision = float(0)
        total_recall = float(0)
        key_count = 0
        _recall_line = ''
        _precision_line = ''
        _header_line = ''
        for chrm in annotation_list.keys():
            if alignment_list.has_key(chrm):
                matches = len(set(annotation_list[chrm].keys()).intersection(set(alignment_list[chrm].keys())))
                _header_line += (str(chrm) + '\t')
                _recall_line += (str(float(matches) / float(max(1, len(annotation_list[chrm].keys())))) + '\t')
                _precision_line += (str(float(matches) / float(max(1, len(alignment_list[chrm].keys())))) + '\t')
                total_precision += (float(matches) / float(max(1, len(alignment_list[chrm].keys()))))
                total_recall += (float(matches) / float(max(1, len(annotation_list[chrm].keys()))))
                non_matched[chrm] = set(alignment_list[chrm].keys()).difference(set(annotation_list[chrm].keys()))
                ### do not include chromosomes with zero values into average
                if matches > 0:
                    key_count += 1
        total_precision /= max(1.0, float(key_count))
        total_recall /= max(1.0, float(key_count))
        _header_line += 'average'
        _recall_line += str(total_recall)
        _precision_line += str(total_precision)
        print >> outf, _header_line
        print >> outf, _recall_line
        print >> outf, _precision_line
        outf.close() 
    
if __name__ == '__main__':
    main()
