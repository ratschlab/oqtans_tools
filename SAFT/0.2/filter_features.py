"""

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2095-2010 Andre Kahles
  Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

  This script takes an intron feature file as input and filters it according
  to given filter criteria.
  
  For detailed usage information type:

    python filter_alignment.py 
"""


import sys
import re

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-i', '--infile', dest='infile', metavar='FILE', help='intron feature file', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    optional.add_option('-X', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismathes [10000]', default=10000)
    optional.add_option('-m', '--min_support', dest='min_support', type='int', metavar='INT', help='minimum support for a feature to be counted [default 1]', default=1)
    optional.add_option('-M', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [100000000]', default='100000000')
    optional.add_option('-o', '--outfile', dest='outfile', metavar='FILE', help='outfile name - default: tagged infile name', default='-')
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def main():
    """Main function extracting intron features."""

    options = parse_options(sys.argv)

    ### set filter tags for filenames
    if options.outfile == '-':
        outfile_base = options.infile
        if options.min_exon_len > 0:
            outfile_base += '_me%s' % options.min_exon_len
        if options.max_mismatches < 10000:
            outfile_base += '_mm%s' % options.max_mismatches
        if options.max_intron_len < 100000000:
            outfile_base += '_mi%s' % options.max_intron_len
        if options.min_support > 1:
            outfile_base += '_mc%s' % options.min_support
        outfile_base = (re.sub('.features', '', outfile_base) + '.features')
    else:
        outfile_base = options.outfile
    
    outfile = open(outfile_base, 'w')

    line_counter = 0
    filter_counter = 0
    for line in open(options.infile, 'r'):
        if options.verbose and line_counter % 10000 == 0:
            print 'processed %i features from %s - filtered %i - kept %i' % (line_counter, options.infile, filter_counter, line_counter - filter_counter)
        line_counter += 1
        sl = line.strip().split('\t')
        (start, stop) = sl[1:3]
        if int(stop) - int(start) > options.max_intron_len:
            filter_counter += 1
            continue

        if int(sl[3]) < options.min_support:
            filter_counter += 1
            continue

        cont_flag = True
        for _sl in sl[5:]:
            key = _sl.split(':')[0]
            if int(key) / 81 >= options.min_exon_len and int(key) % 81 <= options.max_mismatches:
                cont_flag = False
                break
        if cont_flag:
            filter_counter += 1
            continue
        
        print >> outfile, line,

    outfile.close()

if __name__ == "__main__":
    main()
