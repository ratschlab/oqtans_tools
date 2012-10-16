"""

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2009-2010 Andre Kahles
  Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

  This program extracts intron features from a given alignment.
  
  For detailed usage information type:

    python get_intron_featuress.py 

"""

import sys
import re
import subprocess

class Feature(object):
    """Is an intron feature object"""


    def __init__(self, max_mm=80):
        self.alignment_support = 0
        self.submission_support = 1
        self.mm_ex = dict()
        self.max_mm = max_mm + 1

    def add_mm_ex(self, ex, mm):
        """Adds mm ex information"""

        self.alignment_support += 1
        try:
            self.mm_ex[(ex*self.max_mm) + mm] += 1
        except KeyError:
            self.mm_ex[(ex*self.max_mm) + mm] = 1

    def get_feature_string(self):
        """Returns string with mm ex elements."""
        _line = (str(self.alignment_support) + '\t' + str(self.submission_support) + '\t')
        for key in self.mm_ex:
            _line += (str(key) + ':' + str(self.mm_ex[key]) + '\t')

        return _line[:-1]

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--alignment', dest='alignment', metavar='FILE', help='alignment file in sam format', default='-')
    required.add_option('-o', '--outfile', dest='outfile', metavar='FILE', help='intron feature file', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-b', '--bam_input', dest='bam_input', action='store_true', help='input has BAM format - does not work for STDIN', default=False)
    optional.add_option('-s', '--samtools', dest='samtools', metavar='PATH', help='if SAMtools is not in your PATH, provide the right path here (only neccessary for BAM input)', default='')
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

    ### initialize counters

    line_counter = 0
    skipped_lines = 0

    ### prepare data structures

    features = dict()

    if options.outfile == '-':
        outfile = sys.stdout
    else:
        outfile = open(options.outfile, 'w')

    if options.bam_input:
        infile_handles = []
        if options.samtools == '':
            options.samtools = 'samtools'
        else:
            options.samtools = '%s/samtools' % options.samtools

    if options.alignment == '-':
        infiles = [sys.stdin]
    else:
        if options.alignment.find(',') > -1:
            infiles = []
            for inf in options.alignment.strip().split(','):
                if inf == '':
                    continue
                if options.bam_input:
                    infile_handles.append(subprocess.Popen([options.samtools, 'view', inf], stdout=subprocess.PIPE))
                    infiles.append(infile_handles[-1].stdout)
                else:
                    infiles.append(open(inf, 'r'))
        else:       
            if options.bam_input:
                infile_handles.append(subprocess.Popen([options.samtools, 'view', options.alignment], stdout=subprocess.PIPE))
                infiles = [infile_handles[-1].stdout]
            else:
                infiles = [open(options.alignment, 'r')]

    for infile in infiles:
        for line in infile:
            sl = line.strip().split('\t')

            if options.verbose and line_counter % 10000 == 0:
                print '[ lines read: %i / skipped: %i ]\r' % (line_counter, skipped_lines),
            line_counter += 1

            if len(sl) < 9:
                sl = line.strip().split()
            if len(sl) < 9:
                skipped_lines += 1
                continue
            if sl[5].find('N') == -1:
                skipped_lines += 1
                continue
        
            (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
            size = [int(i) for i in size]
            
            start = int(sl[3]) - 1
            chrm = sl[2]
            chrm = chrm.replace('chr','')
            chrm = chrm.replace('Chr','')

            min_ex_len = 100000

            __cig = sl[5]
            __cig = re.sub('[0-9]*[IHS]', '', __cig) 
            for _cig in __cig.strip().split('N'):
                min_ex_len = min(sum([int('0'+i) for i in re.split('[^0-9]', '0' + _cig + 'Z0')][:-2]), min_ex_len)

            mm = -1
            for tag in sl[11:]:
                if tag[:2] == 'NM':
                    mm = int(tag[5:])
                    break

            if mm == -1:    
                print >> sys.stderr, 'Mismatch information is missing in %s' % options.alignment
                sys.exit(-1)

            offset = 0
            for o in range(len(op)):
                if op[o] == 'N':
                    istart = start + offset
                    iend = istart + size[o]
                    try:
                        features[(chrm, istart, iend)].add_mm_ex(min_ex_len, mm)
                    except KeyError:
                        features[(chrm, istart, iend)] = Feature()
                        features[(chrm, istart, iend)].add_mm_ex(min_ex_len, mm)
                if not op[o] in ['I', 'H', 'S']:
                    offset += size[o]

        infile.close()
        if options.bam_input:
            infile_handles[infiles.index(infile)].kill()

    ### store statistics
    for feature in features.keys():
        _line = (str(feature[0]) + '\t' + str(feature[1]) + '\t' + str(feature[2]) + '\t')
        _line += features[feature].get_feature_string()
        print >> outfile, _line

    if options.outfile != '-':
        outfile.close()

if __name__ == "__main__":
    main()
