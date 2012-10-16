"""

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2095-2010 Andre Kahles
  Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

  This script filters any given alignment by given criteria.
  
  For detailed usage information type:

    python filter_alignment.py 
"""

import sys
import re
import subprocess

class Feature(object):
    """Is an intron feature object"""


    def __init__(self, max_mm=80, feature_string=''):
        if feature_string == '':
            self.alignment_support = 0
            self.submission_support = 1
            self.mm_ex = dict()
            self.max_mm = max_mm + 1
        else:
            self.alignment_support = int(feature_string[0])
            self.submission_support = int(feature_string[1])
            self.mm_ex = dict()
            self.max_mm = max_mm + 1
            for _sl in feature_string[2:]:
                (key, value) = _sl.split(':')
                self.mm_ex[key] = int(value)

    def merge_features(self, feature_string):
        """Merges information in feature_string into current feature object"""
        self.alignment_support += int(feature_string[0])
        self.submission_support += int(feature_string[1])
        for _sl in feature_string[2:]:
            (key, value) = _sl.split(':')
            try:
                self.mm_ex[key] += int(value)
            except KeyError:
                self.mm_ex[key] = int(value)

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


    def get_submission_support(self):
        """Returns  submission support"""

        return int(self.submission_support)


def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    io_opts = OptionGroup(parser, 'I/O OPTIONS')
    io_opts.add_option('-a', '--alignment', dest='align', metavar='FILE', help='alignment file in sam format (- for stdin)', default='-')
    io_opts.add_option('-o', '--outfile', dest='outfile', metavar='FILE', help='outfile - default is tagged infile or stdout for stdin processing', default='')
    io_opts.add_option('-b', '--bam_input', dest='bam_input', action='store_true', help='input has BAM format - does not work for STDIN', default=False)
    io_opts.add_option('-s', '--samtools', dest='samtools', metavar='PATH', help='if SAMtools is not in your PATH, provide the right path here (only neccessary for BAM input)', default='')
    filter_crit = OptionGroup(parser, 'FILTER CRITERIA')
    filter_crit.add_option('-R', '--ignore_multireads', dest='multireads', metavar='FILE', help='file containing the multireads to ignore', default='-')
    filter_crit.add_option('-M', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [100000000]', default='100000000')
    filter_crit.add_option('-e', '--min_exon_len', dest='min_exon_len', metavar='INT', type='int', help='minimal exon length [0]', default=0)
    filter_crit.add_option('-X', '--max_mismatches', dest='max_mismatches', metavar='INT', type='int', help='maximum number of allowed mismathes [10000]', default=10000)
    filter_crit.add_option('-i', '--intron_features', dest='intron_features', metavar='FILE', help='intron features file - only spliced reads present in this file are kept', default='-')
    filter_crit.add_option('-c', '--clip_filter', dest='clip_filter', action='store_true', help='filters clipped reads', default=False)
    filter_hand = OptionGroup(parser, 'FILTER HANDLING')
    filter_hand.add_option('-d', '--del_worse', dest='del_worse', action='store_true', \
        help='all alignments of the same read with more mismatches are deleted as well. File must be sorted by read ID!', default=False)
    filter_hand.add_option('-n', '--no_suboptimal', dest='no_suboptimal', action='store_true', \
        help='if suboptimal spliced alignments occur and a batter alignment is available, keep only the best alignment. File must be sorted by read id and flag', default=False)
    filter_hand.add_option('-w', '--window', dest='window', metavar='INT', type='int', help='size of overlap-window to count suboptimal alignments as one stratum [2]', default=2)
    others = OptionGroup(parser, 'OTHERS')
    others.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(io_opts)
    parser.add_option_group(filter_crit)
    parser.add_option_group(filter_hand)
    parser.add_option_group(others)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def filter_lines_by_criteria(curr_lines, multireads, options, filter_counter):
    """Filters the list of lines by given criteria"""

    min_mm = 100000
    _curr_lines = []

    for sl in curr_lines:
        mm = None 
        if options.del_worse:
            try:
                for opt in sl[11:]:
                    if opt[:3] == 'NM:':
                        mm = int(opt[5:])
                        break
                if mm == None:
                    print >> sys.stderr, 'No mismatch information available or read string missing in %s' % options.align
                    sys.exit(1)
            except IndexError:
                print >> sys.stderr, 'No mismatch information available or read string missing in %s' % options.align
                sys.exit(1)

        if sl[5].find('N') > -1 and options.intron_features != '-':
            cont_flag = True
            (op, size) = (re.split('[0-9]*', sl[5])[1:], re.split('[^0-9]', sl[5])[:-1])
            size = [int(i) for i in size]
            offset = 0
            for o in range(len(op)):
                if op[o] in ['M', 'D']:
                    offset += size[o]
                if op[o] == 'N':
                    start = int(sl[3]) - 1 + offset
                    stop = start + size[o]
                    if options.features.has_key((sl[2], str(start), str(stop))):
                        cont_flag = False
                        break
            if cont_flag:
                filter_counter += 1
                if mm != None:
                    min_mm = min(mm, min_mm)
                continue

        if (sl[0], int(sl[1]) & 196) in multireads:
            filter_counter += 1
            if mm != None:
                min_mm = min(mm, min_mm)
            continue

        if options.max_mismatches < 10000:
            cont_flag = False
            if options.del_worse and mm > options.max_mismatches:
                filter_counter += 1
                if mm != None:
                    min_mm = min(mm, min_mm)
                continue
            else:
                try:
                    for opt in sl[11:]:
                        if (opt[:3] == 'NM:' and int(opt[5:]) > options.max_mismatches):
                            cont_flag = True
                            break
                    if cont_flag:
                        filter_counter += 1
                        if mm != None:
                            min_mm = min(mm, min_mm)
                        continue
                except IndexError:
                    print >> sys.stderr, 'No mismatch information available or read string missing in %s' % options.align
                    sys.exit(1)

        if options.min_exon_len > 0:
            cont_flag = False
            __cig = sl[5]
            __cig = re.sub('[0-9]*[IHS]', '', __cig) 
            for _cig in __cig.strip().split('N'):
                if sum([int('0'+i) for i in re.split('[^0-9]', '0' + _cig + 'Z0')][:-2]) < options.min_exon_len:
                    cont_flag = True
                    break
            if cont_flag:
                filter_counter += 1
                if mm != None:
                    min_mm = min(mm, min_mm)
                continue

        if options.max_intron_len < 100000000:
            cont_flag = False
            (op, size) = (re.split('[0-9]*', sl[5])[1:], re.split('[^0-9]', sl[5])[:-1])
            size = [int(i) for i in size]

            for o in range(len(op)):
                if op[o] == 'N' and (size[o] > options.max_intron_len):
                    cont_flag = True
                    break
            if cont_flag:
                filter_counter += 1
                if mm != None:
                    min_mm = min(mm, min_mm)
                continue

        if options.clip_filter and (sl[5].find('H') > -1 or sl[5].find('S') > -1):
            filter_counter += 1
            if mm != None:
                min_mm = min(mm, min_mm)
            continue

        _curr_lines.append(sl)

    return (min_mm, _curr_lines, filter_counter)

def filter_lines_by_min_mm(curr_lines, min_mm, filter_counter):
    """Filters the given lines by given number of max mismatches """
    
    _curr_lines = []
    for sl in curr_lines:
        for opt in sl[11:]:
            if opt[:3] == 'NM:':
                mm = int(opt[5:])
                break

        if mm > min_mm:
            filter_counter += 1
            continue
        _curr_lines.append(sl)

    return (_curr_lines, filter_counter)

def filter_suboptimal(curr_lines, filter_counter, options):
    """Removes suboptimal spliced alignments, iff better alignment exists."""
    

    ### check, if QPALMA Score is available
    have_spliced = False
    scored_lines = []
    for cline in curr_lines:
        if cline[5].find('N') > -1:
            have_spliced = True
        qp_score = False
        for idx in range(12, len(cline)):
            if cline[idx][:5] == 'ZS:f:':
                scored_lines.append((cline, float(cline[idx][5:])))
                qp_score = True
                break
        if not qp_score:
            break

    ### ignore reads that have not any spliced alignment
    if not have_spliced:
        return (curr_lines, filter_counter)

    if not qp_score:
        print >> sys.stderr, 'ERROR: No QPalma Score for sorting available! Bailing out! - This option is only available for PALMapper alignments'
        sys.exit(2)

    ### sort lines by starting position
    sorted_lines = sorted(scored_lines, key=lambda cline: int(cline[0][3]))
    ### sort lines by chromosome, keep starting position sortage since sort is stable
    sorted_lines = sorted(sorted_lines, key=lambda cline: cline[0][2])
   
    ### form strata
    strata = []
    last_start = -1000
    last_stop = -1000
    last_chrm = 'XXXX'
    for sline in sorted_lines:
        start = int(sline[0][3])
        _sl = re.sub('[0-9]+[IHS]', '', sline[0][5])
        stop = start + sum([int(i) for i in re.split('[^0-9]', _sl)[:-1]])
        if ((start >= last_start - options.window and start <= last_start + options.window) or (stop >= last_stop - options.window and stop <= last_stop + options.window)) \
         and sline[0][2] == last_chrm:
            strata[-1].append(sline)
        else:
            strata.append([sline])
        last_start = start
        last_stop = stop
        last_chrm = sline[0][2]

    ### sort strata by QPalma-score and check if we keep it 
    for stratum in strata:
        stratum = sorted(stratum, key=lambda cline: cline[1], reverse=True)
        ### keep spliced alignments only, if they are first in the respective stratum
        for idx in range(len(stratum)):
            if stratum[idx][0][5].find('N') > -1:
                if idx != 0:
                    curr_lines.remove(stratum[idx][0])
                    filter_counter += 1

    return (curr_lines, filter_counter)

def main():
    """Main function for filtering alignments in SAM format"""    
    
    options = parse_options(sys.argv)

    multireads = set()
    if options.multireads != '-':
        multireads = []
        print '\nParsing multireads from file %s' % options.multireads
        print '-----------------------------------------'
        for line in open(options.multireads, 'r'):
            try:
                _l = line.strip().split('\t')
                multireads.append((_l[0], int(_l[1])))
            except IndexError:
                print >> sys.stderr, 'skipped line %s' % line
                continue
        multireads = set(multireads)
    counter = 0
    filter_counter = 0

    if options.bam_input:
        if options.samtools == '':
            options.samtools = 'samtools'
        else:
            options.samtools = '%s/samtools' % options.samtools

    if options.outfile == '':
        if options.align != '-':
            outfile_base = options.align

            ### set filter tags for filenames
            if options.min_exon_len > 0:
                outfile_base += '_me%s' % options.min_exon_len
            if options.max_mismatches < 10000:
                outfile_base += '_mm%s' % options.max_mismatches
            if options.max_intron_len < 100000000:
                outfile_base += '_mi%s' % options.max_intron_len
            if options.multireads != '-':
                outfile_base += '_mrFiltered'
            if options.clip_filter:
                outfile_base += '_clipFiltered'
            if options.no_suboptimal:
                outfile_base += '_noSubopt'
            if options.del_worse:
                outfile_base += '_dW'
            if options.intron_features != '-':
                outfile_base += '_featFile'
            
            if options.bam_input:
                outfile_base = (re.sub('.bam', '', outfile_base) + '.bam')
                outfile_handle = subprocess.Popen([options.samtools, 'view', '-bS', '-o' + outfile_base, '-'], stdin=subprocess.PIPE)
                outfile = outfile_handle.stdin
            else:
                outfile_base = (re.sub('.sam', '', outfile_base) + '.sam')
                outfile = open(outfile_base, 'w')
        else:
            outfile = sys.stdout
    else:
        if options.bam_input:
            outfile_handle = subprocess.Popen([options.samtools, 'view', '-bS', '-o' + options.outfile, '-'], stdin=subprocess.PIPE)
            outfile = outfile_handle.stdin
        else:
            outfile = open(options.outfile, 'w')

    first_line = True
    curr_lines = []

    options.features = dict()
    if options.intron_features != '-':
        line_counter = 0
        for line in open(options.intron_features, 'r'):
            if options.verbose and line_counter % 10000 == 0:
                print 'parsed %i features from %s' % (line_counter, options.intron_features)
            line_counter += 1
            sl = line.strip().split('\t')
            (chrm, start, stop) = sl[:3]
            options.features[(chrm, start, stop)] = Feature(80, sl[3:])

    if options.align != '-':
        if options.bam_input:
            infile_handle = subprocess.Popen([options.samtools, 'view', '-h', options.align], stdout=subprocess.PIPE)
            infile = infile_handle.stdout
        else:
            infile = open(options.align, 'r')
    else:
        infile = sys.stdin

    for line in infile:
        if counter % 10000 == 0 and options.verbose:
            print >> sys.stderr, 'lines read: [ %s (taken: %s / filtered: %s)]' % (counter, counter - filter_counter, filter_counter)
        counter += 1

        sl = line.strip().split('\t')
        if len(sl) < 9:
            sl = line.strip().split(' ')
        if len(sl) < 9:
            print >> outfile, line, 
            continue
            
        sl_id = (sl[0], int(sl[1]) & 196 / 64)

        curr_lines.append(sl)

        if options.del_worse or options.no_suboptimal:
            if first_line:
                first_line = False
                last_id = sl_id
            if sl_id == last_id:
                continue
            else:
                assert (len(curr_lines) >= 2)
                curr_lines = curr_lines[:-1]
                last_id = sl_id

        if options.no_suboptimal:
            (curr_lines, filter_counter) = filter_suboptimal(curr_lines, filter_counter, options)

        (min_mm, curr_lines, filter_counter) = filter_lines_by_criteria(curr_lines, multireads, options, filter_counter)

        if options.del_worse:
            (curr_lines, filter_counter) = filter_lines_by_min_mm(curr_lines, min_mm, filter_counter)
        
        for _sl in curr_lines:
            _line = ''
            for bb in _sl:
                _line += (bb + '\t')
            print >> outfile, _line[:-1]

        if options.del_worse or options.no_suboptimal:
            curr_lines = [sl]
        else:
            curr_lines = []

    if len(curr_lines) > 0:
        if options.no_suboptimal:
            (curr_lines, filter_counter) = filter_suboptimal(curr_lines, filter_counter, options)

        (min_mm, curr_lines, filter_counter) = filter_lines_by_criteria(curr_lines, multireads, options, filter_counter)

        if options.del_worse:
            (curr_lines, filter_counter) = filter_lines_by_min_mm(curr_lines, min_mm, filter_counter)

        for sl in curr_lines:
            _line = ''
            for bb in sl:
                _line += (bb + '\t')
            print >> outfile, _line[:-1]

    if options.align != '-':
        infile.close()
    outfile.close()

if __name__ == '__main__':
    main()
