"""
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2095-2010 Andre Kahles
  Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

  This script finds an optimal parameter set to maximize the performance of a 
  given intronfeature file.
  
  For detailed usage information type:

    python find_optimal_param_set.py 

"""

import sys
import cPickle

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
        """Returns submission support"""

        return int(self.submission_support)
    
    def is_valid(self, mm, ex, mc, options):
        """Returns true, if at least one alignment fulfills the requirements with respect to mm, ex, and mc. False otherwise."""

        if self.alignment_support < mc:
            return False

        is_valid = False
        for key in self.mm_ex.keys():
            _ex = int(key) / (options.max_feat_mismatches + 1)
            _mm = int(key) % (options.max_feat_mismatches + 1)
            if _mm <= mm and _ex >= ex:
                is_valid = True
                break
        
        return is_valid


def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-b', '--best_score', dest='best_scores', metavar='FILE', help='file to store the best scoring parameters', default='-')
    required.add_option('-m', '--matrix', dest='matrix', metavar='FILE', help='file to store the full performance matrix', default='-')
    required.add_option('-f', '--features', dest='features', metavar='FILE', help='alignment intron features', default='-')
    required.add_option('-i', '--annotation_introns', dest='anno_int', metavar='FILE', help='annotation intron list', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-E', '--exclude_introns', dest='exclude_introns', metavar='STRINGLIST', help='list of comma separated intron files to exclude from submitted features', default='-')
    optional.add_option('-I', '--max_intron_len', dest='max_intron_len', metavar='INT', type='int', help='maximal intron length [10000000]', default=10000000)
    optional.add_option('-s', '--ignore_strand', dest='ignore_strand', action='store_true', help='ignore strand information present in annotation', default=False)
    optional.add_option('-X', '--max_feat_mismatches', dest='max_feat_mismatches', metavar='INT', type='int', help='max number of mismatches for feat generation [80] (do only change, if you are absolutely sure!)', default=80)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def get_performance_value(full_features, mm, ex, mc, annotation_list, options):
    """Builds up a filtered intron list from the given alignment features and compares to the annotation."""

    alignment_list = dict()

    for feat in full_features.keys():

        chrm = feat[0]
        intron = (0, int(feat[1]), int(feat[2]))

        ### filter step
        if (intron[2] - intron[1]) > options.max_intron_len:
            continue
        if not full_features[feat].is_valid(mm, ex, mc, options):
            continue

        try:
            alignment_list[chrm][intron] = 0
        except KeyError:
            alignment_list[chrm] = {intron:0}

    ### match intron lists
    total_precision = float(0)
    total_recall = float(0)
    key_count = 0
    for chrm in annotation_list.keys():
        if alignment_list.has_key(chrm):
            matches = len(set(annotation_list[chrm].keys()).intersection(set(alignment_list[chrm].keys())))
            total_precision += (float(matches) / float(max(1, len(alignment_list[chrm].keys()))))
            total_recall += (float(matches) / float(max(1, len(annotation_list[chrm].keys()))))
            ### do not include chromosomes with zero values into average
            if matches > 0:
                key_count += 1
    total_precision /= max(1.0, float(key_count))
    total_recall /= max(1.0, float(key_count))

    return (total_precision, total_recall)


def main():
    """Main function extracting intron features."""

    options = parse_options(sys.argv)

    ### get list of annotated introns
    annotation_list = cPickle.load(open(options.anno_int, 'r'))

    if options.ignore_strand:
        for chrm in annotation_list.keys():
            skiplist = set()
            for intron in annotation_list[chrm].keys():
                if intron[0] == 0:
                    continue
                annotation_list[chrm][(0, intron[1], intron[2])] = annotation_list[chrm][intron]
                skiplist.add(intron)
            for intron in skiplist:
                del annotation_list[chrm][intron]
            del skiplist

    ### filter annotation for max intron length
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
    del skiplist


    full_features = dict() 

    if options.verbose:
        print 'Parsing %s' % options.features

    line_counter = 0
    for line in open(options.features, 'r'):
        if options.verbose and line_counter % 1000 == 0:
            print 'parsed %i features from %s' % (line_counter, options.features)
        line_counter += 1
        sl = line.strip().split('\t')

        (chrm, start, stop) = sl[:3]
        try:
            full_features[(chrm, start, stop)].full_features(sl[3:])
        except KeyError:
            full_features[(chrm, start, stop)] = Feature(80, sl[3:])


    ### filter full feature list for excluded introns
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
                for _intron in _ex_intron[chrm].keys():
                    try:
                        del full_features[(chrm, str(_intron[1]), str(_intron[2]))]
                    except KeyError:
                        continue
            del _ex_intron

    if options.verbose:
        print 'Parsing completed.' 
        print 'parsed %i features from %s' % (line_counter, options.features)

    ### SEARCH SPACE
    ### iterate over different filter dimensions
    #ex_list = [2, 4, 6, 8, 10, 12, 15, 20, 25, 30] # 10
    ex_list = [2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 ] # 15
    mm_list = [0, 1, 2, 3, 4, 5, 6]                # 7
    mc_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]      # 10   ==> 700 combinations

    checked_combs = 0
    # pre rec fsc
    max_pre = (0.0, 0.0, 0.0)
    max_rec = (0.0, 0.0, 0.0)
    max_fsc = (0.0, 0.0, 0.0)
    max_pre_idx = (0, 0, 0)
    max_rec_idx = (0, 0, 0)
    max_fsc_idx = (0, 0, 0)

    matrix_file = open(options.matrix, 'w')
    for ex in ex_list:
        for mm in mm_list:
            for mc in mc_list:
                if options.verbose and checked_combs % 10 == 0:
                    print 'checked %i parameter combinations' % checked_combs
                    print 'best scores so far:\n \tbest fScore: %0.2f, best recall: %0.2f, best precision: %0.2f' % (max_fsc[2], max_rec[1], max_pre[0])
                checked_combs += 1

                (pre, rec) = get_performance_value(full_features, mm, ex, mc, annotation_list, options)    

                if float(rec) + float(pre) > 0:
                    fsc = (2 * float(rec) * float(pre)) / (float(rec) + float(pre)) 
                else:
                    fsc = 0.0

                if pre > max_pre[0]:
                    max_pre = (pre, rec, fsc)
                    max_pre_idx = (ex, mm, mc)
                if rec > max_rec[1]:
                    max_rec = (pre, rec, fsc)
                    max_rec_idx = (ex, mm, mc)
                if fsc > max_fsc[2]:
                    max_fsc = (pre, rec, fsc)
                    max_fsc_idx = (ex, mm, mc)
                
                ### store information
                ### ex mm mc pre rec fsc
                print >> matrix_file, '%s\t%s\t%s\t%s\t%s\t%s' % (ex, mm, mc, pre, rec, fsc)

    matrix_file.close()

    best_file = open(options.best_scores, 'w')
    # best precision 
    print >> best_file, '%s\t%s\t%s\t%s\t%s\t%s' % (max_pre_idx[0], max_pre_idx[1], max_pre_idx[2], max_pre[0], max_pre[1], max_pre[2])
    # best recall 
    print >> best_file, '%s\t%s\t%s\t%s\t%s\t%s' % (max_rec_idx[0], max_rec_idx[1], max_rec_idx[2], max_rec[0], max_rec[1], max_rec[2])
    # best fScore 
    print >> best_file, '%s\t%s\t%s\t%s\t%s\t%s' % (max_fsc_idx[0], max_fsc_idx[1], max_fsc_idx[2], max_fsc[0], max_fsc[1], max_fsc[2])
    best_file.close()

   
if __name__ == "__main__":
    main()

