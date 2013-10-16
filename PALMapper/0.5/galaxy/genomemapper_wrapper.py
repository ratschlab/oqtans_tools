#! /usr/bin/python

"""
Runs GenomeMapper on single-end or paired-end data.
"""

import optparse, os, sys, tempfile

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()
 
def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('', '--threads', dest='threads', help='The number of threads to run')
    parser.add_option('', '--input1', dest='input1', help='The (forward or single-end) reads file in Sanger FASTQ format')
    parser.add_option('', '--input2', dest='input2', help='The reverse reads file in Sanger FASTQ format')
    parser.add_option('', '--output', dest='output', help='The output file')
    parser.add_option('', '--paired', dest='paired', help='Whether the data is single- or paired-end')
    parser.add_option('', '--genomeSource', dest='genomeSource', help='The type of reference provided')
    parser.add_option('', '--ref', dest='ref', help='The reference genome to use or index')
    parser.add_option('', '--indexSettings', dest='index_settings', help='Whether or not indexing options are to be set')
    parser.add_option('', '--params', dest='params', help='Whether to use default or specified parameters')
    parser.add_option('', '--seedlength', dest='seedlength', help='GenomeMapper Index Seed Length')
    parser.add_option('', '--alignseedlength', dest='alignseedlength', help='GenomeMapper Alignment Seed Length')
    parser.add_option('', '--format', dest='format', help='Output format (bed or shore)')
    parser.add_option('', '--maxmismatches', dest='maxmismatches', help='Maximal number of mismatches')
    parser.add_option('', '--maxgaps', dest='maxgaps', help='Maximal number of gaps')
    parser.add_option('', '--maxedits', dest='maxedits', help='Maximal number of edit operations')
    parser.add_option('', '--reportall', dest='reportall', help='Report all hits')
    (options, args) = parser.parse_args()
    
    # index if necessary
    if options.genomeSource == 'history':
        # set up commands
        if options.index_settings =='index_pre_set':
            indexing_cmds = ''
        else:
            try:
                indexing_cmds = '%s ' % \
                                (('','-s %s'%options.seedlength)[options.seedlength!='None' and options.seedlength>=1])
            except ValueError:
                indexing_cmds = ''
                
        # make temp directory for placement of indices and copy reference file there
        tmp_dir = tempfile.gettempdir()
        try:
            os.system('cp %s %s' % (options.ref, tmp_dir))
        except Exception, erf:
            stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
        options.ref = os.path.join(tmp_dir, os.path.split(options.ref)[1])
        cmd1 = 'gmindex -v -i %s %s' % (options.ref, indexing_cmds)
        try:
            os.system(cmd1)
        except Exception, erf:
            stop_err('Error indexing reference sequence\n' + str(erf))
    
    if options.params == 'pre_set':
        aligning_cmds = '-v '
    else:
        try:
            print options
            aligning_cmds = '%s %s %s %s %s %s -v ' % \
                            (('','-f %s' % options.format)[options.format!='None'],
                             ('','-a')[options.reportall!='None'],
                             ('','-M %s' % options.maxmismatches)[options.maxmismatches!='None'],
                             ('','-G %s' % options.maxgaps)[options.maxgaps!='None'],
                             ('','-E %s' % options.maxedits)[options.maxedits!='None'],
                             ('','-l %s' % options.alignseedlength)[options.alignseedlength!='None'])
        except ValueError, erf:
            stop_err('Something is wrong with the alignment parameters and the alignment could not be run\n' + str(erf))

    # prepare actual aligning commands
    if options.paired == 'paired':
        print "Sorry, paired end alignments are not implemented yet"
        return
        #cmd2 = 'genomemapper %s %s -1 %s -2 %s > %s ' % (options.ref, options.input1, options.input2, options.output) 
    else:
        cmd2 = 'genomemapper %s -i %s -q %s -o %s ' % (aligning_cmds, options.ref, options.input1, options.output) 

    # align
    try:
        print cmd2
        os.system(cmd2)
    except Exception, erf:
        stop_err("Error aligning sequence\n" + str(erf))

if __name__=="__main__": __main__()
