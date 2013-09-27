
"""
Runs Palmapper on single-end or paired-end data.
"""

import optparse, os, sys, tempfile, shutil, time, re

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()
 
def __main__():
    os.environ['PATH']=os.environ['PATH']+":"+os.environ['OQTANS_PATH']+ "/palmapper/0.5/"
    stime = time.asctime( time.localtime(time.time()) )
    print '----------------------------------------------'
    print 'PALMapper started on ' + stime
    print '----------------------------------------------' 
    #Parse Command Line
    parser = optparse.OptionParser()

    parser.add_option('', '--logfile', dest='logfile', help='log file')

    #Read files
    parser.add_option('', '--paired', dest='paired', help='Whether the data is single- or paired-end')
    parser.add_option('', '--input1', dest='input1', help='The (forward or single-end) reads file in Sanger FASTQ format')
    parser.add_option('', '--input2', dest='input2', help='The reverse reads file in Sanger FASTQ format')
    parser.add_option('', '--strand', dest='strand', help='Strand information (left or right)')
    parser.add_option('', '--protocol', dest='protocol', help='Protocol used (first or second)')

    #Reference genome and index information
    parser.add_option('', '--indexSource', dest='indexSource', default='array', help='The type of index: bwa or array')
    parser.add_option('', '--genomeSource', dest='genomeSource', help='The type of reference provided')
    parser.add_option('', '--ref', dest='ref', help='The reference genome to use or index')
    parser.add_option('', '--indexSettings', dest='index_settings', help='Whether or not indexing options are to be set')
    parser.add_option('', '--seedlength', dest='seedlength', help='Index Seed Length')

    # Splice site predictions
    parser.add_option('', '--ss-pred', dest='ss_pred', help='use splice site predictions')
    parser.add_option('', '--acc', dest='acc', help='Acceptor SS predictions')
    parser.add_option('', '--don', dest='don', help='Donor SS predictions')

    #Output files
    parser.add_option('', '--format', dest='format', help='Output format (bedx, sam or bam)')
    parser.add_option('', '--unspliced-output', dest='unspliced_output', help='The Bedx output file for unspliced reads')
    parser.add_option('', '--spliced-output', dest='spliced_output', help='The Bedx output file for spliced reads')
    parser.add_option('', '--sam-output', dest='sam_output', help='The SAM output file for both spliced and unspliced reads')
    parser.add_option('', '--bam-output', dest='bam_output', help='The BAM output file for both spliced and unspliced reads')
    parser.add_option('', '--bamsort', dest='bamsorting', help='Type of sorting for BAM output (unsorted, position or read)')
    parser.add_option('', '--include-unmapped', dest='unmapped_included', help='Whether unmapped reads are included in output file (only for SAM and BAM format)')

    parser.add_option('', '--coverage-map', dest='coverage', help='Whether the coverage map should be output')
    parser.add_option('', '--junctions', dest='junctions', help='Whether the intron junction library should be built')
    parser.add_option('', '--coverage-output', dest='coverage_output', help='Coverage map output')
    parser.add_option('', '--junctions-output', dest='junctions_output', help='Intron junctions file')

    #GenomeMapper parameters
    parser.add_option('', '--params', dest='params', help='Whether to use default or specified parameters for GenomeMapper')
    parser.add_option('', '--alignseedlength', dest='alignseedlength', help='Alignment Seed Length')
    parser.add_option('', '--maxmismatches', dest='maxmismatches', help='Maximal number of mismatches')
    parser.add_option('', '--maxgaps', dest='maxgaps', help='Maximal number of gaps')
    parser.add_option('', '--maxedits', dest='maxedits', help='Maximal number of edit operations')
    parser.add_option('', '--seedhitcancel', dest='seedhitcancel', help='Number of hits of a seed that lead to its ignoration')
    parser.add_option('', '--threads', dest='threads', help='The number of threads to run')
    parser.add_option('', '--topalignment', dest='topalignment', help='Number of top alignments to report')
    parser.add_option('', '--reportall', dest='reportall', help='Report all alignments')

    #QPALMA parameters
    parser.add_option('', '--qpalma', dest='qpalma', help='QPALMA parameter file')
    parser.add_option('', '--qpalma-params', dest='qpalma_params', help='Whether to use default or specified parameters for QPALMA')
    parser.add_option('', '--mmtrigger', dest='mmtrigger', help='Mismatch threshold to trigger spliced alignments')
    parser.add_option('', '--gtrigger', dest='gtrigger', help='Gap threshold to trigger spliced alignments')
    parser.add_option('', '--maxalign', dest='maxalign', help='Maximal number of spliced alignments per read')
    parser.add_option('', '--aligntrigger', dest='aligntrigger', help='Minimal length of long hit')
    parser.add_option('', '--alignshorttrigger', dest='alignshorttrigger', help='Minimal length of short hit')
    parser.add_option('', '--aligncombinedtrigger', dest='aligncombinedtrigger', help='Minimal combined length')
    parser.add_option('', '--maxintronlength', dest='maxintronlength', help='Maximal intron length')
    parser.add_option('', '--maxintronnum', dest='maxintronnum', help='Maximal number of introns')
    parser.add_option('', '--qmm', dest='qmm', help='Number of matches required for identifying a splice site')
    parser.add_option('', '--clustertol', dest='clustertol', help='Distance in nt to tolerate between hit and existing hit cluster')
    parser.add_option('', '--qpalma-use-map-max-len', dest='mapmaxlen', help='Up and downstream limit of map extension')
    #parser.add_option('', '--filter_ss_tophit', dest='filter_ss_tophit', help='filter_ss_tophit')
    parser.add_option('', '--report_ss', dest='report_ss', help='Splice site-based alignment regions')
    parser.add_option('', '--reportmappedread', dest='reportmappedread', help='Use mapped unspliced reads for determining alignment regions')
    parser.add_option('', '--reportsplicedread', dest='reportsplicedread', help='Use mapped spliced reads for determining alignment regions')

    parser.add_option('', '--rtrim', dest='rtrim', help='Minimal length of read when trimming the righ side')
    parser.add_option('', '--rtrim-step', dest='rtrim_step', help='Right trimming step')
    parser.add_option('', '--polytrim', dest='polytrim', help='Minimal length of read when trimming polyA or polyT ends')

    parser.add_option('', '--junction-remapping', dest='junction_remapping', help='Intron junctions file for remapping strategy (Gff3 format)')
    parser.add_option('', '--junction-coverage', dest='junction_coverage', help='Minimal intron junction support for remapping strategy')
    parser.add_option('', '--non-consensus-search', dest='non_consensus', help='Whether spliced alignments with non consensus sequences as splice sites are searched')

    (options, args) = parser.parse_args()


if __name__=="__main__": __main__()
