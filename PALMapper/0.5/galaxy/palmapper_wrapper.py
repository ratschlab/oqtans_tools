
"""
Runs Palmapper on single-end or paired-end data.
"""

import optparse, os, sys, tempfile, shutil, time, re

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()
 
def __main__():
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
    parser.add_option('', '--outsuffix', dest='outsuffix', default='', help='Append this suffix to all output file names') 
    
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

    # append output suffix
    options.unspliced_output = options.unspliced_output + options.outsuffix ;
    options.spliced_output = options.spliced_output + options.outsuffix ;
    options.sam_output = options.sam_output + options.outsuffix ;
    options.bam_output = options.bam_output + options.outsuffix ;
    options.junctions_output = options.junctions_output + options.outsuffix ;
    options.coverage_output = options.coverage_output + options.outsuffix ;
    options.logfile = options.logfile + options.outsuffix ;
    
    # index if necessary
    if options.genomeSource == 'history':
        
        if options.indexSource == 'array':
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
            cmd1 = 'pmindex -v -i %s %s' % (options.ref, indexing_cmds)
            try:
                os.system(cmd1)
            except Exception, erf:
                stop_err('Error indexing reference sequence\n' + str(erf))
        else:
            # make temp directory for placement of indices and copy reference file there
            tmp_dir = tempfile.gettempdir()
            try:
                os.system('cp %s %s' % (options.ref, tmp_dir))
            except Exception, erf:
                stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
                options.ref = os.path.join(tmp_dir, os.path.split(options.ref)[1])
            cmd1 = 'bwa index %s' % (options.ref)
            try:
                os.system(cmd1)
            except Exception, erf:
                stop_err('Error indexing reference sequence\n' + str(erf))

    #GenomeMapper parameters
    if options.params == 'pre_set':
        # Auto values for: -M -G -E -z
        # Supporting only one thread
        aligning_cmds = '-l 18 -seed-hit-cancel-threshold 10000 '
    else:
        try:
            aligning_cmds = '%s %s %s %s %s %s %s ' % \
                            (('','-l %s' % options.alignseedlength)[options.alignseedlength!='None'],
                             ('','-M %s' % options.maxmismatches)[options.maxmismatches!='None'],
                             ('','-G %s' % options.maxgaps)[options.maxgaps!='None'],
                             ('','-E %s' % options.maxedits)[options.maxedits!='None'],
                             ('','-seed-hit-cancel-threshold %s' % options.seedhitcancel)[options.seedhitcancel!='None'],
                             #('','-threads %s' % options.threads)[options.threads!='None'],
                             ('','-z %s' % options.topalignment)[options.topalignment!='None'],
                             ('','-a')[options.reportall!='false'])
        except ValueError, erf:
            stop_err('Something is wrong with the alignment parameters and the alignment could not be run\n' + str(erf))
    

    #Index type
    aligning_cmds+=('','-bwa 12 ')[options.indexSource=="bwa"]

    #QPALMA parameters
    if options.qpalma_params == 'pre_set':
        # Auto values: -L -K -C -I -NI -QMM
        qpalma_cmds = '-filter-max-mismatches 1 -filter-max-gaps 0 -SA 10 -CT 10 -qpalma-use-map-max-len 5000 -report-splice-sites 0.9 -report-map-read -report-spliced-read -report-map-region -S '
    else:
        try:
            #print options
            qpalma_cmds = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s -S ' % \
                          (('','-filter-max-mismatches %s' % options.mmtrigger)[options.mmtrigger!='None'],
                           ('','-filter-max-gaps %s' % options.gtrigger)[options.gtrigger!='None'],
                           ('','-SA %s' % options.maxalign)[options.maxalign!='None'],
                           ('','-L %s' % options.aligntrigger)[options.aligntrigger!='None'],
                           ('','-K %s' % options.alignshorttrigger)[options.alignshorttrigger!='None'],
                           ('','-C %s' % options.aligncombinedtrigger)[options.aligncombinedtrigger!='None'],
                           ('','-I %s' % options.maxintronlength)[options.maxintronlength!='None'],
                           ('','-NI %s' % options.maxintronnum)[options.maxintronnum!='None'],
                           ('','-QMM %s' % options.qmm)[options.qmm!='None'],
                           ('','-CT %s' % options.clustertol)[options.clustertol!='None'],
                           ('','-qpalma-use-map-max-len %s' % options.mapmaxlen)[options.mapmaxlen!='None'],
                           ('','-report-splice-sites %s' % options.report_ss)[options.report_ss!='None'],
                           ('','-rtrim %s ' % options.rtrim)[options.rtrim!='None'],
                           ('','-rtrim-step %s ' % options.rtrim_step)[options.rtrim!='None'],
                           ('','-polytrim %s ' % options.polytrim)[options.polytrim!='None'],
                           ('','-junction-remapping %s ' % options.junction_remapping)[options.junction_remapping != 'None'],
                           ('','-junction-remapping-coverage %s ' % options.junction_coverage)[options.junction_remapping != 'None'])

            qpalma_cmds +=('','-report-spliced-read ')[options.reportsplicedread=='true']
            qpalma_cmds +=('','-report-map-read ')[options.reportmappedread=='true']
            qpalma_cmds +=('','-non-consensus-search ')[options.non_consensus=='true']

        except ValueError, erf:
            stop_err('Something is wrong with the QPALMA alignment parameters and the alignment could not be run\n' + str(erf))


    # creating copy of critical files on local tmp file system
    # Reference genome
    index_tmp_dname = tempfile.mkdtemp(suffix='', prefix='gmindex_tmp_', dir=None) ;
    if options.ref[0:5]!='/tmp/':
        try:
            os.system('cp %s* %s' % (options.ref, index_tmp_dname))
        except Exception, erf:
            stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
        options.ref = os.path.join(index_tmp_dname, os.path.split(options.ref)[1])

    #Splice site predictions
    if (options.ss_pred == "true"):
        acc_tmp_dname = tempfile.mkdtemp(suffix='', prefix='acc_', dir=None) ;
        if os.path.isdir(os.path.join(options.acc,'pred')):
            try:
                os.system('cp %s/pred/contig_* %s' % (options.acc, acc_tmp_dname))
            except Exception, erf:
                stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
        else:
            try:
                os.system('cp %s/contig_* %s' % (options.acc, acc_tmp_dname))
            except Exception, erf:
                stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
    
        options.acc = os.path.join(acc_tmp_dname, 'contig_%i%c')

        don_tmp_dname = tempfile.mkdtemp(suffix='', prefix='don_', dir=None) ;
        if os.path.isdir(os.path.join(options.don,'pred')):
            try:
                os.system('cp %s/pred/contig_* %s' % (options.don, don_tmp_dname))
            except Exception, erf:
                stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
        else:
            try:
                os.system('cp %s/contig_* %s' % (options.don, don_tmp_dname))
            except Exception, erf:
                stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
        options.don = os.path.join(don_tmp_dname, 'contig_%i%c') 
        
        ss_cmds = '-acc %s -don %s ' % (options.acc, options.don)

    else:
        ss_cmds = '-no-ss-pred '

    # prepare actual aligning commands
    (report_file, report_fname) = tempfile.mkstemp(suffix='', prefix='report_', dir=None) 
    os.close(report_file)
    try:
        os.unlink(report_fname)
    except:
        pass
    

    ## Output files
    if options.format == 'sam':
        output_cmd='-f sam -o %s ' % options.sam_output
        if options.unmapped_included == 'true':
            output_cmd+='-include-unmapped-reads '
        else:
            (unmapped_tmp_file, unmapped_tmp_fname) = tempfile.mkstemp(suffix='', prefix='unmapped_', dir=None) ;
            os.close(unmapped_tmp_file) ;
            output_cmd+='-u %s ' % unmapped_tmp_file
    elif options.format == 'bam':
        if options.bamsorting == "position":
            output_cmd='-f bamp '
        elif options.bamsorting == "read":
            output_cmd='-f bamn '
        else:
            output_cmd='-f bam '
        output_cmd+='-o %s ' % options.bam_output
        if options.unmapped_included == 'true':
            output_cmd+='-include-unmapped-reads '
        else:
            (unmapped_tmp_file, unmapped_tmp_fname) = tempfile.mkstemp(suffix='', prefix='unmapped_', dir=None) ;
            os.close(unmapped_tmp_file) ;
            output_cmd+='-u %s ' % unmapped_tmp_file
    else: #bedx output
        (unmapped_tmp_file, unmapped_tmp_fname) = tempfile.mkstemp(suffix='', prefix='unmapped_', dir=None) ;
        os.close(unmapped_tmp_file) ;
        output_cmd='-f bedx -o %s -H %s -u %s ' % (options.unspliced_output,options.spliced_output, unmapped_tmp_file)

    if options.coverage == 'true':
        output_cmd+='-report-coverage-wig %s ' % options.coverage_output
    if options.junctions == 'true':
        output_cmd+='-report-junctions %s ' % options.junctions_output

    ## Input files
    if options.paired == 'paired':
        input_cmd='-q1 %s -q2 %s ' % (options.input1, options.input2)
        if options.protocol != 'unstranded':
            input_cmd+='-protocol %s ' % options.protocol
    else:
        assert( options.paired == 'single' )
        input_cmd='-q %s ' % options.input1
        if options.strand != 'unstranded':
            input_cmd+='-strand %s ' % options.strand
            if options.protocol != 'unstranded':
                input_cmd+='-protocol %s ' % options.protocol

        
    cmd2a = 'palmapper %s %s -i %s %s  %s -qpalma %s %s -report %s -threads 4 -qpalma-prb-offset-fix ' % (aligning_cmds, qpalma_cmds, options.ref, input_cmd, output_cmd, options.qpalma, ss_cmds, report_fname) #, options.logfile)
    

    # align
    try:
        #os.environ['LD_LIBRARY_PATH']='/home/galaxy/svn/projects/QPalma/dyn_prog/cpplib/:/home/galaxy/software/shogun/lib/'
#        print re.sub(r'palmapper', r'GALAXY-SOFTWARE-DIR', cmd2a)
        #print re.sub(r'/home/galaxy/software/palmapper-0.4/palmapper', r'GALAXY-SOFTWARE-DIR', cmd2a)
        print cmd2a
        os.system(cmd2a)
    except Exception, erf:
        stop_err("Error aligning sequence\n" + str(erf))

    try:
        shutil.rmtree(index_tmp_dname)
        shutil.rmtree(acc_tmp_dname)
        shutil.rmtree(don_tmp_dname)
        os.unlink(report_fname)
    except:
        pass
	
    etime = time.asctime( time.localtime(time.time()) )
    print '----------------------------------------------'
    print 'PALMapper finished on ' + etime
    print '----------------------------------------------' 

if __name__=="__main__": __main__()
