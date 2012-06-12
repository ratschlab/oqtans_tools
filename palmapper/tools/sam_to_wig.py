#!/usr/bin/env python
# Description: Convert SAM files to Wiggle variableStep. 
# Usage: sam_to_wig.py --input=<SAM_file> --ref_file=<reference_genome_file> --output=<wiggle_ouptut_path> ----expName=<experiment_name>
# samtools required.
import optparse, os, sys, subprocess, tempfile, shutil, gzip, re

def __main__():
    ############################################# 1 ########################################################
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--input', dest='input', help='The input SAM dataset' )
    parser.add_option( '', '--ref_file', dest='ref_file', help='The reference dataset' )
    parser.add_option( '', '--output', dest='output', help='The output Wiggle dataset' )
    parser.add_option( '', '--expName', dest='expName', help='Experiment name' )
    ( options, args ) = parser.parse_args()
    tmp_dir = tempfile.mkdtemp()
    if options.ref_file == 'None':
        sys.stderr.write( 'No reference sequence available. Cannot continue, terminating ....\n' )
        sys.exit(-1)
    else:
        try:
            # Create indexes for history reference ( e.g., ~/database/files/000/dataset_1.dat ) using samtools faidx, which will:
            # - index reference sequence in the FASTA format or extract subsequence from indexed reference sequence
            # - if no region is specified, faidx will index the file and create <ref.fasta>.fai on the disk
            # - if regions are specified, the subsequences will be retrieved and printed to stdout in the FASTA format
            # - the input file can be compressed in the RAZF format.
            fai_index_file_base = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
            # At this point, fai_index_file_path will look something like /tmp/dataset_13.dat
            os.symlink( options.ref_file, fai_index_file_base )
            fai_index_file_path = '%s.fai' % fai_index_file_base
            command = 'samtools faidx %s' % fai_index_file_base
            tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
            tmp_stderr = open( tmp, 'wb' )
            proc = subprocess.Popen( args=command, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
            returncode = proc.wait()
            tmp_stderr.close()
            # get stderr, allowing for case where it's very large
            tmp_stderr = open( tmp, 'rb' )
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read( buffsize )
                    if not stderr or len( stderr ) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stderr.close()
            if returncode != 0:
                raise Exception, stderr 
            if len( open( fai_index_file_path ).read().strip() ) == 0:
                raise Exception, 'Index file empty, there may be an error with your reference file or settings.'
        except Exception, e:
            #clean up temp files
            if os.path.exists( tmp_dir ):
                shutil.rmtree( tmp_dir )
            sys.stderr.write( 'Error creating indexes from reference (%s), %s\n' % ( options.ref_file, str( e ) ) )
    try:
        # Extract all alignments from the input SAM file to BAM format ( since no region is specified, all the alignments will be extracted ).
        tmp_aligns_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
        tmp_aligns_file_name = tmp_aligns_file.name
        tmp_aligns_file.close()
        
        # IMPORTANT NOTE: for some reason the samtools view command gzips the resulting bam file without warning,
        # and the docs do not currently state that this occurs ( very bad ).
        command = 'samtools view -bt %s -o %s %s' % ( fai_index_file_path, tmp_aligns_file_name, options.input )
        tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
        tmp_stderr = open( tmp, 'wb' )
        proc = subprocess.Popen( args=command, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()
         
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr
        if len( open( tmp_aligns_file_name ).read() ) == 0:
            raise Exception, 'Initial BAM file empty'
    except Exception, e:
        #clean up temp files
        if os.path.exists( tmp_dir ):
            shutil.rmtree( tmp_dir )
        sys.stderr.write( 'Error extracting alignments from (%s), %s' % ( options.input1, str( e ) ) )
    try:
        # Sort alignments by leftmost coordinates. File <out.prefix>.bam will be created. This command
        # may also create temporary files <out.prefix>.%d.bam when the whole alignment cannot be fitted
        # into memory ( controlled by option -m ).
        tmp_sorted_aligns_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
        tmp_sorted_aligns_file_name = tmp_sorted_aligns_file.name
        tmp_sorted_aligns_file.close()
        command = 'samtools sort %s %s' % ( tmp_aligns_file_name, tmp_sorted_aligns_file_name )
        tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
        tmp_stderr = open( tmp, 'wb' )
        proc = subprocess.Popen( args=command, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        #clean up temp files
        if os.path.exists( tmp_dir ):
            shutil.rmtree( tmp_dir )
        sys.stderr.write( 'Error sorting alignments from (%s), %s' % ( tmp_aligns_file_name, str( e ) ) )
    # Move tmp_aligns_file_name to our output dataset location
    sorted_bam_file = '%s.bam' % tmp_sorted_aligns_file_name
    if os.path.getsize( sorted_bam_file ) == 0:
        #clean up temp files
        if os.path.exists( tmp_dir ):
            shutil.rmtree( tmp_dir )
        sys.stderr.write( 'Error creating sorted version of BAM file' )
    shutil.move( sorted_bam_file, options.output )
    #clean up temp files
    if os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )
    # check that there are results in the output file
    if os.path.getsize( options.output ) > 0:
        sys.stdout.write( 'SAM file converted to BAM\n' )
    else:
        sys.stdout.write( 'The output file is empty, there may be an error with your input file or settings.' )
    ############################################# 2 ########################################################
    tmpDir = tempfile.mkdtemp()
    tmpf0 = tempfile.NamedTemporaryFile( dir=tmpDir )
    tmpf0_name = tmpf0.name
    tmpf0.close()
    tmpf0bam_name = '%s.bam' % tmpf0_name
    tmpf0bambai_name = '%s.bam.bai' % tmpf0_name
    tmpf1 = tempfile.NamedTemporaryFile( dir=tmpDir )
    tmpf1_name = tmpf1.name
    tmpf1.close()
    tmpf1fai_name = '%s.fai' % tmpf1_name
    ##link bam and bam index to working directory (can't move because need to leave original)
    in_bamfile = os.listdir(options.output)
    inb_file = ''
    for s in in_bamfile:
         if re.search(r'.+\.bam$', s):
             inb_file = s 
             break
    os.symlink( inb_file, tmpf0bam_name )
    #os.symlink( options.bamIndex, tmpf0bambai_name )
    opts = '-M %s' % ( 60 ) # --lastCol=no, --indels=no, --mapCap=60
    cmd = 'samtools pileup %s -f %s %s > %s'
    try:
        # have to nest try-except in try-finally to handle 2.4
        try:
            os.symlink( options.ref_file, tmpf1_name )
            cmdIndex = 'samtools faidx %s' % ( tmpf1_name )
            tmp = tempfile.NamedTemporaryFile( dir=tmpDir ).name
            tmp_stderr = open( tmp, 'wb' )
            proc = subprocess.Popen( args=cmdIndex, shell=True, cwd=tmpDir, stderr=tmp_stderr.fileno() )
            returncode = proc.wait()
            tmp_stderr.close()
            # get stderr, allowing for case where it's very large
            tmp_stderr = open( tmp, 'rb' )
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read( buffsize )
                    if not stderr or len( stderr ) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stderr.close()
            #did index succeed?
            if returncode != 0:
                raise Exception, 'Error creating index file\n' + stderr
            txt_file = tempfile.NamedTemporaryFile( dir=tmpDir ).name
	        f_name = re.search(r'^\/tmp\/.+\/(.+)$', txt_file)
            f_name = f_name.group(1)
            pileup_file = options.output+'/'+f_name+'.txt'
            cmd = cmd % ( opts, options.ref_file, options.output+'/'+inb_file, pileup_file )
            #perform pileup command
            tmp = tempfile.NamedTemporaryFile( dir=tmpDir ).name
            tmp_stderr = open( tmp, 'wb' )
            proc = subprocess.Popen( args=cmd, shell=True, cwd=tmpDir, stderr=tmp_stderr.fileno() )
            returncode = proc.wait()
            tmp_stderr.close()
            #did it succeed?
            # get stderr, allowing for case where it's very large
            tmp_stderr = open( tmp, 'rb' )
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read( buffsize )
                    if not stderr or len( stderr ) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stderr.close()
            if returncode != 0:
                raise Exception, stderr
        except Exception, e:
            sys.stderr.write( 'Error running Samtools pileup tool\n' + str( e ) )
    finally:
        #clean up temp files
        if os.path.exists( tmpDir ):
            shutil.rmtree( tmpDir ) 
    os.system('rm -rf '+options.output+'/'+inb_file) 
    print 'BAM file converted to pileup file'
    #################################################### 3 #############################################
    wig_hand = open (options.output+'/'+f_name+'.wig', 'w+')
    try:
        inh_pileup = open( pileup_file, 'rU')
        wig_hand.write('track type=wiggle_0 name=\"'+options.expName+'\"\n')
    except:
        sys.stderr.write('Cannot open pileup file, Check whether --expName option is enabled. program terminating... \n') 
        os.system('rm -rf '+pileup_file)
        sys.exit(-1)
    chr = ''    
    header_flag = 0
    for l in inh_pileup:
        l = l.strip()
        l = l.split('\t')
        if len(l) != 6:
            sys.stderr.write('Strange pileup format, instead of 6-column. Cannot continue. Program terminating ...\n')
            sys.exit(-1)
        if header_flag == 0 or chr != l[0]:
            wig_hand.write('variableStep chrom='+l[0]+'\n')
            header_flag = 1
            if int(l[3]) > 0:
                wig_hand.write( l[1]+' '+l[3] +'\n')
            chr = l[0]
            continue
        if int(l[3]) > 0:
            wig_hand.write( l[1]+' '+l[3] +'\n')
        chr = l[0]
    inh_pileup.close()
    wig_hand.close()
    os.system('rm -rf '+pileup_file)
    print 'Pileup file converted to Wiggle variableStep file.'

if __name__=="__main__": __main__()
