#############################################################################################
#                                                                                           #
#    This class is part of the MLB-Galaxy package, adding some sequence analysis            #
#    functionality to PSU's Galaxy framework.                                               #
#    Copyright (C) 2008 Cheng Soon Ong <chengsoon.ong@tuebingen.mpg.de>                     #
#    Copyright (C) 2008 Gunnar Raetsch <Gunnar.Raetsch@tuebingen.mpg.de>                    #
#    Copyright (C) 2007 Sebastian J. Schultheiss <sebi@umich.edu>                           #
#                                                                                           #
#    This program is free software; you can redistribute it and/or modify                   #
#    it under the terms of the GNU General Public License as published by                   #
#    the Free Software Foundation; either version 3 of the License, or                      #
#    (at your option) any later version.                                                    #
#                                                                                           #
#    This program is distributed in the hope that it will be useful,                        #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                         #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                           # 
#    GNU General Public License for more details.                                           #
#                                                                                           #
#    You should have received a copy of the GNU General Public License                      # 
#    along with this program; if not, see http://www.gnu.org/licenses                       #
#    or write to the Free Software Foundation, Inc., 51 Franklin Street,                    #
#    Fifth Floor, Boston, MA 02110-1301  USA                                                #
#                                                                                           #
#############################################################################################
#                                                                                           #
#  Original Author: Sebastian J. Schultheiss, version 0.77                                  #
#  Please add a notice of any modifications here:                                           #
#     Gunnar Raetsch: rewrote code for training on sequences to be read from files          #
#     Cheng Soon Ong: Added code for educational toolbox                                    #
#     Gunnar Raetsch: rewrote portions of the code                                          #
#                                                                                           #
#############################################################################################

import sys
import parse
import random
import plots
from numpy.random import randn
from numpy import ones, concatenate, array, transpose

class MotifDataDef(object):
    motif = ''
    numseq = 0
    seqlenmin = 0
    seqlenmax = 0
    posstart = 0
    posend = 0
    mutrate = 0.0

################################################################################
# data generation functions

def motifgen(motif, numseq, seqlenmin, seqlenmax, posstart, posend, mutrate):
    """Generate sequences with a particular motif at a particular location.
    Also allow a possible mutation rate of the motif.
    """

    metadata = 'motifgen(%s,%d,%d,%d,%d,%d,%1.2f)' % (motif, numseq, seqlenmin, seqlenmax, posstart, posend, mutrate)

    acgt='acgt'
    seqlist = []
    for i in xrange(0,numseq):
        str=[] ;
        seqlen=random.randint(seqlenmin,seqlenmax) ;
        for l in xrange(0,seqlen):
            str.append(acgt[random.randint(0,3)])
        pos=random.randint(posstart,posend) ;
        for l in xrange(0,len(motif)):
            if (random.random()>=mutrate) and (pos+l<seqlen) and (pos+l>=0):
                str[pos+l]=motif[l]
        seqlist.append(''.join(str))

    return metadata, seqlist


def cloudgen(numpoint, numfeat, fracpos, width):
    """Generate two Gaussian point clouds, centered around one and minus one."""

    numpos = int(round(fracpos*numpoint))
    numneg = numpoint - numpos

    metadata = 'cloudgen(%d,%d,%d,%3.2f)' % (numpos, numneg, numfeat, width)

    datapos = concatenate((ones((numpos,1)),ones((numpos, numfeat)) + width*randn(numpos, numfeat)),axis=1)
    dataneg = concatenate((-ones((numneg,1)), -ones((numneg, numfeat)) + width*randn(numneg, numfeat)), axis=1)

    pointcloud = list(datapos) + list(dataneg)
    
    return metadata, pointcloud

################################################################################
# ARFF functions

def arffwrite_real(filename, numpoint, numfeat, fracpos=0.5, width=1.0):
    """Write an ARFF file containing a vectorial dataset"""
    import arff

    (metadata, pointcloud) = cloudgen(numpoint, numfeat, fracpos, width)
    alist = [('label',1,[])]
    for ix in xrange(numfeat):
        attname = 'att%d' % ix
        alist.append((attname,1,[]))

    f = open(filename,'w')
    arff.arffwrite(f,alist,pointcloud,name='pointcloud',comment=metadata)
    f.close()
    return (pointcloud, metadata)

    
def arffwrite_sequence(filename,p, n):
    """Write an ARFF file containing a sequence dataset"""
    import arff

    (metadatapos,seqlistpos) = motifgen(p.motif, p.numseq, p.seqlenmin, p.seqlenmax, p.posstart, p.posend, p.mutrate)
    (metadataneg,seqlistneg) = motifgen(n.motif, n.numseq, n.seqlenmin, n.seqlenmax, n.posstart, n.posend, n.mutrate)

    seqlist = zip(ones(len(seqlistpos)),seqlistpos) + zip(-ones(len(seqlistneg)),seqlistneg)
    alist = [('label',1,[]),('sequence',0,[])]
    f = open(filename,'w')
    arff.arffwrite(f,alist,seqlist,name='motif',comment=metadatapos+' '+metadataneg)
    f.close()

def arffread_real(filename):
    """Read an ARFF file containing a vectorial dataset"""
    import arff

    f = open(filename,'r')
    (dataname,issparse,alist,data) = arff.arffread(f)
    f.close()
    if (alist[0][0]!='label'):
        sys.stderr.write('First column of ARFF file needs to be the label\n')
        sys.exit(-1)

    all_labels = [ex[0] for ex in data]
    all_labels = array(all_labels)
    all_examples = [ex[1:] for ex in data]
    all_examples = transpose(array(all_examples))

    print '%d features, %d examples' % all_examples.shape
    print '%d labels' % len(all_labels)

    return all_examples, all_labels


def arffread_sequence(filename):
    """Read an ARFF file containing a sequence dataset"""
    import arff

    f = open(filename,'r')
    (dataname,issparse,alist,data) = arff.arffread(f)
    f.close()
    if (alist[0][0]!='label'):
        sys.stderr.write('First column of ARFF file needs to be the label\n')
        sys.exit(-1)

    all_labels = [ex[0] for ex in data]
    all_labels = array(all_labels)
    all_examples = [ex[1].upper() for ex in data]

    print 'sequence length = %d, %d examples' % (len(all_examples[0]),len(all_examples))
    print '%d labels' % len(all_labels)

    return all_examples, all_labels


def arffread(kernelname,datafilename):
    """Decide based on kernelname whether to read a sequence or vectorial file"""

    if kernelname == 'gauss' or kernelname == 'linear' or kernelname == 'poly' or kernelname == None:
        (examples, labels) = arffread_real(datafilename)
    elif kernelname == 'wd' or kernelname == 'localalign' or kernelname == 'localimprove' or kernelname == 'spec':
        (examples, labels) = arffread_sequence(datafilename)
    else:
        print 'Unknown kernel in arffread'

    return examples, labels

################################################################################
# fasta functions

def fastawrite_sequence(filename,p):
    """Write a FASTA file containing a sequence dataset"""

    (metadata,seqlist) = motifgen(p.motif, p.numseq, p.seqlenmin, p.seqlenmax, p.posstart, p.posend, p.mutrate)

    f = open(filename,'w')
    for i in xrange(len(seqlist)):
        f.write( ">%s_sequence%i\n" % (metadata,i+1) )
        f.write( "%s\n" % seqlist[i] ) 
    f.close()

def fastaread(fnamepos,fnameneg=None):
    """Read two fasta files, the first positive, the second negative"""
    fpos = open(fnamepos,'r')
    (fa1,head1) = read_fasta_list(fpos,fasta=[],headers=[])
    fpos.close()

    if fnameneg is not None:
        fneg = open(fnameneg,'r')
        (fa2,head2) = read_fasta_list(fneg,fasta=[],headers=[])
        fneg.close()

        print 'positive: %d, negative %d' % (len(fa1),len(fa2))
        all_labels = concatenate((ones(len(fa1)),-ones(len(fa2))))
        all_examples = fa1 + fa2
    else:
        all_examples = fa1
        all_labels = ones(len(fa1))

    return all_examples, all_labels

def read_fasta(f, fasta=dict(), id_only=True):
    """ read fasta as dictionary """
    import numpy
    str=f.read()
    idx = str.find('>')
    if idx==-1:
        # file name list?
        files = str.split('\n') ;
        for s in files:
            if len(s)==0: continue
            print s.strip() + '\n'
            fasta = read_fasta(file(s.strip()), fasta) 
    else:
        # real fasta file?
        sequences = str.split('>') 
        for s in sequences:
            if len(s)==0: continue
            header = s[0:s.index('\n')]
            if id_only:
                header = header.split(' ')[0] ;
                header = header.split('\t')[0] ;
            seq = s[s.index('\n')+1:].replace('\n','').upper()
            #print 'has_key', fasta.has_key(header),header
            if fasta.has_key(header):
                print "duplicate key " + header
            assert(not fasta.has_key(header))
            fasta[header]=seq ;

    return fasta

def read_fasta_list(f, fasta=[], headers=[], id_only=True):
    """ read fasta as dictionary """
    import numpy
    str=f.read()
    idx = str.find('>')
    if idx==-1:
        # file name list?
        files = str.split('\n') ;
        for s in files:
            if len(s)==0: continue
            print s.strip() + '\n'
            (fasta,header) = read_fasta_list(file(s.strip()), fasta, headers) 
    else:
        # real fasta file?
        sequences = str.split('>') 
        for s in sequences:
            if len(s)==0: continue
            header = s[0:s.index('\n')]
            if id_only:
                header = header.split(' ')[0] ;
                header = header.split('\t')[0] ;
            seq = s[s.index('\n')+1:].replace('\n','').upper()
            fasta.append(seq) ;
            headers.append(header) ;

    return (fasta,headers)

################################################################################
# main

if __name__ == '__main__':

    #if 1:
    #    sys.stderr.write( "xxx %i %s %s\n" % (len(sys.argv),sys.argv[7], sys.argv[8]) )
    #    sys.exit(-1)

    if len(sys.argv)<3 or (sys.argv[1]=='motif' and sys.argv[2]!='arff' and sys.argv[2]!='fasta') or (sys.argv[1]=='motif' and sys.argv[2]=='fasta' and len(sys.argv)<9) or (sys.argv[1]=='motif' and sys.argv[2]=='arff' and len(sys.argv)<14) or (sys.argv[1]=='cloud' and len(sys.argv)<7) or (sys.argv[1]!='motif') and (sys.argv[1]!='cloud'):
        sys.stderr.write( "usage: %s motif fasta MOTIF numSeq seqLenRange positionRange mutationRate output.fa\n   or: %s motif arff MOTIFPOS numSeq-pos seqLenRange-pos positionRange-pos mutationRate-pos \\\n                             motif-neg numSeq-neg seqLenRange-neg positionRange-neg mutationRange-neg output.arff\n   or: %s cloud numpoints dimensions fractionOfPositives cloudWidth output.arff\n" % (sys.argv[0],sys.argv[0],sys.argv[0]) )
        sys.exit(-1)

    random.seed()

    if sys.argv[1] == 'motif':
        if sys.argv[2]=='fasta':
            # generate sequences in FASTA format
            p = MotifDataDef()
            p.motif = sys.argv[3]
            p.numseq = int(sys.argv[4])
            (p.seqlenmin,p.seqlenmax) = parse.parse_range(sys.argv[5])
            (p.posstart,p.posend) = parse.parse_range(sys.argv[6])
            p.mutrate = float(sys.argv[7])

            filename = sys.argv[8]
            fastawrite_sequence(filename, p)

        else:
            # generate sequences in ARFF format
            assert(sys.argv[2]=='arff')
            p = MotifDataDef()
            p.motif = sys.argv[3]
            p.numseq = int(sys.argv[4])
            (p.seqlenmin,p.seqlenmax) = parse.parse_range(sys.argv[5])
            (p.posstart,p.posend) = parse.parse_range(sys.argv[6])
            p.mutrate = float(sys.argv[7])
            
            n = MotifDataDef()
            n.motif = sys.argv[8]
            n.numseq = int(sys.argv[9])
            (n.seqlenmin,n.seqlenmax) = parse.parse_range(sys.argv[10])
            (n.posstart,n.posend) = parse.parse_range(sys.argv[11])
            n.mutrate = float(sys.argv[12])
            
            filename = sys.argv[13]
            arffwrite_sequence(filename, p, n)
        
    elif sys.argv[1] == 'cloud':
        # generate a data cloud in ARFF format
        numpoint = int(sys.argv[2])
        numfeat = int(sys.argv[3])
        fracpos = float(sys.argv[4])
        width = float(sys.argv[5])

        filename = sys.argv[6]
	(pointcloud,metadata)=arffwrite_real(filename, numpoint, numfeat, fracpos, width)
	
	if len(sys.argv)>=8:
		plots.plotcloud(pointcloud,sys.argv[7],metadata)
    else:
        print 'Unknown option %s\n' % sys.argv[1]
