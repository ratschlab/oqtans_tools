"""
#######################################################################################
#                                                                                     #
#    kmotif.py is a command-line front-end to the KIRMES pipeline                     #
#    BibTeX entries below. Please cite:                                               #
#    Sebastian J. Schulheiss, Wolfgang Busch, Jan U. Lohmann, Oliver Kohlbacher,      #
#    and Gunnar Raetsch (2009) KIRMES: Kernel-based identification of regulatory      #
#    modules in euchromatic sequences. Bioinformatics 16(25):2126-33.                 # 
#                                                                                     #
#    Copyright (C) 2007-20010 Sebastian J. Schultheiss <sebi@umich.edu>               #
#                                                                                     #
#    This program is free software; you can redistribute it and/or modify             #
#    it under the terms of the GNU General Public License as published by             #
#    the Free Software Foundation; either version 3 of the License, or                #
#    (at your option) any later version.                                              #
#                                                                                     #
#    This program is distributed in the hope that it will be useful,                  #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                   #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                     # 
#    GNU General Public License for more details.                                     #
#                                                                                     #
#    You should have received a copy of the GNU General Public License                # 
#    along with this program; if not, see http://www.gnu.org/licenses                 #
#    or write to the Free Software Foundation, Inc., 51 Franklin Street,              #
#    Fifth Floor, Boston, MA 02110-1301  USA                                          #
#                                                                                     #
#######################################################################################
#                                                                                     #
#  Original Author: Sebastian J. Schultheiss, version 0.8.0                           #
#  Please add a notice of any modifications here:                                     #
#                                                                                     #
#                                                                                     #
#######################################################################################
"""

__version__ = "0.8.0"
__license__  = "GNU General Public License"

# bibtex entry
__author__ = """
@article{Schultheiss2009KIRMES,
    author = {Sebastian J. Schultheiss and Wolfgang Busch and Jan U. Lohmann and Oliver Kohlbacher and Gunnar Raetsch},
    title = {{KIRMES}: Kernel-based identification of regulatory modules in euchromatic sequences},
    year = {2009},
    journal = {Bioinformatics},
    publisher = {Oxford Journals},
    volume = {25},
    issue = {16},
    pages = {2126--33},
    month = {April},
    doi = {10.1093/bioinformatics/btp278},
    abstract = {Motivation: Understanding transcriptional regulation is one of the 
                main challenges in computational biology. An important problem is 
                the identification of transcription factor binding sites in promoter 
                regions of potential transcription factor target genes. It is 
                typically approached by position weight matrix-based motif 
                identification algorithms using Gibbs sampling, or heuristics to 
                extend seed oligos. Such algorithms succeed in identifying single, 
                relatively well-conserved binding sites, but tend to fail when it 
                comes to the identification of combinations of several degenerate 
                binding sites, as those often found in cis-regulatory modules.  
                Results: We propose a new algorithm that combines the benefits of 
                existing motif finding with the ones of Support Vector Machines (SVMs) 
                to find degenerate motifs in order to improve the modeling of 
                regulatory modules. In experiments on microarray data from Arabidopsis 
                thaliana, we were able to show that the newly developed strategy 
                significantly improves the recognition of transcription factor targets.  
                Availability: The PYTHON source code (open source-licensed under GPL), 
                the data for the experiments and a Galaxy-based web service are 
                available at http://www.fml.mpg.de/raetsch/projects/kirmes.  
                Contact: sebastian.schultheiss@tuebingen.mpg.de},
    URL = {http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btp278v1}
}
"""

__usage__ = """Usage: 
  To find motifs in 2 FASTA files, supply positives and negatives: 
        %prog -t -p positives.fasta -n negatives.fasta [options]
"""
# system imports
import os
from optparse import OptionParser, OptionValueError
import sys
from shutil import copy

try:
    # own imports
    import kirmes_ini
    from Inclusive import Inclusive
    from Kmers import Kmers, createIMMFile
except ImportError:
    print "ImportError:"
    print "One of the required imports failed, please make sure the files "
    print "kirmes_ini.py or KIRMES_INI.pyc, "
    print "DBTools.py or DBTools.pyc and FileTools.py or FileTools.pyc"
    print "are present in the current directory, or download this program again."
    raise

def check_file(option, opt_str, value, _parser):
    """See if a file exists on the file system, raises an OptionValueError"""
    if value == "None":
        value = None
    elif not os.path.isfile(value):
        raise OptionValueError("Cannot open %s as a file. Please check if it exists." % value)
    setattr(_parser.values, option.dest, value)

def check_pickle(option, opt_str, value, _parser):
    """Check for the kernel file in testing mode"""
    if not _parser.values.train:
        check_file(option, opt_str, value, _parser)    
    else:
        setattr(_parser.values, option.dest, value)

def optionparse(parser):
    """Completes the option parser object, adding defaults and options"""
    parser.add_option("-t", "--type", type = "string",  dest = "type",
                          help = "motif finding strategy to use: MotifSampler (kims), PRIORITY (krgp), k-mer (kkmc), or just the locator, must supply a valid imm file (kiml) [default %default]")
    parser.add_option("-p", "--positives", dest = "positives",
                          action = "callback", callback = check_file, type = "string", 
                          help="path to the fasta file with a positive set of regulatory regions [default %default]")
    parser.add_option("-n", "--negatives", dest = "negatives", type = "string", 
                          action = "callback", callback = check_file, 
                          help="path to the fasta file with a negative set of regulatory regions [default %default]")
    parser.add_option("-i", "--pgff", dest = "pgff", type = "string", 
                          help="path to the output gff file of motif positions from the positive regulatory regions [default %default]")
    parser.add_option("-j", "--ngff", dest = "ngff", type = "string", 
                          help="path to the output gff file of motif positions from the negative regulatory regions [default %default]")
    parser.add_option("-x", "--matrix", dest = "imm", type = "string", 
                          help="path to the input or output imm file of motif motif models as position weight matrices [default %default]")
    parser.add_option("-m", "--motifs", type = "int", dest = "nof_motifs", 
                          help = "number of motifs to consider [default %default]")
    parser.add_option("-l", "--length", type = "int", dest = "motif_length", 
                          help = "length of the motifs to search for [default %default]")
    parser.set_defaults(positives = kirmes_ini.POSITIVES_FILENAME, 
                        negatives = kirmes_ini.NEGATIVES_FILENAME,
                        nof_motifs = kirmes_ini.NOF_MOTIFS,
                        motif_length = kirmes_ini.MOTIF_LENGTH,
                        type = kirmes_ini.SAMPLING_STRATEGY,
                        ngff = kirmes_ini.NGFF_FILENAME,
                        pgff = kirmes_ini.PGFF_FILENAME,
                        imm = kirmes_ini.IMM_FILENAME)

def motifScan(fastafile, matrixfile, gfffile):
    """Search for motifs with existing matrix defintion"""
    ive = Inclusive()
    ive.fastafilename = fastafile
    ive.immfilename = matrixfile
    gff = ive.getMotifGff()
    copy(gff, gfffile)    

def kims(options):
    """Run the MotifSampler Program"""
    ive = Inclusive()
    ive.fastafilename = options.positives
    ive.settings.setMotifLength(options.motif_length)
    ive.settings.setNofMotifs(options.nof_motifs)
    pgff = ive.getMotifGff()
    copy(pgff, options.pgff)
    #second round, find motifs in negative sequences
    ive.fastafilename = options.negatives
    ngff = ive.getMotifGff()
    copy(ngff, options.ngff)
    imm = ive.getMatrixFile()
    copy(imm, options.imm)
        
def krgp(options):
    """Run the Priority Program (falls back to Inclusive)"""
    kims(options)
    
def kkmc(options):
    """Run the kmer counting strategy, search for motifs with KIML"""
    mer = Kmers()
    mer.setFastaFile(options.positives)
    mer.settings.setMotifLength(options.motif_length)
    pkmerdict = mer.countKmers(options.nof_motifs)
    createIMMFile(options.imm, pkmerdict.keys())
    motifScan(options.positives, options.imm, options.pgff)
    motifScan(options.negatives, options.imm, options.ngff)

def kiml(options):
    """Search for motifs with existing matrix defintion"""
    motifScan(options.positives, options.imm, options.pgff)
    if options.negatives:
        motifScan(options.negatives, options.imm, options.ngff)
        
#######################
# main                #
#######################
def main(argv = None):
    """main() block"""     
    if argv is None:         
        argv = sys.argv    
    parser = OptionParser(version = "%prog " + __version__, usage = __usage__)
    optionparse(parser)
    (options, args) = parser.parse_args()
    if options.type == "krgp":
        krgp(options)
    elif options.type == "kkmc":
        kkmc(options)
    elif options.type == "kiml":
        kiml(options)
    else:
        kims(options)

if __name__ == "__main__":
    main()





   
    
    

