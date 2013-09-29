"""
#######################################################################################
#                                                                                     #
#    This class is part of the KIRMES package.                                        #
#    Copyright (C) 2008 - 2010 Sebastian J. Schultheiss <sebi@umich.edu>              #
#    Please cite:                                                                     #
#    Sebastian J. Schulheiss, Wolfgang Busch, Jan U. Lohmann, Oliver Kohlbacher,      #
#    and Gunnar Raetsch (2009) KIRMES: Kernel-based identification of regulatory      #
#    modules in euchromatic sequences. Bioinformatics 16(25):2126-33.                 # 
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

__usage__ = """Usage: 
  To return the example files ending in "gen" for generated, type:
    %prog -e gen -p positives.fasta -g positives.gff
          
  You may leave out any of the files."""

__version__ = "0.8.0"

__desc__ = {"gen": """The generated dataset contains a planted motif, GATTACA in the positive dataset. You
can generate more datasets like these with our Toy Data > Motif Gen (FASTA) tool.""",
"ath" : """The example data from Arabidopsis thaliana contains +1000 base pairs upstream from the 
gene start for those genes that reacted to a heat shock treatment with sudden over- 
or underexpression [1]. The negative dataset contains genes that were relatively highly 
expressed both in the treated and untreated plants. They seem to be unaffected by
the heat shock and are thus probably not regulated by the same mechanism as the 
genes that do react. One can find certain motifs that are common among most of the 
reacting genes, including the published heat shock element [1].

[1] Busch W, Wunderlich M, Schoeffl F. (2004) Identification of novel heat 
shock factor-dependent genes and biochemical pathways in Arabidopsis 
thaliana. Plant J 1(41):1-14."""}


# system imports
import os
from optparse import OptionParser, OptionGroup, OptionValueError
import sys
import shutil

def check_file(option, opt_str, value, _parser):
    """See if a file exists on the file system, raises an OptionValueError"""
    #if not os.path.isfile(value):
    #    raise OptionValueError("Cannot open %s as a file. Please check if it exists." % value)
    setattr(_parser.values, option.dest, value)

def optionparse(parser):
    gffgroup = OptionGroup(parser, "Options for GFF files",
                        "These options contain the paths to the GFF files to be saved.")
    gffgroup.add_option("-i", "--posgff", type = "string", action = "callback", 
                        callback = check_file, dest = "positives.gff", 
                        help="path to the positive GFF example file")
    gffgroup.add_option("-j", "--neggff", type = "string", action = "callback", 
                        callback = check_file, dest = "negatives.gff", 
                        help="path to the negative GFF example file")
    gffgroup.add_option("-u", "--predgff", type = "string", action = "callback", 
                        callback = check_file, dest = "for-prediction.gff", 
                        help="path to the unlabeled GFF example file")    
    parser.add_option_group(gffgroup)
    
    fastagroup = OptionGroup(parser, "Options for FASTA files",
                        "These options contain the paths to the FASTA files to be saved.")
    fastagroup.add_option("-p", "--positives", dest = "positives.fasta",
                          action = "callback", callback = check_file, type = "string", 
                          help="path to the FASTA file with a positive set of regulatory regions")
    fastagroup.add_option("-n", "--negatives", dest = "negatives.fasta", type = "string", 
                          action = "callback", callback = check_file, 
                          help="path to the FASTA file with a negative set of regulatory regions")
    fastagroup.add_option("-q", "--prediction", dest = "for-prediction.fasta", type = "string", 
                          action = "callback", callback = check_file, 
                          help="path to the FASTA file with an unlabeled set of regulatory regions")
    parser.add_option_group(fastagroup)

    immgroup = OptionGroup(parser, "Option for IMM matrix file",
                           "This option contains the paths to the IMM file to be saved.")
    fastagroup.add_option("-x", "--matrix", dest = "matrix",
                          action = "callback", callback = check_file, type = "string", 
                          help="path to the IMM file with a set of motifs found in the FASTAs")
    parser.add_option_group(immgroup)
        
    parser.add_option("-e", "--example", dest = "eset", type = "string", 
                      help = "determines which set of example files to return (ath, gen, ...)")
    parser.add_option("-d", "--exampledir", dest = "edir", type = "string", 
                      help = "directory for the example files [auto-default: %default]")
    parser.set_defaults(edir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../examples/"), 
                        eset = "gen")

#######################
# main                #
#######################
def main(argv = None):     
    if argv is None:         
        argv = sys.argv
    parser = OptionParser(version = "%prog " + __version__, usage = __usage__)
    optionparse(parser)
    (options, args) = parser.parse_args()
    print options.edir
    print __desc__[options.eset]
    for ext in (".fasta", ".gff"):
        for kind in ("positives", "negatives", "for-prediction"):
            if options.__dict__[kind + ext] is not None:
                shutil.copy(os.path.join(options.edir, str(kind + "-" + options.eset + ext)), 
                            str(options.__dict__[kind + ext]))
    shutil.copy(os.path.join(options.edir, "matrix-" + str(options.eset) + ".imm"), str(options.matrix))

if __name__ == '__main__':
    main()
