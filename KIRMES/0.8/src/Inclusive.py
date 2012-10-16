"""#############################################################################################
#                                                                                           #
#    This class is part of the metaTFBS ExpressionAnalysis package.                         #
#    Copyright (C) 2006 - 2010 Sebastian J. Schultheiss <sebi@umich.edu>                    #
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
#  Original Author: Sebastian J. Schultheiss, version 0.7.7                                 #
#  Please add a notice of any modifications here:                                           #
#                                                                                           #
#                                                                                           #
#############################################################################################
"""

from Inputs import GffFile, FastaFile
import MotifFinder
import os

class Inclusive(MotifFinder.MotifFinder):
    """Interface to INCLUSive suite of motif tools"""
    def __init__(self, binpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../bin/")):
        """Default constructor, using this module's path as base for the bin directory"""
        MotifFinder.MotifFinder.__init__(self, binpath = binpath)
        self.immfilename = None
        self.fastafilename = None
        self.bgfilename = None
        self.unique = False

    def __del__(self):
        """Destructor, deleting temporary files"""
        self.rmTempFile([self.bgfilename])
        MotifFinder.MotifFinder.__del__(self)
        
    def makeUniqueFasta(self):
        """Opens the current self.fastafilename, removes non-unique headers and rewrites it"""
        fastafilename = self.tempFile("fasta")
        tfasta = FastaFile()
        tfasta.read(self.fastafilename, "", True)
        tfasta.getFormattedFasta(fastafilename)
        self.fastafilename = fastafilename
        self.unique = True

    def createBackground(self, order = 1, organismname = "organism"):
        """Interface to the CreateBackgroundModel program, you can specify the model order."""
        if self.bgfilename is None:
            self.bgfilename = self.tempFile("bg", "")
        if not self.unique:
            self.makeUniqueFasta()
        os.system(os.path.join(self.binpath, "CreateBackgroundModel") \
                  + " -f " + self.fastafilename + " -b " + self.bgfilename \
                  + " -o " + str(order) + " -n " + organismname + " 2> /dev/null")
        
    def sampleFiles(self, length = None, nof_motifs = None):
        """Interface to the MotifSampler program, you can specify length and number of motifs"""
        self.immfilename = self.tempFile("imm", "")
        gffoutfilename = self.tempFile("gff", "")
        if length is None:
            length = self.settings.getMotifLength()
        if nof_motifs is None:
            nof_motifs = self.settings.getNofMotifs()
        os.system(os.path.join(self.binpath, "MotifSampler") + \
                  " -f " + self.fastafilename + \
                  " -b " + self.bgfilename + " -o " + gffoutfilename + \
                  " -m " + self.immfilename + " -w " + str(length) + \
                  " -n " + str(nof_motifs) + " 2> /dev/null")
        return gffoutfilename
    
    def scanFiles(self, prior = 0.9):
        """Interface to the MotifScanner program, you can specify the prior of finding 1 motif copy"""
        gffoutfilename = self.tempFile("gff")
        os.system(os.path.join(self.binpath, "MotifScanner") +
                  " -f " + self.fastafilename + " -b " + self.bgfilename + 
                  " -m " + self.immfilename + " -o " + gffoutfilename +
                  " -p " + str(prior) + " 2> /dev/null")
        return gffoutfilename
    
    def setFastaFile(self, fastafile):
        """Set the FASTA file in which to search"""
        self.fastafilename = fastafile
    
    def findMotifs(self):
        """Fill a best-results dictionary and return it"""
        if len(self.seqs) == 0 and self.fastafilename is None:
            raise ValueError("Please add some sequences to find motifs in")
        if self.bgfilename is None:
            self.bgfilename = self.tempFile("bg", "")
            self.createBackground(self.fastafilename, self.bgfilename)
        if self.immfilename is None:
            self.sampleFiles()
        gffoutfilename = self.scanFiles(0.9999)
        gff = GffFile()
        gff.read(gffoutfilename, "")
        best_dict = gff.getAllMotifsAllGenes()
        return best_dict

    def getMotifGff(self):
        """Get a GFF with the best motifs for all sequences"""
        gff_dict = self.findMotifs()
        gff = []
        for motif in gff_dict.values():
            gff.extend(motif)
        bestfile = GffFile()
        bestfile.gff = gff
        gfffilename = self.tempFile("gff", bestfile.getFormattedGff())
        return gfffilename

    def getMatrixFile(self):
        """Returns the path name for the INCLUSive Motif Model matrix file 
        with all motifs as position weight matrices"""
        return self.immfilename

    
    
    
    
    
    
    
