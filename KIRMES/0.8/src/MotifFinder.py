"""
#######################################################################################
#                                                                                     #
#    This class is part of the KIRMES package.                                        #
#    Copyright (C) 2008 - 2009 Sebastian J. Schultheiss <sebi@umich.edu>              #
#    Please cite:                                                                     #
#    Sebastian J. Schulheiss, Wolfgang Busch, Jan U. Lohmann, Oliver Kohlbacher,      #
#    and Gunnar Raetsch (2008) KIRMES: Kernel-based identification of regulatory      #
#    modules in euchromatic sequences. In Andreas Beyer and  Michael Schroeder (Eds.) #
#    German Conference on Bioinformatics, 158-167, GI, Springer Verlag Heidelberg.    # 
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

import Inputs
import os
import random

class MotifFinderSettings(object):
    """Holds various settings for the MotifFinder base class, e.g. length of motifs"""
    
    def __init__(self, motif_length = 6, window_size = 20, replacement_character = "A", nof_motifs = 100):
        """Initialize with default values, or set them all"""
        self.setMotifLength(motif_length)
        self.setWindowSize(window_size)
        self.setNofMotifs(nof_motifs)
        self.setReplacementCharacter(replacement_character)
        
    def getMotifLength(self):
        """Returns the motif length member variable"""
        return self.motif_length
    
    def setMotifLength(self, motif_length):
        """Sets the motif length member variable"""
        assert isinstance(motif_length, int)
        self.motif_length = motif_length
        
    def setWindowSize(self, window_size):
        """Sets the window size to cut out around the motif"""
        assert isinstance(window_size, int)
        self.window_size = window_size
        
    def getWindowSize(self):
        """Returns the current window size member variable"""
        return self.window_size

    def getNofMotifs(self):
        """Returns the current number of motifs member variable"""
        return self.nof_motifs
    
    def setNofMotifs(self, nof_motifs):
        """Sets the number of motifs to search for"""
        assert isinstance(nof_motifs, int)
        self.nof_motifs = nof_motifs
        
    def setReplacementCharacter(self, replacement_character):
        """Replace any non-nucleotide character in the input sequence
        with replacement_character, one of A, C, G, T, or R (randomly)"""
        assert len(replacement_character) == 1
        assert isinstance(replacement_character, str)
        if replacement_character not in "ACGTR":
            raise ValueError("replacement_character has to be one of A, C, G, T, or R (random)")
        self.replace = replacement_character
        
    def getReplacementCharacter(self):
        """Returns the current replacement character"""
        return self.replace
    
class MotifFinder(object):
    """Base class for motif finder programs and the Kmer class"""

    def __init__(self, fasta = Inputs.FastaFile(), finder_settings = MotifFinderSettings(), binpath = ""):
        """Initialize with empty objects or pass them upon class instantiation"""
        self.settings = finder_settings
        self.seqs = fasta
        self.results = None
        self.motifs = None 
        self.binpath = binpath #os.path.split(binpath)
        self.filelist = []
        random.seed()
        self.nof_motifs = -1

    def __del__(self):
        self.rmTempFile(self.filelist)
    
    def tempFile(self, extension, contents = ""):
        """Create a temp file name in the style of
        /tmp/motif_tmp_198123498.txt
        on a POSIX system or on NT:
        \tmp\motif_tmp_198123498.txt"""
        tmpint = str(random.randint(100000000, 999999999))
        fname = os.path.join(os.sep, "tmp", "motif_tmp_" + tmpint + os.extsep + extension)
        f = open(fname, 'w')
        f.writelines(contents)
        f. close()
        self.filelist.append(fname)
        return fname   
    
    def rmTempFile(self, filenames = None):
        """Will delete all tmp files created by this class instance.
        Can also delete individual files or lists of files.
        If they have been created by this class instance, they are
        also removed from the file list of this class."""
        if (type(filenames) == type("str")):
            os.remove(filenames)
            return
        if filenames is None:
            filenames = self.filelist
        for f2del in filenames:
            try:
                os.remove(f2del)
            except OSError: # Can't delete this file ...
                pass # ... let's try to clean up the rest anyway       
            try:
                self.filelist.remove(f2del)
            except ValueError: # It's not in the list ...
                pass # ... let's continue anyway
    
    def getMotifLength(self):
        """Returns the motif length of the motif finder settings member"""
        return self.settings.getMotifLength()    
    
    def setFastaFile(self, fastafile):
        """Set the FASTA file in which to search"""
        self.seqs = Inputs.FastaFile(fastafile)
        self.seqs.nucleotideReplace(self.settings.getReplacementCharacter())
        
    def setMotifFinderSettings(self, finder_settings):
        """Pass a new MotifFinderSettings object"""
        self.settings = finder_settings
        
    def getResults(self): #, top_n = -1):
        """Returns the motif windows in a dict and the position vector
        as a list"""
        return self.getMotifWindow()
        
    def findMotifs(self):
        """Should be overwritten in the descending classes"""
        pass
    
    def setMotifs(self, gff):
        """Read a pre-computed GFF file with motif positions"""
        gfffile = Inputs.GffFile()
        gfffile.read(gff)
        self.motifs = gfffile.getAllMotifsAllGenes()
            
    def getMotifWindow(self): #, top_n): # = self.nof_motifs):
        """Returns the cut-out windows for all sequences 
        for the top n motifs"""
        genes = self.getGenes()
        winsize = self.settings.getWindowSize()
        windows = {}
        for motif in self.motifs.keys():
            windows[motif] = [0] * len(genes)
        positions = {}
        for gene in genes:
            positions[gene] = {}
        warn = 0
        for motif, occurrences in self.motifs.iteritems():
            for ge in occurrences:
                center = ge.getStart() + (ge.getStop() - ge.getStart()) / 2
                window = self.seqs.getSequenceWindow(ge.getGene(), center, winsize)
                assert len(window) == winsize 
                positions[ge.getGene()][motif] = center
                #if this fails, it's on the - strand and we need to reverse-complement the window
                if not ge.getStrand():
                    window = Inputs.quickRevComp(window)
                windows[motif][genes.index(ge.getGene())] = window
                assert len(window) == winsize
            # this should not be necessary
            while 0 in windows[motif]:
                k = windows[motif].index(0)
                windows[motif][k] = "A" * winsize
                warn += 1
        if warn > 0:
            print "#Note: Not all motifs found in all genes."
            print """#      Check your GFF: Each motif has to be listed with at least 
#      one location in every sequence of the FASTA file. Occurrences:"""
            print "#     ", str(warn)
        position_vec = self.getPositionVector(positions)
        return windows, position_vec

    def getGenes(self):
        """Return the FASTA identifiers for all involved genes in a list"""
        return self.seqs.getGenes()
    
    def formatVector(self, positions):
        """Formats the vector of positions so the position for every
        kmer and pairwise distance is at the same entry in 
        the positions vector."""
        motifs = positions.values()[0].keys()
        motifs.sort()
        alphabetical = []
        for gene in positions.keys():        
            temp_list = []
            for motif in motifs:
                if motif not in positions[gene]:
                    maxlen = self.seqs.getFastaEntry(gene).getSequenceLength()
                    positions[gene][motif] = random.randint(0, maxlen)
                temp_list.append(positions[gene][motif])
            alphabetical.append(temp_list)
        return alphabetical 

    
    def getPositionVector(self, positions):
        """Return a vector with all pairwise distances and positions of
        kmers within all sequences"""
        for gene, motif_dict in positions.iteritems():
            motifs = motif_dict.keys()
            while len(motifs) > 1:
                for i in xrange(len(motifs) - 1):
                    motif_dict[motifs[0] + "-" + motifs[i + 1]] = abs(motif_dict[motifs[0]] - motif_dict[motifs[i + 1]])
                motifs.pop(0)
            positions[gene] = motif_dict
        alphabetical = self.formatVector(positions)
        return alphabetical
    
