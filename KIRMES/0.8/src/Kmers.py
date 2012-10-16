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

import kirmes_ini
import MotifFinder
import re
from heapq import nlargest, itemgetter

def createIMM(kmer, pseudocount = 0.003, id = ""):
    """Creates an INCLUSive Motif Model string from a kmer sequence, with pseudocounts"""
    kmer = kmer.upper()
    nonnuc = re.findall("[^ACGT]", kmer)
    try:
        assert len(nonnuc) == 0
    except AssertionError:
        print "Non-nucleotide characters found in your kmer sequence:"
        print nonnuc
        raise
    imm = "#INCLUSive Motif Model created from a kmer with KIRMES\n#\n#ID = "
    if id == "":
        imm += "kmer_" + kmer
    else:
        imm += str(id) + "_" + kmer
    imm += "\n#Score = 1.0\n#W = " + str(len(kmer))
    imm += "\n#Consensus = " + kmer + "\n"
    count = 1.0 - 3 * pseudocount
    dic = {"A": 0, "C": 1, "G": 2, "T": 3}
    for character in kmer:
        quad = 4 * [pseudocount]
        quad[dic[character]] = count
        imm += str(quad[0]) + "\t" + str(quad[1]) + "\t" + str(quad[2]) + "\t" + str(quad[3]) + "\n" 
    imm += "\n"
    return imm

def createIMMFile(filename, kmers, pseudocount = 0.003):
    """Write out an IMM file using a list of kmers"""
    imm = ""
    for i, kmer in enumerate(kmers):
        imm += createIMM(kmer, pseudocount, "kmer_" + str(i + 1))
    immfile = open(filename, "w")
    immfile.write(imm)
    immfile.close()
        

class Kmers(MotifFinder.MotifFinder):
    """Implements motif finding through counting k-mers"""
    def __init__(self):
        MotifFinder.MotifFinder.__init__(self)
        self.blur_degree = 1
        
    def getSequences(self):
        """Return the sequence operating on"""
        return self.seqs
        
    def setKmerLength(self, length):
        """Manually set the length of k-mers to count"""
        self.settings.setMotifLength(length)
    
    def getKmerLength(self):
        """Report the currently set length of k-mers to count"""
        return self.settings.getMotifLength()
    
    def findMotifs(self, kmer_dict = None):
        """Fill the self.results object and return"""
        if len(self.seqs) == 0:
            raise ValueError("Please add some sequences to find motifs in")
        self.results = []
        if kmer_dict is None:
            kmer_dict = self.countKmers()
        kmers = kmer_dict.keys()
        occurrences = kmer_dict.values()
        for seqno in xrange(len(self.seqs)):
            indices = self.findKmers(self.seqs.getiSequence(seqno).upper(), kmers)
            for i in xrange(len(kmers)):
                self.results.append(MotifFinder.MotifFinderResult(seqno, kmers[i], occurrences[i], indices[i]))
        return True
    
    def countKmers(self, top_n = -1):
        """Puts all k-mers occurring in a string into a dictionary 
        and returns it, with values containing the number of occurrences."""
        kmers = {}
        for seq in self.seqs.getAllSequences():
            seq = seq.upper()
            for i in xrange(len(seq)):
                kmer = seq[i:i + self.getKmerLength()]
                if len(kmer) == self.getKmerLength():
                    if kmer in kmers:
                        kmers[kmer] += 1
                    else:
                        kmers[kmer] = 1
        if top_n == -1:
            return kmers
        else:
            top_list = nlargest(top_n, kmers.iteritems(), itemgetter(1))
            top_dict = {}
            for key, value in top_list:
                top_dict[key] = value
            return top_dict
    
    def getResults(self, top_n = -1):
        """Returns the n kmers which occur most frequently,
        input: a kmer dictionary with occurrence as values,
        returns: a set of pairs with the n most frequent kmers"""
        if self.results is None:
            self.findMotifs()
        if top_n == -1:
            top_n = len(self.results)
        kmers = {}
        for i in xrange(len(self.results)):
            kmers[self.results[i].getMotif()] = self.results[i].getScore()
        return nlargest(top_n, kmers.iteritems(), itemgetter(1))
    
    def blurKmer(self, kmer, blur_degree = None):
        """Returns a regex pattern string that contains a match-all (.) 
           at each position in the kmer, concatenated with or (|)""" 
        if blur_degree is None:
            blur_degree = self.blur_degree
        if blur_degree > len(kmer):
            blur_degree = len(kmer)
        blurred_kmers = [kmer]
        for deg in xrange(blur_degree):
            temp_kmers = []
            for blur_mer in blurred_kmers:
                for i in range(0, len(blur_mer)):
                    temp_kmers.append(blur_mer[:i] + "." + blur_mer[i + 1:])
            blurred_kmers = temp_kmers
        regex = ""
        for blur_mer in blurred_kmers:
            if regex != "":
                regex += "|"
            regex += blur_mer
        return regex
    
    def misMatch(self, seq, kmer):
        """Returns the index of the closest matching kmer in a string, preferably the core of the kmer"""
        idx = seq.rfind(kmer[1:])
        if idx > -1:
            return idx
        idx = seq.rfind(kmer[:-1])
        if idx > -1:
            return idx
        blur_degree = 1
        while idx < 0:
            regex = self.blurKmer(kmer, blur_degree)
            rmatch = None
            for match in re.finditer(regex, seq):
                rmatch = match
            if rmatch is None:
                blur_degree += 1
            else:
                #print "found with blur_degree " + str(blur_degree) + " at position " + rmatch.group() + " for kmer " + kmer
                idx = rmatch.start()
        return idx
                
    def findKmers(self, seq, kmers):
        """Takes a list of kmers and returns their best matches 
           in the input sequence."""
        indices = []
        edge = (kirmes_ini.MOTIF_WINDOW_WIDTH - kirmes_ini.KMER_LENGTH + 1) / 2 + 1
        for kmer in kmers:
            idx = seq[edge:-edge].rfind(kmer)
            if idx == -1:
                idx = self.misMatch(seq[edge:-edge], kmer)
            indices.append(idx + edge)
        return indices


