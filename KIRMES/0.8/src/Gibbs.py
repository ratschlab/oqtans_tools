#############################################################################################
#                                                                                           #
#    This class is part of the metaTFBS ExpressionAnalysis package.                         #
#    Copyright (C) 2006 - 2007 Sebastian J. Schultheiss <sebi@umich.edu>                    #
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

class Priority(MotifFinder.MotifFinder):
    """Interface to INCLUSive Motif suite of tools"""
    def __init__(self, binpath = "./"):
        MotifFinder.MotifFinder.__init__(self, binpath)

    def __del__(self):
        MotifFinder.MotifFinder.__del__(self)
        

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
