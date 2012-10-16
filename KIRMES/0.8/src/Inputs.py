"""
#############################################################################################
#                                                                                           #
#    This class is part of the metaTFBS ExpressionAnalysis package, it manages input files. #
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
#  Original Author: Sebastian J. Schultheiss, version 0.77                                  #
#  Please add a notice of any modifications here:                                           #
#                                                                                           #
#                                                                                           #
#############################################################################################
"""

import os
from re import sub
from random import randint

def quickRevComp(dna):
    """A simple reverse complement function for DNA sequences; no error checking"""
    table = {"A":"T", "T":"A", "G":"C", "C":"G"}
    ret = ""
    for c in dna.upper():
        if c in table:
            ret += table[c]
        else:
            ret += c
    return ret[::-1]

class Genelist:
    """Defines handling of gene lists and tools"""
    
    def __init__(self):
        """Initialize the only member, a list of genes"""
        self.list = []
        
    def addList(self, listfile):
        """Extend the current list with another, if it's a filename, open and read"""
        if (isinstance(listfile, file)):
            addlist = listfile.readlines()
        else:
            lfile = open(listfile)
            addlist = lfile.readlines()
        for entry in addlist:
            if (len(entry)> 2):
                self.list.append(entry.strip())
    
    def getList(self):
        """Return the list member"""
        return self.list
    
    def add(self, entry):
        """Adds an entry to the Genelist"""
        if (len(entry) > 2 and isinstance(entry, str)):
            self.list.append(entry)
        else:
            raise TypeError("Could not add whatever you gave me:\n" + entry)
    
    def reset(self, relist):
        """Replaces the current Genelist list member with the supplied relist"""
        self.list = relist
        
    def clear(self):
        """Sets the list member to an empty list"""
        self.list = []
        
        
class FastaEntry:
    """An entry in a FastaFile, header and sequence"""
    
    def __init__(self, header = "", sequence = ""):
        """Default constructor with empty strings for header and sequence"""
        self.head = header
        self.seq = sequence.upper()
        
    def getSequence(self):
        """Returns the sequence member"""
        return self.seq
    
    def getHeader(self):
        """Returns the header member"""
        return self.head
    
    def setHeader(self, header):
        """Set the header member"""
        self.head = header
        
    def setSequence(self, sequence):
        """Set the sequence member"""
        self.seq = sequence.upper()
        
    def getiSequence(self, start, stop):
        """Get a fragment of the sequence, from start to stop"""
        return self.seq[start:stop]
    
    def getSequenceWindow(self, center, windowsize):
        """Get a fragment of the sequence intelligently, guaranteed to be of length windowsize,
        approximately centered around center"""
        if len(self.seq) < windowsize:
            raise ValueError("Sequence is shorter than the requested window size")
        begin = max(0, center - windowsize / 2)
        end = min(self.getSequenceLength(), center + windowsize / 2)
        if end - begin < windowsize:
            if begin == 0:
                end = windowsize
            else:
                begin = end - windowsize
        return self.seq[begin:end]
        
    def getFasta(self):
        """Returns a properly formatted FASTA entry with > and one newline"""
        return (">" + self.head + "\n" + self.seq)
    
    def setFasta(self, fasta):
        """Parse the supplied fasta string for a header and use the rest as sequence"""
        self.head = self.parseForHeader(fasta)
        self.seq = self.parseForSequence(fasta).upper()
    
    def getSequenceLength(self):
        """Return lenght of the sequence member"""
        return len(self.seq)
        
    def parseForHeader(self, fasta):
        """Parse a string for a typical >.\n fasta header"""
        if (not isinstance(fasta, str)):
            head = fasta[0].strip()
        else:
            head = fasta.strip()
        if (head.find("\n") != -1):
            head = head[:head.find("\n")]
        if (head[0] == '>'):
            return head[1:]
        else:
            raise TypeError("Not a FastA Sequence:\n" + fasta)

    def parseForSequence(self, fasta):
        """Parse a string for a typical fasta sequence behind a header"""
        if (not isinstance(fasta, str)):
            if (fasta[0][0] != ">"):
                raise TypeError("Not a FastA Sequence:\n" + fasta)
            seqs = fasta[1:]
            seq = ""
            for s in seqs:
                seq = seq + s.strip()
            return seq
        else:
            if (fasta[0] != ">"):
                raise TypeError("Not a FastA Sequence:\n" + fasta)
            seq = fasta.strip()
            seq = seq[seq.find("\n"):]
            seq.replace("\n", "")
            return seq   
    
    def nucleotideReplace(self, replacement_character):
        """Replace non-nucleotide characters in self.seq with
        one of A, C, G, T or R (i.e., randomly)"""
        assert len(replacement_character) == 1
        assert isinstance(replacement_character, str)
        if replacement_character not in "ACGTR":
            raise ValueError("replacement_character has to be one of A, C, G, T, or R (random)")
        if replacement_character == "R":
            replacement_character = "ACGT"[randint(0,3)]
        seq2 = sub("[^ACGT]", replacement_character, self.seq.upper())
        assert len(self.seq) == len(seq2)
        self.seq = seq2
 
    
class FastaFile:
    """A class containing several FastaEntries"""
    def __init__(self, fastafile = "", directory = ""):
        """Default constructor, optionally supply a fastafile to be read"""
        self.fasta = []
        if fastafile != "":
            self.read(fastafile, directory)
        self.filename = os.path.join(directory, fastafile)
        self.unique = False
        
    def __len__(self):
        """Used for len(self), returns the number of FastaEntries"""
        return len(self.fasta)
    
    def nucleotideReplace(self, replacement_character):
        """Replace non-nucleotide characters in the FastaEntries with replacement_character"""
        for fastaentry in self.fasta:
            fastaentry.nucleotideReplace(replacement_character)
    
    def getUniqueFasta(self):
        """Returns a uniqified version of this FastaFile, i.e. headers are unique"""
        ret = FastaFile()
        ret.unique = True
        ret.fasta = self.getSubset(set(self.getGenes())).fasta
        return ret
    
    def isUnique(self):
        """Returns the unique member, i.e. headers are unique"""
        return self.unique
    
    def getFasta(self):
        """Returns the list of FastaEntries"""
        return self.fasta
    
    def getiFasta(self, i):
        """Returns the ith FastaEntry"""
        return self.fasta[i]
    
    def getLength(self):
        """Returns the number of FastaEntries"""
        return len(self)
    
    def getFormattedFasta(self, filename):
        """Returns a string repesenting a regular FASTA formatted file"""
        fastafile = open(filename, 'w')
        for f in self.fasta:
            fastafile.write(f.getFasta() + "\n")
        fastafile.close()
    
    def getGenes(self):
        """Returns the FASTA headers"""
        ret = []
        for fa in self.fasta:
            ret.append(fa.getHeader())
        return ret
       
    def getFastaEntry(self, header):
        """Returns the FastaEntry object with header = header"""
        for entry in self.fasta:
            if entry.getHeader().startswith(header):
                return entry
        raise IndexError("No entry named " + header)
    
    def getiSequence(self, i):
        """Returns the sequence as a string of the ith FastaEntry"""
        return self.fasta[i].getSequence()
    
    def getAllSequences(self):
        """Returns all sequences in a list"""
        seqlist = []
        for fastaEntry in self.fasta:
            seqlist.append(fastaEntry.getSequence())
        return seqlist
                
    def getSequenceWindow(self, header, center, windowsize):
        """Returns a window cut-out from the sequence, centered approximately around center,
        of length windowsize. If windowsize is larger than sequence, an error occurs."""
        return self.getFastaEntry(header).getSequenceWindow(center, windowsize)
            
    def addFastaEntry(self, fastaentry):
        """Append a FastaEntry object"""
        assert isinstance(fastaentry, FastaEntry)
        self.fasta.append(fastaentry)
        
    def addFasta(self, fastastr):
        """Add a string as a FastaEntry by parsing it for a header and a sequence"""
        fastaentry = FastaEntry()
        self.fasta.append(fastaentry.setFasta(fastastr))
            
    def read(self, fastafile, directory = '', unique = True):    
        """Read a FASTA formatted file and add its entries to this FastaFile object"""
        lfile = open(os.path.join(directory, fastafile))
        self.filename = os.path.join(directory, fastafile)
        parsee = lfile.readlines()
        #parsee.append("END")
        fe = FastaEntry()
        if unique:
            headers = {}
            self.unique = True
        for i in xrange(len(parsee)):
            if (len(parsee[i].strip()) < 2):
                continue
            if (parsee[i].startswith(">")): #new record
                if i > 0: 
                    self.addFastaEntry(fe)
                    fe = FastaEntry()
                fe.setHeader(parsee[i][1:].strip())
                if unique:
                    thisheader = fe.getHeader()
                    if thisheader.find(" ") > 1:
                        #this truncates the header at a space for the motif finder
                        thisheader = thisheader[:thisheader.find(" ")]
                        fe.setHeader(thisheader)
                    if thisheader.find("\t") > 1:
                        #this truncates the header at a tab for the motif finder
                        thisheader = thisheader[:thisheader.find("\t")]
                        fe.setHeader(thisheader)
                    if thisheader in headers:
                        fe.setHeader(thisheader + "_" + str(headers[thisheader]))
                        headers[thisheader] += 1
                    else:
                        headers[thisheader] = 2
            else: #continue with same record
                fe.setSequence(fe.getSequence() + parsee[i].strip())
        assert fe not in self.fasta
        self.addFastaEntry(fe)

    def getSubset(self, genelist):
        """returns a subset of the FastaFile that contains only the genes passed in list"""
        ret = FastaFile()
        if isinstance(list, Genelist):
            genelist = genelist.getList()
        for item in genelist:
            entry = self.getFastaEntry(item)
            if entry != 0:
                ret.addFastaEntry(entry)
            else:
                ### error msg
                print "Unmatched entry: " + item
        return ret
             

class GffFile:
    """wraps several GffEntry instances to one file"""
    def __init__(self):
        self.gff = []
    
    def getGff(self):
        """Returns the list of GffEntries"""
        return self.gff
    
    def getiGff(self, i):
        """Returns the ith GffEntry"""
        return self.gff[i]
    
    def __len__(self):
        """Used for the len() command, returns the number of GffEntries"""
        return len(self.gff)        
    
    def getLength(self):
        """Returns the number of GffEntries"""
        return len(self)
        
    def getFormattedGff(self):
        gffstring = "#INCLUSive GFF File generated by KIRMES.web\n"
        for line in self.gff:
            gffstring = gffstring + line.getGffLine()
        return gffstring
    
    def read(self, gffFile, dir = ""):    
        if (isinstance(gffFile, file)):
            parsee = gffFile.readlines()
        else:
            lfile = open(os.path.join(dir, gffFile))
            #lfile.readline() #skip the first line
            parsee = lfile.readlines()
        for line in parsee:
            if line.startswith("#"):
                continue
            ge = GffEntry()
            ge.setFromLine(line)
            self.add(ge)
    
    def add(self, ge):
        self.gff.append(ge)
        
    def getHitsPerGene(self):
        #score = 0
        genedict = {}
        for ge in self.gff:
            if (genedict.has_key(ge.getGene())):
                genedict[ge.getGene()] = genedict[ge.getGene()] + 1
            else:
                genedict[ge.getGene()] = 0
        return genedict
    
    def getBestHit(self):
        max = 0
        bestline = ""
        for line in self.gff:
            if max < line.getScore():
                max = line.getScore()
                bestline = line
        return bestline
    
    def getGenesBestHit(self, gene):
        max = 0
        bestline = ""
        for line in self.gff:
            if line.getGene() == gene:
                if max < line.getScore():
                    max = line.getScore()
                    bestline = line
        return bestline
    
    def getAllBestHits(self):
        best = {}
        for line in self.gff:
            if best.has_key(line.getGene()):
                max = best[line.getGene()]
                if max.getScore() < line.getScore():
                    best[line.getGene()] = line
            else:
                best[line.getGene()] = line
        return best
    
    def getAllMotifsAllGenes(self):
        genedict = {}
        for ge in self.gff:
            if ge.getMatrixFile() not in genedict:
                genedict[ge.getMatrixFile()] = {}
            if ge.getGene() not in genedict[ge.getMatrixFile()]:
                genedict[ge.getMatrixFile()][ge.getGene()] = ge
            elif genedict[ge.getMatrixFile()][ge.getGene()].getScore() > ge.getScore():
                genedict[ge.getMatrixFile()][ge.getGene()] = ge
        motifs = {}
        for k,v in genedict.iteritems():
            motifs[k] = v.values()
        return motifs
    
    def getGenes(self):
        ret = []
        for i in xrange(len(self.gff)):
            ret.append(self.gff[i].getGene())
        return ret

class GffEntry:
    """Representing one line in a gff File"""
    def __len__(self):
        return abs(self.start - self.stop)
    
    def __init__(self):
        self.gene = ""
        self.program = "" 
        self.feature = ""
        self.start = 0
        self.stop = 0
        self.length = 0
        self.score = 0.0
        self.strand = True
        self.dot = "."
        self.id = ""
        self.site = ""
        self.matrixFile = ""
    
    def setFromLine(self, line):
        array = line.split("\t")
        self.setGene(array[0])
        self.setProgram(array[1])
        self.setFeature(array[2])
        self.setStartStop(array[3], array[4])
        self.setScore(array[5])
        self.setStrand(array[6])
        self.setDot(array[7])
        self.setId(array[8])
        self.siteFromId()
        
    def setGene(self, line):
        self.gene = line.strip()
        
    def setProgram(self, line):
        self.program = line.strip()
        
    def setFeature(self, line):
        self.feature = line.strip()
        
    def setStartStop(self, s, e):
        self.start = int(s.strip())
        self.stop = int(e.strip())
        self.length = self.start - self.stop
        if (self.length < 0):
            self.length = -self.length
    
    def setScore(self, line):
        self.score = line.strip()    
        
    def setStrand(self, line):
        if (line.strip() == "+"):
            self.setStrandFromBool(True)
        elif (line.strip() == "-"):
            self.setStrandFromBool(False)
    
    def setStrandFromBool(self, bool):
        self.strand = bool
        
    def setDot(self, line):
        self.dot = line.strip()
        
    def setId(self, line):
        self.id = line.strip()
        
    def siteFromId(self): #this writes the motif found at this position in a separate var
        pos = self.id.find("site")    
        self.site = self.id[pos + 6 : -2] # subtracts 'site "' and '";'
        self.matrixFile = self.id[4 : pos - 3] # subtracts the additional chars around the matrix file name
 
    def strandString(self):
        if (self.strand):
            return "+"
        else:
            return "-"
        
    def getGffLine(self):
        return (self.gene + '\t' + self.program + '\t' 
                + self.feature + '\t' + str(self.start) + '\t' 
                + str(self.stop) + '\t' + str(self.score) + '\t' 
                + self.strandString() + '\t' + self.dot 
                + '\t' + self.id + '\n')
        
    def getGene(self):
        return self.gene
    
    def getProgram(self):
        return self.program
    
    def getFeature(self):
        return self.feature
    
    def getStart(self):
        return self.start
    
    def getStop(self):
        return self.stop
    
    def getLength(self):
        return self.length
    
    def getScore(self):
        return self.score
    
    def getStrand(self):
        return self.strand
    
    def getStrandString(self):
        return self.strandString()
    
    def getDot(self):
        return self.dot
    
    def getId(self):
        return self.id
    
    def getSite(self):
        return self.site
        
    def getMatrixFile(self):
        return self.matrixFile
    