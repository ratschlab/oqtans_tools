#!/usr/bin/env python2.6
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Written (W) 2011 Christian Widmer
# Copyright (C) 2011 Max-Planck-Society

"""
Created on 23.05.2011
@author: Christian Widmer
@summary: Create structured output training set for prokaryotic genome

"""

import os
import numpy
import scipy.io


class GFFEntry(object):
    """
    simple container for GFF file entry
    """

    
    def __init__(self, row):
        """
        parse from line
        """
        
        tokens = row.strip().split("\t")

        self.contig = tokens[0]
        self.event = tokens[1]
        self.type = tokens[2]
        self.start = int(tokens[3])
        self.stop = int(tokens[4]) 
        self.strand = tokens[6]
        self.id = tokens[8]


    def __repr__(self):
        """
        pretty print
        """
        
        return "%s %s %s %i %i %s %s" % (self.contig, self.event, self.type, self.start, self.stop, self.strand, self.id)



def parse_gff_file(input_fn):
    """
    parse gff file with data description
    """

    print "loading", input_fn

    data_file = file(input_fn)
    
    #genome_handler = GenomeHandler()
    
    entries = []
    
    for line in data_file:
        if not line.startswith("#"):
            entries.append(GFFEntry(line))
    
    
    print "parsed %i entries" % (len(entries))
    
    return entries



class GenomeHandlerSingleContig:
    """
    class to load sequences for one-contig genomes of prokaryotes
    """
    
    def __init__(self, fasta_fn):
        
        self.fasta_fn = fasta_fn
        self.seqs = {}

        from Bio import SeqIO
        seq_io = SeqIO.parse(file(fasta_fn), "fasta")
        
        for record in seq_io:

                print "loading chromosome %s" % (record.id)
                
                self.seq = record
                

    def get_codon(self, pos):
        
        return self.seq.seq[pos:pos + 3]
    
    
    def get_window(self, start, stop):
                       
        return self.seq.seq[start:stop]
    
    
    def get_length(self):
        
        return len(self.seq)


class CodonMapper(object):
    """
    helper to take care of two way mapping between codons and idx
    """

    def __init__(self):
        """
        small helper to construct codon map
        """
        
        self.codons = []
        self.codon_to_idx = {}
        self.idx_to_codon = {}
    
        self.alphabet = ["A", "C", "G", "T"]
        self.valid_start_codons = ["ATG", "GTG", "TTG", "ATT", "CTG"]
        self.valid_stop_codons = ["TAG", "TAA", "TGA"]
        
        counter = 0

        
        for a1 in self.alphabet:
            for a2 in self.alphabet:
                for a3 in self.alphabet:
    
                    # construct codon and create two-way mapping
                    codon = a1 + a2 + a3
                    self.codons.append(codon)
                    self.codon_to_idx[codon] = counter
                    self.idx_to_codon[counter] = codon
                    
                    counter += 1
        
        # put start codons to front in integer encoding
        for (idx, start_codon) in enumerate(self.valid_start_codons):
            start_idx = self.codon_to_idx[start_codon]
            if start_idx != idx:
                tmp_codon = self.idx_to_codon[idx]
                self.codon_to_idx[start_codon] = idx
                self.codon_to_idx[tmp_codon] = start_idx
                
                self.idx_to_codon[idx] = start_codon
                self.idx_to_codon[start_idx] = tmp_codon
        
        
        
        # put stop codons to back in integer encoding
        for (idx, stop_codon) in enumerate(self.valid_stop_codons):
            
            new_idx = 64 - idx - 1  
            stop_idx = self.codon_to_idx[stop_codon]
            if stop_idx != new_idx:
                tmp_codon = self.idx_to_codon[new_idx]
                self.codon_to_idx[stop_codon] = new_idx
                self.codon_to_idx[tmp_codon] = stop_idx
                
                self.idx_to_codon[new_idx] = stop_codon
                self.idx_to_codon[stop_idx] = tmp_codon
        
        
         
    def get_codon(self, idx):
        """
        get codon with index
        """
        
        return self.idx_to_codon[idx]
    
    
    def get_idx(self, codon):
        """
        get idx of condon
        """
        
        return self.codon_to_idx[str(codon)]
        


def create_example(gff_entry, genome_handler, offset_left, offset_right):
    """
    takes coordinate from gff, fetches sequence
    from genome handler
    
    does not do much sanity checking ATM (overlapping genes and so forth)
    
    Start codon:
    In addition to AUG, alternative start codons, mainly GUG and UUG 
    are used in prokaryotes. For example E. coli uses 83% ATG 
    (AUG), 14% GTG (GUG), 3% TTG (UUG) and one or two others (e.g., ATT and CTG).
    
    Stop codon:
    TAG
    TAA
    TGA
    """


    # GFF is one based --> zero based
    start = gff_entry.start - offset_left -1
    stop = gff_entry.stop + offset_right -1
    length = stop - start
    length_cds = gff_entry.stop - gff_entry.start 
    
    
    
    assert(start >= 0)
    assert(stop < genome_handler.get_length())
    assert(length > 0)
    
    cm = CodonMapper()
    
    seq_window = str(genome_handler.get_window(start, stop))
    codons = [str(genome_handler.get_codon(pos)) for pos in xrange(start, stop-2)]
    codons_idx = numpy.array([cm.get_idx(codon) for codon in codons])

    # prokaryotes use alternative start codons

    
    # sanity checks
    if not str(codons[offset_left]) in cm.valid_start_codons:
        print "BAD START CODON:", codons[offset_left]
        print seq_window
        print "length", length
        print "valid codons", cm.valid_start_codons
        assert False, "bad start codon: %s" % (str(codons[offset_left])) 
    
    # stop codon not included in CDS
    stop_codon = str(codons[offset_left + length_cds + 1])
    if not stop_codon in cm.valid_stop_codons:
        print "BAD STOP CODON:", stop_codon
        print seq_window
        print "length", length
        print "valid codons", cm.valid_stop_codons
        assert False, "bad stop codon"
    
    
    print "start codon", str(codons[offset_left]), "stop codon", stop_codon, not str(codons[offset_left]) in cm.valid_start_codons
    
    # create label
    mask = numpy.array([0.0]*offset_left + [1.0]*(length_cds+1) + [0.0]*(offset_right-1)) # include stop codon
    interval = (offset_left, offset_left + length_cds + 1) 


    #TODO: sicherstellen dass kein stopcodon inframe vorkommt
    assert (len(mask) == length)
    assert (len(mask) == len(codons) + 2)
    assert (len(mask) == len(codons_idx) + 2)
    assert (interval[1] - interval[0] == length_cds + 1)


    return seq_window, codons, codons_idx, mask, interval



def create_dataset_wg(gff_entries, genome_handler, save_path):
    """
    fetch sequences from genome based on CDS from gff 
    and save to datastructure that containes whole genome
    and mask
    """
    

    feat = numpy.zeros((3, genome_handler.get_length()))

    print "shape", feat.shape

    # skip first and last
    for (idx, gff) in enumerate(gff_entries):
        print "processing entry:", idx

        print gff.start. gff.stop
        
    
    return feat


 def create_dataset(gff_entries, genome_handler, save_path):
    """
    fetch sequences from genome based on CDS from gff 
    and saves one matlab file per example
    """
    

    num_skipped = 0
    counter = 0

    # skip first and last
    for (idx, tmp_entry) in enumerate(gff_entries):
        print "processing entry:", idx
        
        # fix corner cases
        if idx == 0:
            previous_end = 0
        else:
            previous_end = gff_entries[idx-1].stop
        if idx == len(gff_entries) - 1:
            following_start = genome_handler.get_length() 
        else:
            following_start = gff_entries[idx+1].start
        
        # respect distance between genes
        offset_left = (gff_entries[idx].start - previous_end) / 2
        offset_right = (following_start - gff_entries[idx].stop) / 2
        
        # set max offset to 500
        #offset_left = min(offset_left, 500)
        #offset_right = min(offset_right, 500)
        
        # skip if offset is too small
        if offset_left <= 3 or offset_right <= 3:
            print "skipping example %i, offset_left=%i, offset_right=%i" % (idx, offset_left, offset_right)
            num_skipped += 1
            continue
        
        try:
            seq_window, codons, codons_idx, mask, interval = create_example(tmp_entry, genome_handler, offset_left, offset_right)
        except Exception:
            print "problem encountered, skipping example %i" % (idx)
            num_skipped += 1
            continue
        
        ret = {}
        ret["sequence"] = seq_window
        ret["codons"] = codons
        ret["codons_idx"] = codons_idx
        ret["mask"] = mask
        ret["interval"] = interval
        ret["example_id"] = counter

        mat_fn = "%sexample_%i.mat" % (save_path, counter)

        
        assert codons_idx[interval[0]] <= 4, "start codon idx %i %s %s %i %i %s" % (codons_idx[interval[0]], codons[interval[0]], seq_window, interval[0], interval[1], seq_window[interval[0]:interval[1]])
        assert codons_idx[interval[1]] >= 61, "stop codon idx %i < %i" % (codons_idx[interval[1]], 61)

        # save to matlab format
        scipy.io.savemat(mat_fn, ret)
        
        counter += 1
        
    
    print "number of skipped examples %i, remaining examples %i" % (num_skipped, counter)
    
    
    return ret

     
   
   
def generate_dataset(save_path, fna_fn, gff_fn):
    """
    load gff and genome, invoke generation routine
    """ 
    
    # load gff
    gff_entries = parse_gff_file(gff_fn)
    coding_seqs_pos = [entry for entry in gff_entries if entry.type == "CDS" and entry.strand == "+"]
    #coding_seqs_pos = [entry for entry in gff_entries if entry.type == "CDS" and entry.strand == "+"][1:20] #skip first and last
    
    print "num of CDS:", len(coding_seqs_pos)
    
    # load genome
    genome_handler = GenomeHandlerSingleContig(fna_fn)
    
    # create data structure with relevant information
    data = create_dataset(coding_seqs_pos, genome_handler, save_path)
    
    return data



def main():

    
    #"Deinococcus_radiodurans_R1_uid57665",#Bacteria; Deinococcus-Thermus; Deinococci; Deinococcales; Deinococcaceae; Deinococcus; Deinococcus radiodurans; Deinococcus radiodurans R1
    #"Listeria_monocytogenes_EGD_e_uid61583", #Bacteria; Firmicutes; Bacillales; Listeriaceae; Listeria; Listeria monocytogenes; Listeria monocytogenes serotype 4b str. F2365
    #"Streptomyces_coelicolor_A3_2__uid57801", #Bacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Streptomycineae; Streptomycetaceae; Streptomyces; Streptomyces coelicolor; Streptomyces coelicolor A3(2)
    
    data_path = "/fml/ag-raetsch/home/cwidmer/Documents/phd/projects/so-mtl/sotool/data/"
    data_path = "../../data/"

    org_list = [
        "Escherichia_coli_BW2952_uid59391", #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Escherichia; Escherichia coli; 
        "Escherichia_fergusonii_ATCC_35469_uid59375", #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Escherichia; Escherichia fergusonii; Escherichia fergusonii ATCC 35469
        "Enterobacter_638_uid58727", #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Enterobacter; Enterobacter sp. 638
        "Klebsiella_pneumoniae_MGH_78578_uid57619", #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Klebsiella; Klebsiella pneumoniae; Klebsiella pneumoniae subsp. pneumoniae MGH 78578
        "Salmonella_enterica_serovar_Heidelberg_SL476_uid58973", #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Salmonella; Salmonella enterica; Salmonella enterica subsp. enterica serovar Typhi; Salmonella enterica subsp. enterica serovar Typhi str. CT18
        
        "Agrobacterium_tumefaciens_C58_uid57865", #Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Rhizobiaceae; Rhizobium/Agrobacterium group; Agrobacterium; Agrobacterium tumefaciens; Agrobacterium tumefaciens str. C58
        "Rhizobium_leguminosarum_bv__viciae_3841_uid57955", #Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Rhizobiaceae; Rhizobium/Agrobacterium group; Rhizobium; Rhizobium leguminosarum; Rhizobium leguminosarum bv. viciae 3841
            
        "Helicobacter_pylori_26695_uid57787", #Bacteria; Proteobacteria; Epsilonproteobacteria; Campylobacterales; Helicobacteraceae; Helicobacter; Helicobacter pylori; Helicobacter pylori 26695
        
        
        
        "Mycobacterium_tuberculosis_H37Rv_uid57777", #Bacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Corynebacterineae; Mycobacteriaceae; Mycobacterium; Mycobacterium tuberculosis complex; Mycobacterium tuberculosis; Mycobacterium tuberculosis H37Rv
        "Bifidobacterium_longum_NCC2705_uid57939", #Bacteria; Actinobacteria; Actinobacteridae; Bifidobacteriales; Bifidobacteriaceae; Bifidobacterium; Bifidobacterium longum; Bifidobacterium longum subsp. infantis ATCC 15697
        
        "Bacillus_anthracis_Ames_uid57909", #Bacteria; Firmicutes; Bacillales; Bacillaceae; Bacillus; Bacillus cereus group; Bacillus anthracis; Bacillus anthracis str. 'Ames Ancestor'
        "Bacillus_subtilis_168_uid57675", #Bacteria; Firmicutes; Bacillales; Bacillaceae; Bacillus; Bacillus subtilis; Bacillus subtilis subsp. natto BEST195
        
        "Bacteroides_thetaiotaomicron_VPI_5482_uid62913", #Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; Bacteroides thetaiotaomicron; Bacteroides thetaiotaomicron VPI-5482
        "Flavobacterium_psychrophilum_JIP02_86_uid61627", #Bacteria; Bacteroidetes; Flavobacteria; Flavobacteriales; Flavobacteriaceae; Flavobacterium; Flavobacterium psychrophilum; Flavobacterium psychrophilum JIP02/86
        
        "Mycoplasma_genitalium_G37_uid57707", #Bacteria; Tenericutes; Mollicutes; Mycoplasmataceae; Mycoplasma; Mycoplasma genitalium; Mycoplasma genitalium G37
        "Mycoplasma_pneumoniae_M129_uid57709", #Bacteria; Tenericutes; Mollicutes; Mycoplasmataceae; Mycoplasma; Mycoplasma pneumoniae; Mycoplasma pneumoniae M129
        
        "Clostridium_botulinum_A_ATCC_19397_uid58927", #Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium; Clostridium botulinum; Clostridium botulinum A str. ATCC 3502
        "Streptococcus_thermophilus_CNRZ1066_uid58221" #Bacteria; Firmicutes; Lactobacillales; Streptococcaceae; Streptococcus; Streptococcus thermophilus; Streptococcus thermophilus CNRZ1066
                
        "Acidobacterium_capsulatum_ATCC_51196_uid59127", #Bacteria; Acidobacteria; Acidobacteriales; Acidobacteriaceae; Acidobacterium; Acidobacterium capsulatum; Acidobacterium capsulatum ATCC 51196   

        "Fusobacterium_nucleatum_ATCC_25586_uid57885", #Bacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae; Fusobacterium; Fusobacterium nucleatum; Fusobacterium nucleatum subsp. nucleatum ATCC 25586
    ]

    """
    org_list = [
      ## "Sulfolobus_islandicus_M_14_25_uid58849",
      ## "Methanobrevibacter_smithii_ATCC_35061_uid58827"
    	"Escherichia_coli_BW2952_uid59391"
		]
    """


    for org in org_list:
        
        print "processing organism", org
        
        org_path = data_path + org + "/"
        save_path = org_path + "mat/"
        
        gff_files = [gff for gff in os.listdir(org_path) if gff.endswith(".gff")]
        gff_file_sizes = [os.stat(org_path + gff).st_size  for gff in gff_files] 
    
        print gff_file_sizes
    
        # pick largest one (in case we have several contigs
        largest_idx = numpy.argmax(gff_file_sizes) 
        
        # determine file names
        selected_gff = org_path + gff_files[largest_idx]
        selected_fna = selected_gff.replace(".gff", ".fna")
        
        # create output directory if needed
        if not os.path.exists(save_path):
            print "creating save path", save_path
            os.mkdir(save_path)
    
        # invoke generation procedure
        generate_dataset(save_path, selected_fna, selected_gff)
        
   


if __name__ == '__main__':
    
    main()
    
