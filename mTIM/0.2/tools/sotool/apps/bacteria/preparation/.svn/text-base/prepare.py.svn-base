#!/usr/bin/env python2.6
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Written (W) 2011 Christian Widmer
#             2011 Nico Goernitz
# Copyright (C) 2011 Max-Planck-Society

"""
Created on 05.10.2011
@author: Christian Widmer, Nico Goernitz
@summary: Create structured output training set for prokaryotic genome

"""

import os
import numpy
import scipy.io
import re

import taxonomy

class PTTEntry(object):
    """
    simple container for PTT file entry
    """

    
    def __init__(self, row):
        """
        parse from line
        """
        
        tokens = row.strip().split("\t")

        self.location = tokens[0]
        self.strand = tokens[1]
        self.length = int(tokens[2])
        self.pid = tokens[3]
        self.gene = tokens[4] 
        self.synonym = tokens[5]
        self.code = tokens[6]
        self.cog = tokens[7]
        if self.cog.startswith("COG"):

            # eg: AFGT1234M --> 1234
            p = re.compile("\D+(\d+)\D*")
            cog, = p.findall(self.cog)
            self.cog = int(cog)

        else:
            self.cog = -1

        self.product = tokens[8]

        self.start = int(self.location.split("..")[0])

        # disregard stop codon
        self.stop = int(self.location.split("..")[1]) - 3


    def __repr__(self):
        """
        pretty print
        """
        
        return "%s %s" % (self.synonym, self.cog)


def parse_ptt_file(fn):
    """
    read in ptt file and construct relevant data structure

    Location    Strand  Length  PID     Gene    Synonym     Code    COG     Product
    0           1       2       3       4       5           6       7       8
    """

    cogs = set()
    entries = []

    for line_num, line in enumerate(file(fn)):

        # skip header
        if line_num <= 2:
            continue

        # add current line
        entry = PTTEntry(line)
        entries.append(entry)

        # sanity check
        if entry.cog in cogs:
            print "paralog found:", entry.cog

        if entry.cog != -1:
            cogs.add(entry.cog)
        

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
        
    
        self.alphabet = ["A", "C", "G", "T"]
        self.valid_start_codons = ["ATG", "GTG", "TTG", "ATT", "CTG"]
        self.valid_stop_codons = ["TAG", "TAA", "TGA"]

        self.codons = []
        self.codon_to_idx = {}
        self.idx_to_codon = ["XXX"]*(len(self.alphabet)**3)
        
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
       
        # get idx of start and stop codons
        self.valid_start_codons_idx = [self.codon_to_idx[codon] for codon in self.valid_start_codons]
        self.valid_stop_codons_idx = [self.codon_to_idx[codon] for codon in self.valid_stop_codons]
 
         
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
    #codons_idx = numpy.array([cm.get_idx(codon) for codon in codons])
    codons_idx = ([cm.get_idx(codon) for codon in codons])

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
    
    start_codon = str(codons[offset_left]) 
    print "start codon", start_codon, "stop codon", stop_codon, not str(codons[offset_left]) in cm.valid_start_codons
    
    # create label
    #  mask = numpy.array([0.0]*offset_left + [1.0]*(length_cds+1) + [0.0]*(offset_right-1)) # include stop codon
    mask = ([0]*offset_left + [1]*(length_cds+1) + [0]*(offset_right-1)) # include stop codon
    interval = (offset_left, offset_left + length_cds + 1) 


    #TODO: sicherstellen dass kein stopcodon inframe vorkommt
    assert (len(mask) == length)
    assert (len(mask) == len(codons) + 2)
    assert (len(mask) == len(codons_idx) + 2)
    assert (interval[1] - interval[0] == length_cds + 1)


    return start_codon, stop_codon, seq_window, codons, codons_idx, mask, interval

     
   

class DataHandler(object):
    """
    class to create list of all relevant data
    """

    def __init__(self, tax):
        """
        set up lists
        """
        
        # setup self.return value
        self.ret = {}
        self.ret["sequence"] = []
        self.ret["codons"] = []
        self.ret["codons_idx"] = []
        self.ret["mask"] = []
        self.ret["relative_interval"] = []
        self.ret["interval"] = []
        self.ret["example_id"] = []
        self.ret["strand"] = []
        self.ret["cog"] = []
        self.ret["org_id"] = []
        self.ret["org_interval"] = []
        self.ret["start_codon"] = []
        self.ret["stop_codon"] = []

        # taxonomy information
        self.ret["org_names"] = ""
        self.ret["taxonomy"] = []

        # codon maps
        cm = CodonMapper()
        self.ret["idx_to_codon"] = cm.idx_to_codon
        self.ret["start_codons_idx"] = cm.valid_start_codons_idx
        self.ret["stop_codons_idx"] = cm.valid_stop_codons_idx

        print self.ret

        # set up mapping task_to_id
        self.task_to_id = {}
        for (task_id, leaf) in enumerate(tax.get_leaves()):
            assert not self.task_to_id.has_key(leaf.name), "duplicate key %s" % (leaf.name)
            self.task_to_id[leaf.name] = task_id

        # set up mapping node_to_id and concat node names
        self.node_to_id = {}
        for (node_id, node) in enumerate(tax.get_all_nodes()):
            assert not self.node_to_id.has_key(node.name), "duplicate key %s" % (node.name)
            self.node_to_id[node.name] = node_id
            self.ret["org_names"] += node.name + "\t"

        # process tax
        self.taxonomy = tax
        self.add_taxonomy(tax)




    def add_organism(self, org_id, fna_fn, ptt_fn):
        """
        fetch relevant information from genome based on ptt file
        """


        # load ptt
        entries = parse_ptt_file(ptt_fn)

        # split by strand
        coding_seqs_pos = [entry for entry in entries if entry.strand == "+"]
        coding_seqs_neg = [entry for entry in entries if entry.strand == "-"]
        
        print "num of CDS (only on the + strand):", len(coding_seqs_pos)
        print "num of CDS (only on the - strand):", len(coding_seqs_neg)

        # load genome
        genome_handler = GenomeHandlerSingleContig(fna_fn)
        

        num_skipped = 0
        counter = 0


        # keep track of the intervals for the examples
        # starting with 1 (matlab conform)
        current_pos = 1;

        # 
        print "warning: using only positive strand"
        entries = coding_seqs_pos

        # skip first and last
        #for (idx, tmp_entry) in enumerate(coding_seqs_pos):
        for (idx, tmp_entry) in enumerate(entries):
            print "processing entry:", idx
            
            # fix corner cases
            if idx == 0:
                previous_end = 0
            else:
                previous_end = entries[idx-1].stop
            if idx == len(entries) - 1:
                following_start = genome_handler.get_length() 
            else:
                following_start = entries[idx+1].start
            
            # respect distance between genes
            offset_left = (entries[idx].start - previous_end) / 2
            offset_right = (following_start - entries[idx].stop) / 2
           
            # skip if offset is too small
            if offset_left <= 3 or offset_right <= 3:
                print "skipping example %i, offset_left=%i, offset_right=%i" % (idx, offset_left, offset_right)
                num_skipped += 1
                #import ipdb
                #ipdb.set_trace()
                continue
            
            try:
                start_codon, stop_codon, seq_window, codons, codons_idx, mask, interval = create_example(tmp_entry, genome_handler, offset_left, offset_right)
            except Exception:
                print "problem encountered, skipping example %i" % (idx)
                num_skipped += 1
                continue
            
            seq_length = len(seq_window)    # the length of the current sequence we need for tracking the 
                                            # intervals for the examples

            self.ret["sequence"] += seq_window
            self.ret["codons"] += codons + ['XXX', 'YYY']    # add filling codons such that length of codons(_idx) is the same
            self.ret["codons_idx"] += codons_idx + [-1,-1]   # as sequence length
            self.ret["mask"] += mask
            self.ret["relative_interval"] += [interval]    # the position of the cds gene within the current region (example)
            self.ret["interval"] += [(current_pos,current_pos+seq_length-1)]  # all nukleotides within this region belong to 
                                                                         # the current examples (inlcuding current_pos and
                                                                         # current_pos+seq_length-1)
            self.ret["example_id"] += [counter]
            self.ret["strand"] += [tmp_entry.strand]
            #import ipdb
            #ipdb.set_trace()
            self.ret["start_codon"] += [start_codon]
            self.ret["stop_codon"] += [stop_codon]
            self.ret["cog"] += [tmp_entry.cog]
            self.ret["org_id"] += [org_id]

            counter += 1

            #        if counter >= 50:
            #            break

            current_pos += seq_length    # should point to the first position after the last entry
         

        # encode org_interval as tuple
        if len(self.ret["org_interval"]) > 0:
            left = self.ret["org_interval"][-1][1]
        else:
            left = 0

        self.ret["org_interval"].append( (left, counter) ) 
        print "number of skipped examples %i, remaining examples %i" % (num_skipped, counter)
        

    def save_to_file(self, save_path):
        """
        save output
        """

        # save to matlab file
        mat_fn = "%scds_genes.mat" % save_path

        #print self.ret["codons_idx"]
        scipy.io.savemat(mat_fn, self.ret)
        
        return self.ret


    def add_taxonomy(self, tax):
        """
        hierarchical training procedure
        """


        # enqueue root node 
        grey_nodes = [tax]

        # set up matrix encoding
        #TODO think about safe mapping
        num_leaves = len(tax.get_leaves())
        num_nodes = len(tax.get_all_nodes())
        matrix_encoding = numpy.zeros( (num_nodes, num_leaves+1) )

        print "matrix.shape", matrix_encoding.shape

        
        #####################################################
        #         top-down processing of taxonomy           #
        #####################################################
     
        while len(grey_nodes)>0:
        
            node = grey_nodes.pop(0) # pop first item

            # set up name mapping
            node_id = self.node_to_id[node.name]
        
            # enqueue children
            if node.children != None:
                grey_nodes.extend(node.children)
            
            # get data below current node
            data_keys = node.get_data_keys()
            
            # fill row
            for key in data_keys:
                task_id = self.task_to_id[key]
                matrix_encoding[node_id, task_id] = 1

            # store parent information in last column
            if node.is_root():
                # no parent at root node
                matrix_encoding[node_id, -1] = -1
            
            else:
                # regularize against parent predictor
                matrix_encoding[node_id, -1] = self.node_to_id[node.parent.name]

    
        print matrix_encoding
        self.ret["taxonomy"] = matrix_encoding

        return matrix_encoding




def main():
    """
    generate dataset
    """

    data_path = "../../../data/"

    # set up container for data
    tax = taxonomy.create_tax_four()
    #tax = taxonomy.create_tax_many()
    dh = DataHandler(tax)


    for org in tax.get_data_keys():

        org_id = dh.task_to_id[org] 
        print "processing organism %s (id=%i)" % (org, org_id)
        
        org_path = data_path + org + "/"
        save_path = org_path
        
        ptt_files = [ptt for ptt in os.listdir(org_path) if ptt.endswith(".ptt")]
        print("Gff files: %s",ptt_files)
        ptt_file_sizes = [os.stat(org_path + ptt).st_size  for ptt in ptt_files] 
    
        print ptt_file_sizes
    
        # pick largest one (in case we have several contigs)
        largest_idx = 0;
        if (len(ptt_file_sizes)>1):
            print 'There are multiple .ptt files available for this organism.'
            print 'I only use the biggest one.'
            largest_idx = numpy.argmax(ptt_file_sizes)
            print 'Largest file is: ', largest_idx
        
        # determine file names
        selected_ptt = org_path + ptt_files[largest_idx]
        selected_fna = selected_ptt.replace(".ptt", ".fna")
       
        # invoke generation procedure
        dh.add_organism(org_id, selected_fna, selected_ptt)
    
        import ipdb
        ipdb.set_trace()
  
    # store final result
    save_path = "/tmp/"
    dh.save_to_file(save_path)
   


if __name__ == '__main__':
    main()
    
