#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  Written (W) 2005-2009 Gunnar Raetsch, Gabriele Schweikert, Jonas Behr, Soeren Sonnenburg, Alexander Zien, Georg Zeller, Andre Noll, Cheng Soon Ong, Petra Philips
#  Copyright (C) 2005-2009 Max Planck Society and Fraunhofer Institute FIRST 
#

#TODO define interperter
"""
make_gio builds a genome information object for an input fasta file

usage:

make_gio(fasta_file, gio_dir)
make_gio(fasta_file, gio_dir, info_file)

takes 2 or 3 arguments:
   <fasta_file> is the input file in fasta format
   <gio_dir> is the directory to which the genome information object will be written
   <info_file> is the file name to which some information about the fasta file will be written
"""
# import std packages
import os
import sys
import getopt


def make_gio(in_file_name, gio_path, info_file):
    """
    make_gio(fasta_file, gio_dir, info_file)
    takes 2 or 3 arguments:
    <fasta_file> is the input file in fasta format
    <gio_dir> is the directory to which the genome information object will be written to
    <info_file> is the file name to which some information about the fasta file will be written
    """
    try:
        f_in = file(in_file_name, "r")
    except Exception, msg:
        print msg
        print "cannot open infile '" + in_file_name + "'"
        sys.exit(1)
        
    flat_path = os.path.join(gio_path, "genome")
    try:
	if os.path.exists(flat_path):
		print "directory " + flat_path + " exists already."
	else:
        	os.makedirs(flat_path)

    except Exception, msg:
        print msg
        print "cannot create path '" + flat_path + "'"
        sys.exit(1)

    f_out = None
    contig_list = []
    num_bases = 0
    
    for line in f_in:

        if line.isspace():
            print "warning: wrong format. ignoring empty line in file '" + in_file_name + "'"
            continue

        if line[0].isspace():
            print "wrong format: leading white space in file '" + in_file_name + "'"
            sys.exit(1)
    
        if line.startswith(">"):
            
            if f_out != None:
                f_out.close()

	    contig_list.append(line[1:-1].split()[0])
            out_name = os.path.join(flat_path, contig_list[-1] + ".flat")
            try:
                f_out = file(out_name, "w")
                print "creating file '" + out_name + "'"
            except Exception, msg:
                print msg
                print "cannot open file '" + out_name + "'"
                sys.exit(1)
                
        else:
            try:
                f_out.write(line[0:-1].lower())
            except Exception, msg:
		if f_out != None:
                	print msg
	                print "cannot write to file '" +out_name + "'"
        	        sys.exit(1)
                else:
			print "improper input format. No header in first line"
			sys.exit(1)

            num_bases += len(line)-1

    f_out.close()

    try:
        print "creating file '" + os.path.join(gio_path, "genome.config") + "'"
        f_conf = file(os.path.join(gio_path, "genome.config"),"w")
        f_conf.write("BASEDIR " +  os.path.abspath(gio_path) +"\n\n")
        f_conf.write("CONTIGS " +  str(len(contig_list)) +"\n")
        for c in contig_list:
            f_conf.write(c + "\tgenome/" + c + ".flat\tgenome/" + c + ".dna\n")
        f_conf.write("\nALPHABET acgt\n\n")
        f_conf.write("ESTFILES 0\n\n")
        f_conf.write("CDNAFILES 0\n\n")
        f_conf.write("ANNOTATIONFILES 0\n")
        f_conf.close()
    except Exception, msg:
        print msg
        print "cannot create file '" + os.path.join(gio_path, "genome.config") + "'"
        sys.exit(1)

    if info_file == None:
        sys.exit(0)

    try:
        print "creating file '" + info_file + "'"
        f_inf = file(info_file, "a")
        f_inf.write("Genome properties:\n")
        f_inf.write(" * " + str(len(contig_list)) + " contigs\n")
        f_inf.write(" * " + str(int(round(num_bases/1000))) + "kbs total length\n\n")
        f_inf.write("Contig list:\n")
        for c in contig_list:
            f_inf.write("  " + c + "\n")
        f_inf.close()
    except Exception, msg:
        print msg
        print "cannot create file '" + info_file + "'"
        sys.exit(1)
        
def main():

    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs=", ["help", "sizes="])

    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)

    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)

        if o in ("-s", "--sizes"):
            sizes = a.split(",")


    if len(args) < 2:
        print "need at least two arguments"
        sys.exit(2)

    if len(args) > 3:
        print "need at most three arguments"
        sys.exit(2)    

    if len(args) < 3:
        make_gio(args[0], args[1], None)
    else:
        make_gio(args[0], args[1], args[2])
        
if __name__ == "__main__":

    main()

