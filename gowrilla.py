#!/bin/python
"""
#######################################################################################
#                                                                                     #
#    gowrilla.py is a wrapper for the http://cbl-gorilla.cs.technion.ac.il/           #
#    website submission form to create an ontology from a ranked list                 #
#    of proteins or genes. For GOrilla, please cite:                                  #
#    Eden E, Navon R, Steinfeld I, Lipson D, Yakhini Z (2009) GOrilla: A tool         #
#    for discovery and visualization of enriched GO terms in ranked gene lists.       #
#    BMC Bioinformatics 10:48.                                                        #
#                                                                                     #
#    Copyright (C) 2011 Sebastian J. Schultheiss <sebi@umich.edu>                     #
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
#  Original Author: Sebastian J. Schultheiss, version 0.1.0                           #
#                                                                                     #
#  Please add a notice of any modifications here:                                     #
#                                                                                     #
#                                                                                     #
#######################################################################################
"""

__version__ = "0.1.0"
__license__  = "GNU General Public License"
__author__ = "Sebastian J. Schultheiss <sebi@umich.edu>"
__usage__ = """Usage: 
    Supply one ranked list of genes or an unranked and a background list,
    give information about the organism, and you will get an HTML result page
    from a run of the GOrilla web service."""

from urllib import urlopen, urlencode
import sgmllib 
import os, sys, shutil, os.path, math
from optparse import OptionValueError, OptionParser
from spider import Spider
from math import log
from time import sleep

class ResultParser(sgmllib.SGMLParser):
    """Parsing a GOrilla result document, retaining the resulting image"""

    def parse(self, result):
        """Parse the result file"""
        self.feed(s)
        self.close()

    def __init__(self, verbose = 0):
        """Initialize the SGMLParser Superclass"""
        sgmllib.SGMLParser.__init__(self, verbose)
        self.links = []
    


organisms = "ARABIDOPSIS_THALIANA, SACCHAROMYCES_CEREVISIAE, CAENORHABDITIS_ELEGANS, DROSOPHILA_MELANOGASTER, DANIO_RERIO, HOMO_SAPIENS, MUS_MUSCULUS, RATTUS_NORVEGICUS"
ontologies = "proc, func, comp, all"
gorillaurl = "http://cbl-gorilla.cs.technion.ac.il/servlet/GOrilla"
galaxydatasetindexhtml = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE html    PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta HTTP-EQUIV="REFRESH" content="0; url=index.html">
</head>
<body>
</body>
</html>
"""

def check_organism(option, opt_str, value, _parser):
    """Checks if the organism name is formatted correctly, raises an OptionValueError"""
    if value not in organisms.split(", "):
        raise OptionValueError("It looks like organism '%s' is not one of %s, please check and try again." % (value, organisms))
    setattr(_parser.values, option.dest, value)

def check_ontology(option, opt_str, value, _parser):
    """Checks if the ontology area is formatted correctly, raises an OptionValueError"""
    if value not in ontologies.split(", "):
        raise OptionValueError("It looks like '%s' is not in the list of ontologies (%s), please check and try again." % (value, ontologies))
    setattr(_parser.values, option.dest, value)

def check_file(option, opt_str, value, _parser):
    """See if a file exists on the file system, raises an OptionValueError"""
    if not os.access(value, os.R_OK):
        raise OptionValueError("Cannot read file %s. Please check permissions and whether it exists." % value)
    setattr(_parser.values, option.dest, value)

def check_dir(option, opt_str, value, _parser):
    """See if a directory exists on the file system, raises an OptionValueError"""
    if not os.access(value, os.W_OK):
        #raise OptionValueError("Cannot write to directory %s. Please check permissions and whether it exists." % value)
        pass
    value = os.path.abspath(value)
    setattr(_parser.values, option.dest, value)

def round_pvalue(option, opt_str, value, _parser):
    """Rounds the pvalue option to the nearest available 1E-x, raises an OptionValueError"""
    try:
        pvalue = float(value)
    except ValueError:
        raise OptionValueError("Please enter a p-value threshold between 1E-2 (0.01) and 1E-11 (0.0000000001)")
    ex = int(round(math.log(value, 10)))
    if ex > -2:
        ex = -2
    elif ex < -11:
        ex = -11
    value = "%.*f" % (abs(ex), 10**ex)
    setattr(_parser.values, option.dest, value)

def optionparse(parser):
    """Creates an option dictionary from the argv"""
    parser.add_option("-o", "--organism", action = "callback", callback = check_organism,  
                      dest = "organism", type = "string",
                      help = "One of " + organisms + " [default %default]")    
    parser.add_option("-p", "--pvalue", action = "callback", callback = round_pvalue, 
                      dest = "pvalue", type = "float",
                      help = "p-value threshold p, 1E-2 > p > 1E-11 [default %default]")    
    parser.add_option("-l", "--list", dest = "list_file", action = "callback", 
                      callback = check_file,  type = "string", 
                      help = "path to the (target) gene list file [default %default]")
    parser.add_option("-d", "--directory", dest = "directory", action = "callback", 
                      callback = check_dir,  type = "string", 
                      help = "path to the html output directory [default '%default']")
    parser.add_option("-b", "--background", dest = "background_file", type = "string", 
                      action = "callback", callback = check_file, 
                      help = "path to the background list file, only supply this if the target list is unranked")
    parser.add_option("-a", "--area", dest = "area", type = "string", action = "callback",
                      callback = check_ontology,  
                      help="ontology area, one of proc(ess), func(tion), comp(onent), all [default %default]")
    parser.add_option("-g", "--galaxyfile", dest = "galaxyfile", type = "string", action = "store",
                      help = "Galaxy-specific: path to the dataset file, refers to index.html")
    parser.set_defaults(organism = organisms.split(", ")[0], pvalue = "0.001", 
                        list_file = "genes.txt", area = "all", background_file = None,
                        directory = ".", galaxyfile="")

def readgenes(filename):
    """Reads a text file's contents and returns a single string"""
    fileh = open(filename)
    contents = fileh.read()
    fileh.close()
    return contents

def dataencode(options):
    """Encode the supplied command line parameters into HTTP POST data"""
    targets = readgenes(options.list_file)
    # set mode and data for 2-file version (target, background)
    if options.background_file is not None:
        background = readgenes(options.background_file)
        data = urlencode({"species" : options.organism, "target_set" : targets, 
                          "ontology" : options.area, "pvalue_thresh" : options.pvalue,
                          "background_set" : background, "run_mode" : "hg",
                          "target_file_name" : "", "background_file_name" : "",
                          "run_gogo_button" : True, "output_excel" : False, 
                          "output_unresolved" : False, "output_revigo" : False,
                          "fast_mode" : True})    
    #set mode and data for ranked list version
    else:
        data = urlencode({"species" : options.organism, "target_set" : targets, 
                          "ontology" : options.area, "pvalue_thresh" : options.pvalue,
                          "run_mode" : "mhg", "target_file_name" : "",
                          "background_set" : "", "background_file_name" : "",
                          "run_gogo_button" : True, "output_excel" : False, 
                          "output_unresolved" : False, "output_revigo" : False,
                          "fast_mode" : True})
    return data
    
def spider(data, directory):    
    """Submits GOrilla form, spiders result into directory"""
    # submit query to the webserver
    conn = urlopen(gorillaurl, data)
    resulturl = conn.geturl()
    conn.close()
    # wait 15 secs for the server to finish and retrieve all result files
    sleep(15)
    sp = Spider()
    sp.webmirror(root = directory, base = resulturl, t = 1, depth = 1)
    sleep(5)
    # consolidate the downloaded files, delete unneccesary ones
    subdir1 = os.path.join(directory, "GOrilla")
    subdir = os.path.join(subdir1, resulturl.split("=")[1])
    shutil.move(os.path.join(directory, "GOResults.html"), os.path.join(subdir, "index.html"))
    shutil.rmtree(os.path.join(directory, "pics"))
    for delfile in ("example.html", "help.html", "index.html", "news.html", "VantVeerMoreLess5.txt"):
        os.remove(os.path.join(directory, delfile))
    resfiles = os.listdir(subdir)
    for resfile in resfiles:
        shutil.move(os.path.join(subdir, resfile), os.path.join(directory, resfile))
    shutil.rmtree(subdir1)
    # finally, remove the back link from the top frame
    topfile = open(os.path.join(directory, "top.html"))
    top = topfile.readlines()
    topfile.close()
    topfile = open(os.path.join(directory, "top.html"), "w")
    for line in top:
        if not line.startswith('<td align="center" class="link" id="Back"><h2><a href="http://cbl-gorilla.cs.technion.ac.il/"'):
            topfile.write(line)
    topfile.close()

def galaxyindex(dataset):
    """Writes a referring html page to index.html for Galaxy's result folder structure"""
    try:
        datasetf = open(dataset, "w")
        datasetf.write(galaxydatasetindexhtml)
        datasetf.close()
    except:
        print "Could not write to galaxy index dataset -- check if file exists"
     

def main(argv = None):
    """Main routine called by default with optional, user-definable argv arguments"""     
    # set up command line options
    if argv is None:         
        argv = sys.argv    
    parser = OptionParser(version = "%prog " + __version__, usage = __usage__)
    optionparse(parser)
    (options, args) = parser.parse_args()
    # put together the HTTP POST data object
    data = dataencode(options)
    # execute and download results, put into results dir
    spider(data, options.directory)
    if not options.galaxyfile == "":
        galaxyindex(options.galaxyfile)
        
if __name__ == "__main__":
    main()


