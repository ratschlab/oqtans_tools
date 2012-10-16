"""
#######################################################################################
#                                                                                     #
#    kirmes.py is a command-line front-end to the KIRMES pipeline                     #
#    BibTeX entries below. Please cite:                                               #
#    Sebastian J. Schulheiss, Wolfgang Busch, Jan U. Lohmann, Oliver Kohlbacher,      #
#    and Gunnar Raetsch (2008) KIRMES: Kernel-based identification of regulatory      #
#    modules in euchromatic sequences. In Andreas Beyer and  Michael Schroeder (Eds.) #
#    German Conference on Bioinformatics, 158-167, GI, Springer Verlag Heidelberg.    # 
#                                                                                     #
#    Copyright (C) 2007-2009 Sebastian J. Schultheiss <sebi@umich.edu>                #
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

__version__ = "0.8.0"
__license__  = "GNU General Public License"

# bibtex entry
__author__ = """
@article{Schultheiss2009KIRMES,
    author = {Sebastian J. Schultheiss and Wolfgang Busch and Jan U. Lohmann and Oliver Kohlbacher and Gunnar Raetsch},
    title = {{KIRMES}: Kernel-based identification of regulatory modules in euchromatic sequences},
    year = {2009},
    journal = {Bioinformatics},
    publisher = {Oxford Journals},
    month = {April},
    pages = {btp278},
    doi = {10.1093/bioinformatics/btp278},
    abstract = {Motivation: Understanding transcriptional regulation is one of the 
                main challenges in computational biology. An important problem is 
                the identification of transcription factor binding sites in promoter 
                regions of potential transcription factor target genes. It is 
                typically approached by position weight matrix-based motif 
                identification algorithms using Gibbs sampling, or heuristics to 
                extend seed oligos. Such algorithms succeed in identifying single, 
                relatively well-conserved binding sites, but tend to fail when it 
                comes to the identification of combinations of several degenerate 
                binding sites, as those often found in cis-regulatory modules.  
                Results: We propose a new algorithm that combines the benefits of 
                existing motif finding with the ones of Support Vector Machines (SVMs) 
                to find degenerate motifs in order to improve the modeling of 
                regulatory modules. In experiments on microarray data from Arabidopsis 
                thaliana, we were able to show that the newly developed strategy 
                significantly improves the recognition of transcription factor targets.  
                Availability: The PYTHON source code (open source-licensed under GPL), 
                the data for the experiments and a Galaxy-based web service are 
                available at http://www.fml.mpg.de/raetsch/projects/kirmes.  
                Contact: sebastian.schultheiss@tuebingen.mpg.de},
    URL = {http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btp278v1}
}
"""

__usage__ = "%prog " + __version__ + """
  For a training run, supply positives and negatives: 
        %prog -t -p positives.fasta -n negatives.fasta [options]
          
  For a testing or evaluation run, supply the saved training file and a query file: 
        %prog -e -q query.fasta -r trained_kirmes.pck [options]

  You may combine training and evaluating in a single command line, if 
  you specify the -e option and supply all necessary files."""


# system imports
import os
from optparse import OptionParser, OptionGroup, OptionValueError
import sys
#from cPickle import Pickler, Unpickler

try:
    from numpy import ones, int32, array, append, float64, argsort
except ImportError:
    print "ImportError:"
    print "Importing numpy failed. The numerical python package is required."
    raise

try:
    # own imports
    import kirmes_ini
    from MotifFinder import MotifFinderSettings, MotifFinder
except ImportError:
    print "ImportError:"
    print "One of the required imports failed, please make sure all KIRMES files "
    print "are present in the current directory, or download this program again."
    raise

try:
    import EasySVM
    # shogun imports
    from shogun.Kernel import GaussianKernel, WeightedDegreePositionStringKernel, CombinedKernel
    from shogun.Features import Labels, CombinedFeatures
    from shogun.Classifier import LibSVM
except ImportError:
    print "ImportError:"
    print "Importing the large-scale machine learning toolbox SHOGUN failed."
    print "Please make sure that the shogun installation is valid and"
    print "that shogun is installed for python using the configure option, e.g."
    print "./configure --interfaces=libshogun,libshogunui,python,python_modular"
    print "or refer to instructions on http://www.shogun-toolbox.org/"
    raise

def check_file(option, opt_str, value, _parser):
    """See if a file exists on the file system, raises an OptionValueError"""
    if not os.path.isfile(value):
        raise OptionValueError("Cannot open %s as a file. Please check if it exists." % value)
    setattr(_parser.values, option.dest, value)

def optionparse(parser):
    """Adds all options to a command-line parser object"""
    testgroup = OptionGroup(parser, "Options for prediction and evaluation",
                        "These options are required during a prediction run. ")
    testgroup.add_option("-e", "--evaluate", help = "do a prediction run with KIRMES", 
                      action = "store_true", dest = "eval", default = False)
    testgroup.add_option("-q", "--query", type = "string", action = "callback", callback = check_file, dest = "query", 
                      help="path to the fasta file with a set of regulatory regions to predict [default %default]")
    testgroup.add_option("-u", "--qgff", type = "string", action = "callback", callback = check_file, dest = "qgff", 
                      help="path to the gff file with motif positions [default %default]")
    parser.add_option_group(testgroup)
    
    traingroup = OptionGroup(parser, "Options for training",
                        "These options are required during a training run. ")
    traingroup.add_option("-p", "--positives", dest = "positives",
                          action = "callback", callback = check_file, type = "string", 
                          help="path to the fasta file with a positive set of \
                                regulatory regions [default %default]")
    traingroup.add_option("-n", "--negatives", dest = "negatives", type = "string", 
                          action = "callback", callback = check_file, 
                          help="path to the fasta file with a negative set of \
                                regulatory regions [default %default]")
    traingroup.add_option("-s", "--sigma", type = "float", dest = "sigma", 
                          help = "sigma value for the RBF kernel [default %default]")
    traingroup.add_option("-w", "--width", type = "int", dest = "window_width", 
                          help = "width of the cut-out windows of the WDS kernel [default %default]")
    traingroup.add_option("-m", "--motif-contrib", dest = "contrib", action="store_true",
                          help = "show a summary of motif contributions to classification",
                          default = False)
    traingroup.add_option("-l", "--sequence-logos", dest = "logos", action="store_true",
                          help = "show sequence logos for each motif kernel", default = False)
    traingroup.add_option("-i", "--pgff", dest = "pgff", type = "string", 
                          help="path to the gff file of motif positions from the \
                                positive regulatory regions [default %default]")
    traingroup.add_option("-j", "--ngff", dest = "ngff", type = "string", 
                          help="path to the gff file of motif positions from the \
                                negative regulatory regions [default %default]")
    parser.add_option_group(traingroup)
        
    parser.add_option("-r", "--replace", dest = "replace", type = "string", 
                      help = "Replace any non-nucleotide characters with either A, C, G, T, or randomly [default %default]")
    parser.add_option("-c", "--conservation", dest = "conservation_file", type = "string", 
                      help = "path to a file with conservation information in bed format [default %default]")
    parser.add_option("-k", "--kernel", dest = "kernel_file", type = "string", 
                      help = "path to the kernel parameters, output of a KIRMES training run [default %default]")
    parser.add_option("-o", "--output", dest = "output_html", type = "string", 
                      help = "path to the output html file with a summary [default %default]")
    parser.add_option("-y", "--img-dir", dest = "output_path", type = "string",
                      help = "image directory for motif images")
    parser.set_defaults(positives = kirmes_ini.POSITIVES_FILENAME, 
                        negatives = kirmes_ini.NEGATIVES_FILENAME,
                        query = kirmes_ini.QUERY_FILENAME, 
                        kernel_file = kirmes_ini.KERNEL_FILENAME,
                        sigma = kirmes_ini.RBF_SIGMA, 
                        nof_motifs = kirmes_ini.NOF_MOTIFS,
                        motif_length = kirmes_ini.MOTIF_LENGTH,
                        window_width = kirmes_ini.MOTIF_WINDOW_WIDTH,
                        ngff = kirmes_ini.NGFF_FILENAME,
                        pgff = kirmes_ini.PGFF_FILENAME,
                        qgff = kirmes_ini.QGFF_FILENAME,
                        conservation_file = kirmes_ini.CONSERVATION_FILENAME,
                        output_file = kirmes_ini.OUTPUT_FILENAME,
                        replace = kirmes_ini.REPLACEMENT_CHARACTER)
    
def htmlize(res_dict, outfilename):
    """Makes the result into an HTML string"""
    outstr = "<html>\n<head>\n  <title>KIRMES Result File -- Classification</title>\n</head>\n<body>\n"
    outstr += "<h1>KIRMES Result File</h1>\n<table>\n"
    if "query" in res_dict:
        res = res_dict["query"]
        outstr += "<tr>\n  <th>Gene FASTA ID</th>\n  <th>Classification</th>\n  <th>+/-1</th>\n</tr>\n"
        resarray = res.split("\n")
        for line in resarray:
            tds = line.split("\t")
            if len(tds) < 3: 
                continue
            outstr += "<tr>\n  <td>" + tds[0]
            outstr += "</td><td>" + tds[1] + "</td><td align='right'>"
            outstr += tds[2] + "</td>\n</tr>"
        outstr += "\n</table>"
    if "contrib" in res_dict:
        outstr += res_dict["contrib"]
    if "poims" in res_dict:
        outstr += res_dict["poims"]
    outstr += "\n</body>\n</html>\n"
    outfile = open(outfilename, "w")
    outfile.write(outstr)
    outfile.close()

def training_run(options):
    """Conduct a training run and return a trained SVM kernel"""
    settings = MotifFinderSettings(kirmes_ini.MOTIF_LENGTH, 
                                   options.window_width, 
                                   options.replace)
    positives = MotifFinder(finder_settings = settings)
    positives.setFastaFile(options.positives)
    positives.setMotifs(options.pgff)
    pmotifs, ppositions = positives.getResults()
    negatives = MotifFinder(finder_settings = settings)
    negatives.setFastaFile(options.negatives)
    negatives.setMotifs(options.ngff)
    nmotifs, npositions = negatives.getResults()
    
    wds_kparams = kirmes_ini.WDS_KERNEL_PARAMETERS
    wds_svm = EasySVM.EasySVM(wds_kparams)
    num_positives = len(pmotifs.values()[0])
    num_negatives = len(nmotifs.values()[0])
    #Creating Kernel Objects    
    kernel = CombinedKernel()
    features = CombinedFeatures()
    kernel_array = []
    motifs = pmotifs.keys()
    motifs.sort()
    #Adding Kmer Kernels
    for motif in motifs:
        all_examples = pmotifs[motif] + nmotifs[motif]
        motif_features = wds_svm.createFeatures(all_examples)
        wds_kernel = WeightedDegreePositionStringKernel(motif_features, motif_features, \
                                                        wds_kparams['degree'])
        wds_kernel.set_shifts(wds_kparams['shift'] * ones(wds_kparams['seqlength'], dtype = int32))
        features.append_feature_obj(motif_features)
        kernel_array.append(wds_kernel)
        kernel.append_kernel(wds_kernel)
    rbf_svm = EasySVM.EasySVM(kirmes_ini.RBF_KERNEL_PARAMETERS)
    positions = array(ppositions + npositions, dtype = float64).T
    position_features = rbf_svm.createFeatures(positions)
    features.append_feature_obj(position_features)
    motif_labels = append(ones(num_positives), -ones(num_negatives))
    complete_labels = Labels(motif_labels)
    rbf_kernel = GaussianKernel(position_features, position_features, \
                                kirmes_ini.RBF_KERNEL_PARAMETERS['width'])
    kernel_array.append(rbf_kernel)
    kernel.append_kernel(rbf_kernel)
    #Kernel init
    kernel.init(features, features)
    kernel.set_cache_size(kirmes_ini.K_CACHE_SIZE)
    svm = LibSVM(kirmes_ini.K_COMBINED_C, kernel, complete_labels)
    svm.parallel.set_num_threads(kirmes_ini.K_NUM_THREADS)
    #Training
    svm.train()
    if not os.path.exists(options.output_path):
        os.mkdir(options.output_path)
    html = {}
    if options.contrib:
        html["contrib"] = contrib(svm, kernel, motif_labels, kernel_array, motifs)
    if options.logos: 
        html["poims"] = poims(svm, kernel, kernel_array, motifs, options.output_path) 
    if options.query:
        html["query"] = evaluate(options, svm, kernel, features, motifs)
    htmlize(html, options.output_html)

def contrib(svm, kernel, motif_labels, kernel_array, motifs):  
    """Calculate the contribution of each kernel"""
    (contrib_str, baseline_roc, baseline_prc, mean_acc) = EasySVM.evaluate(svm.classify().get_labels(),
                                                                           None, motif_labels)
    html = "\n\n<h2>Training Error Evaluation</h2>\n<pre>" + contrib_str + "</pre>\n"
    rocs = []
    prcs = []
    for i in xrange(kernel.get_num_subkernels()):
        kernel_array[i].set_combined_kernel_weight(0.0)
        training_less_1 = svm.classify().get_labels()
        kernel_array[i].set_combined_kernel_weight(1.0)
        (res_str, mean_roc, mean_prc, mean_acc) = EasySVM.evaluate(training_less_1, None, motif_labels)
        rocs.append(baseline_roc - mean_roc)
        prcs.append(baseline_prc - mean_prc)
    #output formatting
    html += "<table>\n<tr><th>Motif Kernel</th>\n<th>PRC</th>\n</tr>\n"
    motifs.append("positional information (RBF)")
    order = argsort(prcs)
    for i in xrange(len(prcs) - 1, -1, -1):
        html += "<tr><td>Kernel for " + motifs[order[i]] + "</td>\n"
        #html += "<td>" + str(rocs[order[i]]) + "</td>\n"
        html += "<td>" + str(prcs[order[i]]) + "</td></tr>\n"
    html += "</table>\n"
    motifs.pop()
    return html
    
def poims(svm, kernel, kernel_array, motifs, path):
    """Plot POIMs and sequence logos from WD kernels"""
    html = "\n\n<h2>POIMs and Sequence Logos</h2>\n"
    html += "<table>\n<tr><th>Motif</th>\n<th>POIM</th>\n<th>Sequence Logo</th>\n</tr>\n"
    for i in xrange(kernel.get_num_subkernels() - 1):
        (poim, max_poim, diff_poim, poim_totalmass) = \
            EasySVM.computePOIMs(svm, kernel_array[i], kirmes_ini.MOTIF_LENGTH, 
                                 kirmes_ini.MOTIF_WINDOW_WIDTH)
        motif_filename = os.path.join(path, motifs[i])
        EasySVM.plotPOIMs(motif_filename + "_poim.png", 
                          poim, max_poim, diff_poim, poim_totalmass, 
                          kirmes_ini.MOTIF_LENGTH + 1, kirmes_ini.MOTIF_WINDOW_WIDTH)
        EasySVM.weblogoPOIM(motif_filename + ".png", poim, kirmes_ini.MOTIF_WINDOW_WIDTH)
        html += "<tr><td>" + motifs[i] + "</td>\n<td>"
        rel_src = os.path.split(motif_filename)[1]
        html += "<img src='" + rel_src + "_poim.png' /></td>\n<td>"
        html += "<img src='" + rel_src + ".png' width='350' /></td></tr>\n"
    html += "</table>\n"
    return html

def evaluate(options, svm, kernel, features, motifs):
    """Evaluate examples using a trained kernel"""
    query = MotifFinder(finder_settings = MotifFinderSettings(kirmes_ini.MOTIF_LENGTH, options.window_width))
    query.setFastaFile(options.query)
    query.setMotifs(options.qgff)
    qmotifs, qpositions = query.getResults()
    feats_query = CombinedFeatures()
    wds_svm = EasySVM.EasySVM(kirmes_ini.WDS_KERNEL_PARAMETERS)
    try:
        assert set(qmotifs.keys()).issuperset(set(motifs))
    except AssertionError:
        print "The motif positions in the query sequence are incomplete, there are no positions for:"
        print set(motifs).difference(qmotifs.keys())
        raise
    for motif in motifs:
        feats_query.append_feature_obj(wds_svm.createFeatures(qmotifs[motif]))
    query_positions = array(qpositions, dtype = float64)
    query_positions = query_positions.T
    rbf_svm = EasySVM.EasySVM(kirmes_ini.RBF_KERNEL_PARAMETERS)
    feats_query.append_feature_obj(rbf_svm.createFeatures(query_positions))
    kernel.init(features, feats_query)
    out = svm.classify().get_labels()
    qgenes = query.getGenes()
    ret_str = ""
    print "#example\toutput\tsplit"
    for i in xrange(len(out)):
        if out[i] >= 0:
            classif = "\tpositive\t"
        else:
            classif = "\tnegative\t"
        ret_str += qgenes[i] + classif + str(out[i]) + "\n"
        print str(i) + "\t" + str(out[i]) + "\t0"
    return ret_str
    
def main(argv = None):     
    """    
    #######################
    # main function       #
    #######################"""
    if argv is None:         
        argv = sys.argv    
    parser = OptionParser(version = "%prog " + __version__, usage = __usage__)
    optionparse(parser)
    (options, args) = parser.parse_args()
    training_run(options)

if __name__ == "__main__":
    main()





   
    
    

