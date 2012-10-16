"""
#############################################################################################
#                                                                                           #
#    This class is part of the MLB-Galaxy package, adding some sequence analysis            #
#    functionality to PSU's Galaxy framework.                                               #
#    Copyright (C) 2008 Cheng Soon Ong <chengsoon.ong@tuebingen.mpg.de>                     #
#    Copyright (C) 2008 Gunnar Raetsch <Gunnar.Raetsch@tuebingen.mpg.de>                    #
#    Copyright (C) 2007, 2009 Sebastian J. Schultheiss <sebi@umich.edu                      #
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
#  Original Author: Sebastian J. Schultheiss, version 0.8.0                                 #
#  Please add a notice of any modifications here:                                           #
#     Gunnar Raetsch: rewrote code for training on sequences to be read from files          #
#     Cheng Soon Ong: Added code for educational toolbox                                    #
#     Sebastian J. Schultheiss: Class-ified for KIRMES use in 02/2009                       #
#     Sebastian J. Schultheiss: Updated for shogun-0.9.1 in  02/2010                        #
#     Sebastian J. Schultheiss: Tweaks for Galaxy Integration in 12/2010                    #
#                                                                                           #
#############################################################################################
"""
__version__ = "EasySVM.py 0.8.0"

import os
import random
import shutil
import warnings
import numpy
from numpy import ones

from shogun.Kernel import GaussianKernel, WeightedDegreePositionStringKernel
from shogun.Kernel import LinearKernel, PolyKernel, LocalAlignmentStringKernel, LocalityImprovedStringKernel, CommWordStringKernel
from shogun.Features import RealFeatures, Labels, StringCharFeatures, DNA, StringWordFeatures
from shogun.Classifier import LibSVM
from shogun.PreProc import SortWordString
from shogun.Evaluation import PerformanceMeasures


class EasySVM(object):
    """A wrapper around shogun SVM objects like the kernel, features ..."""
    
    def __init__(self, kparam = {}, kernel_file = None):
        "Initialize with default parameters (None)"
        self.kparam = kparam
        self.kernel = None
        if kernel_file is not None:
            self.setKernel(kernel_file)
        
    def getKernel(self, kernel_file = None):
        "Writes the kernel to a file"
        return pickle.dump(self.kernel, kernel_file)
    
    def setKernel(self, kernel_file):
        "Load a kernel from a file"
        self.kernel = pickle.load(kernel_file)
    
    def getKernelName(self):
        "Return the kernel parameter \'Kernel name\'"
        return self.kparam['name']
    
    def setKernelName(self, name):
        "Set the kernel name/type"
        self.kparam['name'] = name
        
    def getKernelParameter(self, parameter_name = None):
        """Returns any of the following kernel parameters:
            kernelname, 
            width, 
            modelsel_name, modelsel_params,
            degree, scale,
            inhomogene, normal
            shift, seqlength, 
            indeg, outdeg,
            C,
            poimdegree, ...
        """
        if parameter_name is None:
            return self.kparam
        else:
            return self.kparam[parameter_name]
        
    def setKernelParameter(self, parameter_name = None, parameter_value = None):
        """Set an arbitrary kernel parameter, 
        or set the whole kparam dictionary if parameter_name is empty"""
        if parameter_name is None and parameter_value is not None:
            self.kparam = parameter_value
        elif parameter_name is not None:
            self.kparam[parameter_name] = parameter_value
            
    def parseKernelParameters(self, parameter_string, model_selection = False):
        """Parse the arguments for a particular kernel"""
        parameters = parameter_string.split(" ")
        kernelname = parameters[0] 
        self.kparam = {}
        self.kparam["name"] = kernelname
        self.kparam["modelsel_name"] = None
        self.kparam["modelsel_params"] = None
        
        if kernelname == 'gauss':
            if len(parameters) < 2:
                raise ValueError('Not enough arguments for a Gauss-type kernel.\nUsage: gauss <width>\n')
            if model_selection:
                self.kparam['width'] = None
                self.kparam["modelsel_name"] = "width"
                self.kparam["modelsel_params"] = parseFloatList(parameters[1])
            else:
                self.kparam['width'] = float(parameters[1])
        elif kernelname == 'linear':
            self.kparam['scale'] = 1
            # no parameters
        elif kernelname == 'poly':
            if len(parameters) < 4:
                raise ValueError('Not enough arguments for a polynomial kernel.\nUsage: poly <degree> <true|false> <true|false>\n')
            if model_selection:
                self.kparam['degree'] = None
                self.kparam["modelsel_name"] = "degree"
                self.kparam["modelsel_params"] = parseIntList(parameters[1])
            else:
                self.kparam['degree'] = int(parameters[1])
            self.kparam['inhomogene'] = (parameters[2] == 'true')
            self.kparam['normal'] = (parameters[3] == 'true')
        elif kernelname == 'wd':
            if len(parameters) < 3:
                raise ValueError('Not enough arguments for a WD kernel.\nUsage: wd <degree> <shift>\n')
            if model_selection:
                self.kparam['degree'] = None
                self.kparam["modelsel_name"] = "degree"
                self.kparam["modelsel_params"] = parseIntList(parameters[1])
            else:
                self.kparam['degree'] = int(parameters[1])
            if model_selection and len(self.kparam["modelsel_params"]) == 1:
                self.kparam['degree'] = self.kparam["modelsel_params"][0]
                self.kparam['shift'] = None
                self.kparam["modelsel_name"] = "shift"
                self.kparam["modelsel_params"] = parseIntList(parameters[2])
            else:
                self.kparam['shift'] = int(parameters[2])
        elif kernelname == 'spec':
            if len(parameters) < 2:
                raise ValueError('Not enough arguments for a Spectrum kernel.\nUsage: spec <degree>\n')
            if model_selection:
                self.kparam['degree'] = None
                self.kparam["modelsel_name"] = "degree"
                self.kparam["modelsel_params"] = parseIntList(parameters[1])
            else:
                self.kparam['degree'] = int(parameters[1])
        elif kernelname == 'localalign':
            # no parameters
            pass
        elif kernelname == 'localimprove':
            if len(parameters) < 4:
                raise ValueError('Not enough arguments for a localimprove kernel.\nUsage: localimprove <length> <indegree> <outdegree>\n')
            self.kparam['length'] = int(parameters[1])
            if model_selection:
                self.kparam['width'] = None
                self.kparam["modelsel_name"] = "indeg"
                self.kparam["modelsel_params"] = parseIntList(parameters[2])
            else:
                self.kparam['indeg'] = int(parameters[2])
            self.kparam['outdeg'] = int(parameters[3])
        else:
            raise ValueError('Unknown kernel name \"' + kernelname + '\" in the parameter_string\n')
            
    def setC(self, C):
        "Set the oft-used kernel parameter C"
        self.setKernelParameter('C', C)
    
    def getC(self):
        "Return the current value for kernel parameter C"
        return self.getKernelParameter('C')
        
    def createFeatures(self, examples):
        """Converts numpy arrays or sequences into shogun features"""
        if self.kparam['name'] == 'gauss' or self.kparam['name'] == 'linear' or self.kparam['name'] == 'poly':
            examples = numpy.array(examples)
            feats = RealFeatures(examples)
            
        elif self.kparam['name'] == 'wd' or self.kparam['name'] == 'localalign' or self.kparam['name'] == 'localimprove':
            #examples = non_atcg_convert(examples, nuc_con)
            feats = StringCharFeatures(examples, DNA)
        elif self.kparam['name'] == 'spec':
            #examples = non_atcg_convert(examples, nuc_con)
            feats = StringCharFeatures(examples, DNA) 
       
            wf = StringUlongFeatures( feats.get_alphabet() )
            wf.obtain_from_char(feats, kparam['degree']-1, kparam['degree'], 0, kname=='cumspec')
            del feats
    
            if train_mode:
                preproc = SortUlongString()
                preproc.init(wf)
            wf.add_preproc(preproc)
            ret = wf.apply_preproc()
            feats = wf 
    
        else:
            print 'Unknown kernel %s' % self.kparam['name']
            raise ValueError
        
        return feats

    def createKernel(self, feats_train):
        """Call the corresponding constructor for the kernel"""
    
        if self.kparam['name'] == 'gauss':
            kernel = GaussianKernel(feats_train, feats_train, self.kparam['width'])
        elif self.kparam['name'] == 'linear':
            kernel = LinearKernel(feats_train, feats_train, self.kparam['scale'])
        elif self.kparam['name'] == 'poly':
            kernel = PolyKernel(feats_train, feats_train, self.kparam['degree'], 
                                self.kparam['inhomogene'], self.kparam['normal'])
        elif self.kparam['name'] == 'wd':
            kernel = WeightedDegreePositionStringKernel(feats_train, feats_train, self.kparam['degree'])
            kernel.set_shifts(self.kparam['shift'] * numpy.ones(self.kparam['seqlength'], dtype=numpy.int32))
        elif self.kparam['name'] == 'spec':
            kernel = CommWordStringKernel(feats_train, feats_train)
        elif self.kparam['name'] == 'localalign':
            kernel = LocalAlignmentStringKernel(feats_train, feats_train)
        elif self.kparam['name'] == 'localimprove':
            kernel = LocalityImprovedStringKernel(feats_train, feats_train, self.kparam['length'], \
                                                  self.kparam['indeg'], self.kparam['outdeg'])
        else:
            print 'Unknown kernel %s' % self.kparam['name']
            raise ValueError
        self.kernel = kernel
        return kernel

    def __str__(self):
        """Generates a short string describing the model parameters"""
        if self.kparam["modelsel_name"] is None or len(self.kparam["modelsel_params"]) == 1:
            string = "\tC=%1.1f" % self.kparam['C']
        else:
            string = "\tC=%1.1f\t%s=%i" % (self.kparam['C'], self.kparam["modelsel_name"])
        return string
    
    def model2str(self, C, kp):
        """Generates a string describing the model parameters"""
        if self.kparam["modelsel_name"] == None or len(self.kparam["modelsel_params"]) == 1:
            string = "\tC=%1.1f" % C
        else:
            if type(kp) == type(int(0)):
                string = "\tC=%1.1f\t%s=%i" % (C, self.kparam["modelsel_name"], kp)
            else:
                string = "\tC=%1.1f\t%s=%1.2f" % (C, self.kparam["modelsel_name"], kp)
        return string

    def train(self, trainexamples, trainlabels):
        """Trains a SVM with the given kernel"""
        kernel_cache_size = 500
        num_threads = 6

        feats_train = self.createFeatures(trainexamples)
        if self.kparam['name'] == 'wd':
            self.kparam['seqlength'] = len(trainexamples[0])
        self.createKernel(feats_train)
        
        self.kernel.io.disable_progress()
        self.kernel.set_cache_size(int(kernel_cache_size))
    
        labels = Labels(numpy.array(trainlabels, numpy.double))
    
        svm = LibSVM(self.getC(), self.kernel, labels)
        svm.parallel.set_num_threads(num_threads)
        svm.io.disable_progress()
        svm.train()
    
        return (svm, feats_train)
    
    def trainAndTest(self, trainexamples, trainlabels, testexamples):
        """Trains a SVM with the given kernel, and predict on the test examples"""
    
        (svm, feats_train) = self.train(trainexamples, trainlabels) #,C,kname,kparam)
        feats_test = self.createFeatures(testexamples)
        self.kernel.init(feats_train, feats_test)
        output = svm.classify().get_labels()
        return output

    def crossvalidation(self, all_examples, all_labels, xval = 5):
        """Perform cross validation using an SVM
        xval -- the number of folds    
        """
        print 'Using %i-fold crossvalidation' % xval
        
        partitions = getPartitionedSet(len(all_labels), xval)
        all_outputs = [0.0] * len(all_labels)
        all_split = [-1] * len(all_labels)
    
        for repetition in xrange(xval):
            XT, LT, XTE, LTE = getCurrentSplit(repetition, partitions, all_labels, all_examples)
            del LTE
            svmout = self.trainAndTest(XT, LT, XTE)
            
            for i in xrange(len(svmout)):
                all_outputs[partitions[repetition][i]] = svmout[i]
                all_split[partitions[repetition][i]] = repetition
        return (all_outputs, all_split)
#


################################################################################
# main functions

    def crossvalidationSVM(self, xval, examples, labels):
        """A top level member function to run cross validation"""
        # run cross-validation
        (all_outputs, all_split) = self.crossvalidation(xval, examples, labels)
        res_str = '#example\toutput\tsplit\n'
        for ix in xrange(len(all_outputs)):
            res_str += '%d\t%2.7f\t%d\n' % (ix, all_outputs[ix], all_split[ix])
        return res_str

    def modelSelectionSVM(self, xval, examples, labels, Crange):
        """A top level member function to run model selection"""
        # run cross-validation
        mean_rocs = [] 
        mean_prcs = []
        mean_accs = [] 
        all_Cs = [] 
        all_kparam = [] 
    
        if self.kparam["modelsel_name"] == None:
            for C in Crange:
                self.setC(C)
                (all_outputs, all_split) = self.crossvalidation(xval, examples, labels)
                (res_str, mean_roc, mean_prc, mean_acc) = self.evaluate(all_outputs, all_split, labels)
                del res_str
                mean_rocs.append(mean_roc) 
                mean_prcs.append(mean_prc) 
                mean_accs.append(mean_acc) 
                all_Cs.append(C) 
                all_kparam.append(None) 
        else: # also optimize one kernel parameter
            for C in Crange:
                for kp in self.kparam["modelsel_params"]:
                    self.kparam[self.kparam["modelsel_name"]] = kp 
                    self.setC(C)
                    (all_outputs, all_split) = self.crossvalidation(xval, examples, labels)
                    (res_str, mean_roc, mean_prc, mean_acc) = self.evaluate(all_outputs, all_split, labels)
                    del res_str
                    mean_rocs.append(mean_roc) 
                    mean_prcs.append(mean_prc) 
                    mean_accs.append(mean_acc) 
                    all_Cs.append(C) 
                    all_kparam.append(kp)
    
        max_roc = numpy.max(numpy.array(mean_rocs)) 
        max_prc = numpy.max(numpy.array(mean_prcs)) 
        max_acc = numpy.max(numpy.array(mean_accs)) 
    
        if self.kparam["modelsel_name"] == None or len(self.kparam["modelsel_params"]) == 1:
            detail_str = "\tC\tROC\tPRC\tAccuracy (at threshold 0)\n"
        else:
            detail_str = "\tC\t%s\tROC\tPRC\tAccuracy (at threshold 0)\n" % self.kparam["modelsel_name"]
    
        best_roc_str = ''
        best_prc_str = ''
        best_acc_str = ''
        for i in xrange(len(all_Cs)):
            # determine the best parameter combinations
            if mean_rocs[i] == max_roc:
                rocsym = '+'
                best_roc_str += self.model2str(all_Cs[i], all_kparam[i])+'\n'
            else:
                rocsym = ' '
            if mean_prcs[i] == max_prc:
                prcsym = '+'
                best_prc_str += self.model2str(all_Cs[i], all_kparam[i])+'\n'
            else:
                prcsym = ' '
            if mean_accs[i] == max_acc:
                accsym = '+'
                best_acc_str += self.model2str(all_Cs[i], all_kparam[i])+'\n'
            else:
                accsym = ' '
            detail_str += self.model2str(all_Cs[i], all_kparam[i], False)+'\t'
            if self.kparam["modelsel_name"] == None or len(self.kparam["modelsel_params"]) == 1:
                detail_str += '%c%2.1f%%\t%c%2.1f%%\t%c%2.1f%%\n' % (rocsym, 100*mean_rocs[i], prcsym, 100*mean_prcs[i], accsym, 100*mean_accs[i])
            else:
                detail_str += '%c%2.1f%%\t%c%2.1f%%\t%c%2.1f%%\n' % (rocsym, 100*mean_rocs[i], prcsym, 100*mean_prcs[i], accsym, 100*mean_accs[i])
    
        detail_str = ('Best model(s) according to ROC measure:\n%s' % best_roc_str) + detail_str 
        detail_str = ('\nBest model(s) according to PRC measure:\n%s' % best_prc_str) + detail_str
        detail_str = ('\nBest model(s) according to accuracy measure:\n%s' % best_acc_str) + detail_str
    
        detail_str = ('\nDetailed results:\n') + detail_str
        return detail_str

    def predictSVM(self, trainexamples, trainlabels, testexamples):
        """A top level script to parse input parameters and train and predict"""
        # run training and testing
        svmout = self.trainAndTest(trainexamples, trainlabels, testexamples)    
        # write output file
        res_str = '#example\toutput\n'
        for ix in xrange(len(svmout)):
            res_str += str(ix) + '\t' + str(svmout[ix]) + '\n'

        return res_str

    def evaluateSVM(self, trainexamples, trainlabels, prediction_file, roc_or_prc = None):
        """A top level script to parse input parameters and evaluate"""
        (predictions, splitassignments) = parsePrediction(prediction_file)
        roc_fname = None
        prc_fname = None
        if roc_or_prc is not None:
            if roc_or_prc.startswith('roc'):
                roc_fname = roc_or_prc
            elif roc_or_prc.startswith('prc'):
                prc_fname = roc_or_prc
        # run training and testing
        (res_str, mean_roc, mean_prc, mean_acc) = evaluate(predictions, splitassignments, trainlabels, roc_fname, prc_fname)
        del mean_acc
        del mean_prc
        del mean_roc
        return res_str
    
    def poimSVM(self, examples, labels, poimfile):
        """A top level script to parse input parameters and plot poims"""    
        # train svm and compute POIMs
        (svm, feats_train) = self.train(examples, labels)
        del feats_train
        (poim, max_poim, diff_poim, poim_totalmass) = computePOIMs(svm, self.kernel, self.kparam['poimdegree'], len(examples[0]))
    
        # plot poims
        plotPOIMs(poimfile, poim, max_poim, diff_poim, poim_totalmass, self.kparam['poimdegree'], len(examples[0]))
    
# independent functions .............
        
def computePOIMs(svm, kernel, poimdegree, max_len):
    """For a trained SVM, compute Position Oligomer Importance Matrices"""

    distr = ones((max_len, 4))/4
    kernel.prepare_POIM2(distr)
    
    kernel.compute_POIM2(poimdegree, svm)
    poim = kernel.get_POIM2()
    kernel.cleanup_POIM2() 

    (poim, max_poim, diff_poim) = reshapeNormalizeContribs(poim, poimdegree, max_len)
    (poim_weightmass, poim_totalmass) = computeWeightMass(poim, poimdegree, max_len)
    del poim_weightmass

    poim_totalmass = poim_totalmass/numpy.sum(poim_totalmass)

    return (poim, max_poim, diff_poim, poim_totalmass)

def computeWeightMass(C, maxOrder, seqLen):
    """POIM Function"""
    mass = numpy.zeros((maxOrder, seqLen), numpy.double)
    for i in xrange(0, maxOrder):
        mass[i, :] = sum(numpy.abs(C[i]))
    total = numpy.sum(mass)
    return (mass, total)

def reshapeNormalizeContribs(C, maxOrder, seqLen): #, background): #opts = {}):
    """POIM Function"""
    alphabetSize = 4
    Contribs = []
    l = 0
    for i in xrange(0, maxOrder):
        L = l + (alphabetSize**(i + 1)) * seqLen
        vec = C[l:L].copy() 
        Contribs.append(vec.reshape(seqLen, alphabetSize**(i + 1) ).T)
        l = L
    assert(l == len(C))
    maxContribs = numpy.zeros((maxOrder, seqLen), numpy.double)
    maxp_str = numpy.zeros((maxOrder, seqLen), numpy.int)
    for i in xrange(0, maxOrder):
        con = numpy.abs(Contribs[i]) 
        maxContribs[i, :] = numpy.max(con, axis = 0)
        maxp_str[i, :] = numpy.argmax(con, axis = 0)
    diffmaxContribs = numpy.zeros((maxOrder, seqLen), numpy.double)
    for k in xrange(1, maxOrder ):
        numsy = 4**(k + 1)
        for l in  xrange(0, seqLen-k):
            km = maxp_str[k, l] 
            A = numpy.abs(Contribs[k - 1][numpy.floor(km/4), l])
            B = numpy.abs(Contribs[k - 1][numpy.mod(km, numsy/4), l + 1])
            correction = numpy.max([A, B])
            diffmaxContribs[k, l] = maxContribs[k, l] - correction
    return (Contribs, maxContribs, diffmaxContribs)             
    
def plotROC(output, LTE, draw_random = False, figure_fname = "", roc_label = 'ROC'):
    """Uses matplotlib to plot the area under 
    the ROC curve into a supplied figure_fname file"""
    from matplotlib import use, font_manager
    use("Agg")  # matplotlib save without display
    from pylab import figure, plot, xticks, yticks, xlabel, ylabel, legend, savefig, axis
    
    figure(1, dpi = 150, figsize = (4, 4))

    pm = PerformanceMeasures(Labels(numpy.array(LTE)), Labels(numpy.array(output)))

    points = pm.get_ROC()
    points = numpy.array(points).T # for pylab.plot
    plot(points[0], points[1], 'b-', label = roc_label)
    if draw_random:
        plot([0, 1], [0, 1], 'r-', label = 'random guessing')
    axis([0, 1, 0, 1])
    ticks = numpy.arange(0., 1., .1, dtype = numpy.float64)
    xticks(ticks, size = 10)
    yticks(ticks, size = 10)
    xlabel('1 - specificity (false positive rate)', size = 10)
    ylabel('sensitivity (true positive rate)', size = 10)
    legend(loc = 'lower right', prop = font_manager.FontProperties('tiny')) 

    if figure_fname != None:
        warnings.filterwarnings('ignore', 'Could not match*')
        tempfname = figure_fname + '.png'
    savefig(tempfname)
    shutil.move(tempfname, figure_fname)
    
    auROC = pm.get_auROC()
    return auROC

def plotPRC(output, LTE, figure_fname = "", prc_label = 'PRC'):
    """Plots a precision recall curve into the supplied
    figure_fname file"""
    from matplotlib import use
    use("Agg")  # matplotlib save without display
    from pylab import figure, plot, axis, xticks, yticks, ylabel, xlabel, legend, savefig

    figure(2, dpi = 150, figsize = (4, 4))

    pm = PerformanceMeasures(Labels(numpy.array(LTE)), Labels(numpy.array(output)))

    points = pm.get_PRC()
    points = numpy.array(points).T # for pylab.plot
    plot(points[0], points[1], 'b-', label = prc_label)
    axis([0, 1, 0, 1])
    ticks = numpy.arange(0., 1., .1, dtype = numpy.float64)
    xticks(ticks, size = 10)
    yticks(ticks, size = 10)
    xlabel('sensitivity (true positive rate)', size = 10)
    ylabel('precision (1 - false discovery rate)', size = 10)
    legend(loc = 'lower right')
    
    if figure_fname != None:
        warnings.filterwarnings('ignore', 'Could not match*')
        tempfname = figure_fname + '.png'
    savefig(tempfname)
    shutil.move(tempfname, figure_fname)
    
    auPRC = pm.get_auPRC()
    return auPRC 

def weblogoPOIM(logofile, poim, max_len):
    """instead of plotting the POIM heatmap, create a weblogo from the 1st-degree poim"""
    warnings.filterwarnings('ignore', ' This call to matplotlib.use()*')
    from  corebio.data import rna_letters, dna_letters, amino_acid_letters
    from weblogolib import LogoData, LogoOptions, LogoFormat, classic, png_print_formatter
    #print "WEBLOGO!"
    #print "Writing ", logofile
    #print poim[0]
    positive_logo = []
    negative_logo = []
    for i in xrange(len(poim[0])):
        positive_logo.append([])
        negative_logo.append([])
        for j in xrange(len(poim[0][i])):
            if poim[0][i][j] < 0:
                positive_logo[i].append(0)
                negative_logo[i].append(poim[0][i][j] * -10000)
            else:
                negative_logo[i].append(0)
                positive_logo[i].append(poim[0][i][j] * 1000)
    #print "Positive logo: ", positive_logo
    #print "Negative logo: ", negative_logo
    pos_data = LogoData.from_counts('ACGT', numpy.array(positive_logo).T, None)
    neg_data = LogoData.from_counts("ACGT", numpy.array(negative_logo).T, None)
    neg_opt = LogoOptions()
    neg_opt.fineprint += " from KIRMES POIM data"
    #logoopt.number_interval = 5
    neg_opt.small_fontsize = 4 
    neg_opt.title_fontsize = 8
    neg_opt.scale_width = False
    title = os.path.split(logofile)[1]
    title = title[:title.rfind(".")]
    if "_" in title:
        title = title[title.rfind("_") + 1:]
    neg_opt.logo_title = title + " Negative Logo"
    neg_format = LogoFormat(neg_data, neg_opt)
    pos_opt = LogoOptions()
    #pos_opt.show_ends = True
    pos_opt.scale_width = False
    pos_opt.logo_title = title + " Positive Sequence Logo"
    pos_opt.show_fineprint = False
    pos_opt.color_scheme = classic
    pos_format = LogoFormat(pos_data, pos_opt)
    neg_logo = open(logofile + "n.png", 'w')
    png_print_formatter(neg_data, neg_format, neg_logo)
    neg_logo.close()
    pos_logo = open(logofile + "p.png", 'w')
    png_print_formatter(pos_data, pos_format, pos_logo)
    pos_logo.close()
    concatPNG(logofile, (logofile + "p.png", logofile + "n.png"))
    os.remove(logofile + "n.png")
    os.remove(logofile + "p.png")

def concatPNG(outfilename, infilenames):
    """Vertically concatenates a list of PNG files and writes them to outfilename.
    This function uses the width of the first image supplied as"""    
    from PIL import Image
    total_height = 0
    infiles = []
    total_imgs = 0
    for infilename in infilenames:
        try:
            infiles.append(Image.open(infilename))
        except:
            print "Error loading image " + infilename 
            raise
        total_height += infiles[total_imgs].size[1]
        total_imgs += 1
    im = Image.new("RGB", (infiles[0].size[0], total_height))
    insert_at = 0
    for i in range(total_imgs):
        im.paste(infiles[i], (0, insert_at))
        insert_at += infiles[i].size[1]
    im.save(outfilename)
    
def plotPOIMs(poimfile, poim, max_poim, diff_poim, poim_totalmass, poimdegree, max_len):
    """Plot a summary of the information in poims"""
    warnings.filterwarnings('ignore', 'Module pytz was already imported*')
    warnings.filterwarnings('ignore', ' This call to matplotlib.use()*')
    from matplotlib import use
    use("Agg")  # matplotlib save without display
    from pylab import figure, savefig, subplot, title, pcolor, colorbar, yticks, ylabel
    from pylab import axis, plot, xlabel, xticks, subplots_adjust, clf
    #import matplotlib

    figure(3, dpi = 150, figsize = (8, 3.5))

    # summary figures
    fontdict = dict(family = "cursive", weight = "bold", size = 12, y = 1.05) 
    #subplot(3, 2, 1)
    #title('Total POIM Mass', fontdict)
    #plot(poim_totalmass) 
    #ylabel('weight mass', size = 5)
    #colorbar()

    subplot(1, 2, 1)
    title('POIMs', fontdict)
    pcolor(max_poim, shading = 'flat') 
    subplots_adjust(wspace = 0.3) 
    colorbar()
    
    #subplot(3, 2, 5)
    #title('Differential POIMs', fontdict)
    #pcolor(diff_poim, shading = 'flat')

    #for plot in [3, 5]:
    #    subplot(3, 2, 3)
    ticks = numpy.arange(1., poimdegree + 1, 1, dtype = numpy.float64)
    ticks_str = [] 
    for i in xrange(0, poimdegree):
        ticks_str.append("%i" % (i + 1))
        ticks[i] = i + 0.5 
    yticks(ticks, ticks_str)
    ylabel('degree', size = 9)

    # per k-mer figures
    fontdict = dict(family = "cursive", weight = "bold", size = 12, y = 1.04) 

    # 1-mers
    #subplot(3, 2, 2)
    #title('1-mer Positional Importance', fontdict)
    #pcolor(poim[0], shading = 'flat') 
    #ticks_str = ['A', 'C', 'G', 'T'] 
    #ticks = [0.5, 1.5, 2.5, 3.5]
    #yticks(ticks, ticks_str, size = 5)
    #axis([0, max_len, 0, 4])

    # 2-mers
    subplot(1, 2, 2)
    title('2-mer Positional Importance', fontdict)
    pcolor(poim[1], shading = 'flat') 
    i = 0 
    ticks = [] 
    ticks_str = [] 
    for l1 in ['A', 'C', 'G', 'T']:
        for l2 in ['A', 'C', 'G', 'T']:
            ticks_str.append(l1 + l2) 
            ticks.append(0.5 + i) 
            i += 1 
    yticks(ticks, ticks_str, fontsize = 9)
    axis([0, max_len, 0, 16])

    # 3-mers
    #subplot(3, 2, 6)
    #title('3-mer Positional Importance', fontdict)

    #if poimdegree > 2:
    #    pcolor(poim[2], shading = 'flat')
    #    i = 0
    #    ticks = []
    #    ticks_str = []
    #    for l1 in ['A', 'C', 'G', 'T']:
    #        for l2 in ['A', 'C', 'G', 'T']:
    #            for l3 in ['A', 'C', 'G', 'T']:
    #                if numpy.mod(i, 4) == 0:
    #                    ticks_str.append(l1 + l2 + l3) 
    #                    ticks.append(0.5 + i) 
    #                i += 1 
    #    yticks(ticks, ticks_str, fontsize = 5)
    #    axis([0, max_len, 0, 64])

    # x-axis on last two figures
    #for plot in [5, 6]:
    #    subplot(3, 2, plot)
    xlabel('sequence position', size = 9)


    # finishing up
    for plot in xrange(1, 3): # 6):
        subplot(1, 2, plot)
        xticks(fontsize = 9)

    for plot in [1]: #, 3, 5]:
        subplot(1, 2, plot)
        yticks(fontsize = 9)

    

    # write to file
    warnings.filterwarnings('ignore', 'Could not match*')
    #savefig('/tmp/temppylabfig.png')
    savefig(poimfile)
    clf()

def getPartitionedSet(total, crossval_repeat):
    """Partitions a number of samples into crossvalidation bins"""
    size = int(total / crossval_repeat)
    mod = total % crossval_repeat
    
    splits = []
    for i in range(0, crossval_repeat):
        if i < mod:
            splits.append(size + 1)
        else:
            splits.append(size)

    random.seed()
    ipartition = random.sample(xrange(0, total), total) # random sampling
    
    index = 0
    partitions = []

    for size in splits:
        partitions.append(ipartition[index:index + size])
        index += size
    
    return partitions
    
def getCurrentSplit(repetition, partitions, labels, seqs):
    """Returns the correct features & labels for this partition
    for this repetition"""
    X = []; Y = []; XT = []; YT = []
    for i in range(0, len(partitions)):
        if type(seqs) == type(list([])):
            for j in range(0, len(partitions[i])):
                if repetition != i:
                    X.append(seqs[partitions[i][j]])
                    Y.append(labels[partitions[i][j]])
                else:            
                    XT.append(seqs[partitions[i][j]])
                    YT.append(labels[partitions[i][j]])
        else:
            if repetition != i:
                if len(X) == 0:
                    X = seqs.take(partitions[i], axis = 1)
                    Y = labels.take(partitions[i])
                else:
                    X = numpy.concatenate((X, seqs.take(partitions[i], axis = 1)), axis = 1)
                    Y = numpy.concatenate((Y, labels.take(partitions[i])))
            else:
                XT = seqs.take(partitions[i], axis = 1)
                YT = labels.take(partitions[i])

    return X, Y, XT, YT

def saveSVM(pickle_filename, svm, kernel):
    """Pickles a Shogun SVM object to a file by saving its settings"""
    from cPickle import Pickler
    pickle_file = open(pickle_filename, 'wb')
    pck = Pickler(pickle_file)
    pck.dump((__version__, \
              svm.get_num_support_vectors(), \
              kernel.get_name(), \
              svm.get_bias(), \
              svm.get_alphas(), \
              svm.get_support_vectors()))
    pickle_file.close()

def loadSVM(pickled_svm_filename, C, labels):
    """Loads a Shogun SVM object which was pickled by saveSVM"""
    from cPickle import Unpickler, PickleError
    from shogun.Kernel import CombinedKernel 
    pickle_file = open(pickled_svm_filename, 'rb')
    unpck = Unpickler(pickle_file)
    (version, num_sv, name, bias, alphas, svs) = unpck.load()
    if (version == __version__):
        svm = LibSVM(num_sv) # same as .create_new_model(num_sv)
        svm.set_bias(bias)
        svm.set_alphas(alphas)
        svm.set_support_vectors(svs)
        kernel = CombinedKernel() #not sure if this is even required
        kernel.set_name(name) # maybe not required
        svm.set_kernel(kernel)
    else: 
        print "File was pickled by another version of EasySVM.py or is not a kernel:"
        print "Received from ", pickled_svm_filename, ": ", version, "    expected: ", __version__
        raise PickleError
    return svm

def confusionMatrix(labels_test, labels_predicted):
    """Calculates the complete confusion matrix from true/false positives/negatives"""
    if len(labels_test) != len(labels_predicted):
        return 0
    TP = 0; FP = 0; TN = 0; FN = 0
    for i in range(0, len(labels_test)):
        if labels_test[i] == 0 or labels_predicted[i] == 0:
            return 0
        if labels_test[i] > 0:
            if labels_predicted[i] > 0: TP += 1
            else: FN += 1
        else:
            if labels_predicted[i] > 0: FP += 1
            else: TN += 1
    return (TP, TN, FP, FN)

def accuracy(output, labels_test):
    """Calculates the accurracy from true/false positives/negatives"""
    TP, TN, FP, FN = confusionMatrix(labels_test, numpy.sign(output))
    return float(TP + TN) / (TP + TN + FP + FN)

def calcROC(output, LTE):
    """Uses shogun functions to calculate the area under the ROC curve"""
    pm = PerformanceMeasures(Labels(numpy.array(LTE)), Labels(numpy.array(output)))

    auROC = pm.get_auROC()
    return auROC

def calcPRC(output, LTE):
    """Uses shogun functions to calculate the precision recall curve"""
    pm = PerformanceMeasures(Labels(numpy.array(LTE)), Labels(numpy.array(output)))

    auPRC = pm.get_auPRC()
    return auPRC

def parseRange(string):
    """Parses a dash-separated string of ints into a tuple"""
    splitarray = string.split("-")

    if len(splitarray) == 1:
        return (int(splitarray[0]), int(splitarray[0]))
    if len(splitarray) == 2:
        return (int(splitarray[0]), int(splitarray[1]))
    raise ValueError("Cannot parse range " + string)

def parseFloatList(string):
    """Parses a comma-separated string of floats into a list"""
    splitarray = string.split(",")
    float_list = []
    for elem in splitarray:
        float_list.append(float(elem))
    return float_list

def parseIntList(string):
    """Parses a comma-separated string of ints into a list"""
    splitarray = string.split(",")
    int_list = [] 
    for elem in splitarray:
        int_list.append(int(elem))
    return int_list

def parsePrediction(prediction_file):
    """Returns the original output and split assignments
    of a prediction run, from a prediction_file"""
    outputs = []
    splitassignments = []

    f = open(prediction_file)
    lines = f.readlines()
    num = 0 
    for line in lines:
        if len(line) > 0 and line[0] != '#':
            elems = line.split('\t') 
            assert(len(elems) > 1)
            assert(int(elems[0]) == num) 
            num += 1 
            if len(elems) == 2:
                outputs.append(float(elems[1]))
            else:
                assert(len(elems) == 3)
                outputs.append(float(elems[1]))
                splitassignments.append(float(elems[2]))
    f.close() 
    if len(splitassignments) == 0:
        splitassignments = None

    return (outputs, splitassignments)

def non_atcg_convert(seq, nuc_con):
    """ Converts Non ATCG characters from DNA sequence """
    
    if nuc_con == '':
        sys.stderr.write("usage: Provide a choice for non ACGT nucleotide conversion [T|A|C|G|random] at last\n")
        sys.exit(-1)
    flag = 0    
    if len(nuc_con)>1:
        if nuc_con != 'random':
            flag = 1
    else:
        if re.match(r'[^ATCG]', nuc_con, re.IGNORECASE):
            flag = 1
    if flag == 1:        
        sys.stderr.write("usage: Conversion nucleotide choice -"+ nuc_con +"- failed. pick from [T|A|C|G|random]\n")
        sys.exit(-1)
    
    nuc_con = nuc_con.upper()
    mod_seq = []
    for i in range(len(seq)):
        if re.search(r'[^ACTG]', seq[i], re.IGNORECASE):
            if nuc_con == 'RANDOM':
                nucleotide = 'ATCG'
                line = ''
                for single_nuc in seq[i]:
                    if re.match(r'[^ACGT]', single_nuc, re.IGNORECASE):
                        single = random.choice(nucleotide)
                        line += single
                    else:
                        single_nuc = single_nuc.upper()
                        line += single_nuc
                mod_seq.append(line)       
            else:
                seq[i] = re.sub(r'[^ATCG|actg]', nuc_con, seq[i])
                seq[i] = seq[i].upper()
                mod_seq.append(seq[i])
        else:
            seq[i] = seq[i].upper()
            mod_seq.append(seq[i])
    return mod_seq

def non_aminoacid_converter(seq, amino_con):
    """ Converts Non amino acid characters from protein sequence """  

    if amino_con == '':
        sys.stderr.write("usage: Provide a choice for replacing non amino acid characters\n")
        sys.exit(-1)
    flag = 0
    if len(amino_con)>1:
        if amino_con != 'random':
            flag = 1
    else:        
        if re.match(r'[^GPAVLIMCFYWHKRQNEDST]', amino_con, re.IGNORECASE):
            flag = 1
    if flag == 1:        
        sys.stderr.write("usage: Replace aminoacid chioce -"+ amino_con +"- failed. Pick a valid aminoacid single letter code/random\n")
        sys.exit(-1)

    amino_con = amino_con.upper()
    opt_seq = []
    for i in range(len(seq)):
        if re.search(r'[^GPAVLIMCFYWHKRQNEDST]', seq[i], re.IGNORECASE):
            if amino_con == 'RANDOM':
                aminoacid = 'GPAVLIMCFYWHKRQNEDST'
                line = ''
                for single_amino in seq[i]:
                    if re.match(r'[^GPAVLIMCFYWHKRQNEDST]', single_amino, re.IGNORECASE):
                        r_amino = random.choice(aminoacid)
                        line += r_amino
                    else:
                        single_amino = single_amino.upper()
                        line += single_amino
                opt_seq.append(line)         
            else:
                seq[i] = re.sub(r'[^GPAVLIMCFYWHKRQNEDST|gpavlimcfywhkrqnedst]', amino_con, seq[i])
                seq[i] = seq[i].upper()
                opt_seq.append(seq[i])
        else:
            seq[i] = seq[i].upper()
            opt_seq.append(seq[i])
    return opt_seq 



def evaluate(predictions, splitassignments, labels, roc_fname = None, prc_fname = None):
    """Evaluate prediction results"""
    res_str = ""
    xval = 1
    if splitassignments != None:
        for split in splitassignments:
            if split + 1 > xval:
                xval = int(split + 1)
    if xval > 1:
        res_str = "Evaluating on %i examples in %i splits\n" % (len(labels), xval)
    else:
        res_str = "Evaluating on %i examples\n" % len(labels)

    output_splits = xval * [[]]
    label_splits = xval * [[]]
    for i in xrange(xval):
        label_splits[i] = [] 
        output_splits[i] = [] 

    for i in xrange(0, len(labels)):
        if xval > 1:
            split = int(splitassignments[i])
        else:
            split = 0
        output_splits[split].append(predictions[i])
        label_splits[split].append(labels[i])

    sum_accuracy = 0.0
    sum_roc = 0.0
    sum_prc = 0.0

    for split in xrange(xval):
        if xval > 1:
            res_str += 'Split %d\n' % (split + 1)

        LTE = label_splits[split] 
        svmout = output_splits[split]

        numpos = 0 
        for l in LTE:
            if l == 1:
                numpos += 1 
        istwoclass = numpos > 0 and numpos < len(LTE)
        if xval > 1:
            res_str += '   number of positive examples = %i\n' % numpos
        if xval > 1:
            res_str += '   number of negative examples = %i\n' % (len(LTE)-numpos)
        if istwoclass:
            auROC = calcROC(svmout, LTE)
            if xval > 1:
                res_str += '   Area under ROC curve        = %2.1f %%\n' % (100.0 * auROC)
            sum_roc += auROC
            if roc_fname != None:
                if split != xval - 1:
                    plotROC(svmout, LTE, split == xval - 1, None, "ROC curve of SVM, split %i" % (split + 1))
                else:
                    plotROC(svmout, LTE, split == xval - 1, roc_fname, "ROC curve of SVM, split %i" % (split + 1))
            auPRC = calcPRC(svmout, LTE)
            if xval > 1:
                res_str += '   Area under PRC curve        = %2.1f %%\n' % (100.0 * auPRC)
            sum_prc += auPRC
            if prc_fname != None:
                if split != xval - 1:
                    plotPRC(svmout, LTE, None, "PRC curve of SVM, split %i" % (split + 1))
                else:
                    plotPRC(svmout, LTE, prc_fname, "PRC curve of SVM, split %i" % (split + 1))

        acc = accuracy(svmout, LTE)
        if xval > 1:
            res_str += '   accuracy (at threshold 0)   = %2.1f %% \n' % (100.0 * acc)
        sum_accuracy += acc

    numpos = 0 
    for l in labels:
        if l == 1:
            numpos += 1 

    mean_roc = sum_roc/xval
    mean_prc = sum_prc/xval
    mean_acc = sum_accuracy/xval

    res_str += 'Averages\n'
    res_str += '   number of positive examples = %i\n' % round(numpos/xval)
    res_str += '   number of negative examples = %i\n' % round((len(labels) - numpos)/xval)
    res_str += '   Area under ROC curve        = %2.1f %%\n' % (100.0 * mean_roc) 
    res_str += '   Area under PRC curve        = %2.1f %%\n' % (100.0 * mean_prc) 
    res_str += '   accuracy (at threshold 0)   = %2.1f %% \n' % (100.0 * mean_acc) 

    return (res_str, mean_roc, mean_prc, mean_acc)