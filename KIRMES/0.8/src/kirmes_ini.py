"""Short configuration file. 
Most parameters can be overridden 
using command line parameters
for kirmes.py"""

# Filenames
POSITIVES_FILENAME = "positives.fasta"
NEGATIVES_FILENAME = "negatives.fasta"
NGFF_FILENAME = "negatives.gff"
PGFF_FILENAME = "positives.gff"
QUERY_FILENAME = "query.fasta"
QGFF_FILENAME = "query.gff"
IMM_FILENAME = "matrix.imm"
KERNEL_FILENAME = "kirmes_resume.pck"
CONSERVATION_FILENAME = "cons.bed"
OUTPUT_FILENAME = "result.html"

# Program Parameters, all command-line overridable
CROSSVALIDATION_ROUNDS = 5
RBF_SIGMA = 5
MOTIF_WINDOW_WIDTH = 20
MOTIF_LENGTH = 6
NOF_MOTIFS = 50
SAMPLING_STRATEGY = "kkmc"
REPLACEMENT_CHARACTER = "A"

#Shogun SVM Kernel Parameters
K_RBF_WIDTH = 0.5
K_RBF_C = 10.0

K_WDS_WIDTH = 9.0
K_WDS_C = 90.0
K_SEQLENGTH = MOTIF_WINDOW_WIDTH
K_POIMDEGREE = 6
K_SHIFT = 10
K_DEGREE = 10


WDS_KERNEL_PARAMETERS = {}
WDS_KERNEL_PARAMETERS['name'] = 'wd'
WDS_KERNEL_PARAMETERS['shift'] = K_SHIFT
WDS_KERNEL_PARAMETERS['width'] = K_WDS_WIDTH
WDS_KERNEL_PARAMETERS['C'] = K_WDS_C
WDS_KERNEL_PARAMETERS['seqlength'] = K_SEQLENGTH
WDS_KERNEL_PARAMETERS['degree'] = K_DEGREE

RBF_KERNEL_PARAMETERS = {}
RBF_KERNEL_PARAMETERS['name'] = 'gauss'
RBF_KERNEL_PARAMETERS['width'] = K_RBF_WIDTH
RBF_KERNEL_PARAMETERS['C'] = K_RBF_C

K_COMBINED_C = 70
K_CACHE_SIZE = 500
K_NUM_THREADS = 6