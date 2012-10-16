#!/usr/bin/env python 
"""Wrapper script for topGO Ontology program. 

Usage:
python topgo_main.py organism_name diff_gene_exp_data.in p-value fisher_test ontology_term

Example:
python topgo_main.py Arabidopsis_thaliana data/galaxy_import/gene_exp_col0.txt 0.01 fisher MF

This program requires: 
    rpy2
"""
import sys, re 
import collections 
import rpy2.robjects as robjects 

def get_exp_data(diff_fname):
    """Parse genes with differential expression score.
       The file expects genes with a score in a tab delimited format.
       AT1G01010       0.005
    """
    gene_exp_rate=dict()
    dfh=open(diff_fname, "rU")
    for line in dfh:
        parts=line.strip('\n\r').split('\t')
        try:
            float(parts[1])
        except:continue
        gene_exp_rate[parts[0]]=float(parts[1])
    dfh.close()
    return gene_exp_rate

def get_go_map_info(genes_w_pval, go_fname):
    """Parse GO annotation file for the requested organism. 
       The file expects Gene ID and GO terms in a tab delimited format.
       AT1G01010       GO:0006355
    """
    gene2go=collections.defaultdict(list)
    go_fh=open(go_fname, 'rU')
    for line in go_fh:
        parts=line.strip('\n\r').split('\t')
        if parts[0] in genes_w_pval:
            gene2go[parts[0]].append(parts[1])
    go_fh.close()
    return dict(gene2go)

def _vector_space_R(init_dict):
    """Create a named vector in R"""
    return robjects.r.c(**init_dict)

def exe_topGO(total_genes, gene_to_go, gene_pval, st_test, go_term):
    """Wrap the R execution in python with rpy2.
    """
    robjects.r('''library(topGO)''')
    robjects.r('''topDiffgenes = function(allScore) {return (allScore < %s)}''' % gene_pval)
    params = {"ontology" : go_term,
                "annot" : robjects.r["annFUN.gene2GO"],
                "geneSelectionFun" : robjects.r["topDiffgenes"],
                "allGenes" : _vector_space_R(total_genes),
                "gene2GO" : _vector_space_R(gene_to_go)}
    go_data = robjects.r.new("topGOdata", **params)
    results_fis = robjects.r.runTest(go_data, algorithm="classic",
                statistic=st_test)
    scores = robjects.r.score(results_fis)
    table_entries = min(10, len(scores.names)) ## setting an upper limit to 10 nodes 
    #results_table = robjects.r.GenTable(go_data, classic=results,
    #             topNodes=table_entries)
    #print results_table
    robjects.r('''library(Rgraphviz)''')
    robjects.r.printGraph(go_data, results_fis, firstSigNodes=table_entries, useInfo = "def", fn="resultFile", pdfSW="TRUE")

def __main__():
    """main function 
    """
    try:
        organism = sys.argv[1]
        diff_exp_fname = sys.argv[2]
        thd_pval = float(sys.argv[3])
        test_type = sys.argv[4]
        ontology_term = sys.argv[5] ## molecular function (MF)
    except:
        print __doc__
        sys.exit(-1)

    ## ? GO term to Gene ID Mapping file for organisms.
    go_map_org=dict(Arabidopsis_thaliana="/mnt/galaxyTools/tools/oqtans/topGO/data/ATH_GO_TAIR.txt",
		Drosophila_melanogaster="/mnt/galaxyTools/tools/oqtans/topGO/data/DM_GO_FB.txt",
		Mus_musculus="/mnt/galaxyTools/tools/oqtans/topGO/data/MM_GO_MGI.txt",
                Saccharomyces_cerevisiae="/mnt/galaxyTools/tools/oqtans/topGO/data/SCV_GO_SGD.txt")

    ## ? get the geneList based on differential expression data file 
    all_genes=get_exp_data(diff_exp_fname)

    ## ? GO mapping file (Custom GO annotations) 
    gene_to_go=get_go_map_info(all_genes, go_map_org[organism])
    if len(gene_to_go)==0:
        raise ValueError("No GO terms match to provide genes list, check for gene identifier mix-up.")
    
    ## ? run topGO 
    exe_topGO(all_genes, gene_to_go, thd_pval, test_type, ontology_term)

if __name__=="__main__":
    __main__()
