#!/usr/bin/env python
"""
Program to plot the difference in gene expression based on gene 
families in A. thaliana.

Using RNA-Seq data, programs like rDiff or DESeq will generate 
the differential gene expression scores. Based on different gene
families of an organism, the distribution of gene expression may
also vary. Here this program consider three main types: expressed, 
non-expressed and differentially expressed.

Usage: 
python diff_genes.py in.gene_family_name in.rdiff_out 
based on the number of replicates, the gene distribution plots 
will be created in SVG format.

This script requires:
matplotlib: http://matplotlib.sourceforge.net/ 
scipy: http://www.scipy.org/
"""
from __future__ import division
import sys, re
import matplotlib.pyplot as plt 
from matplotlib import mpl 
import scipy as sp 

def	__main__():
    try:
        gene_fam_fname=sys.argv[1]
        diff_exp_fname=sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
    # get family name for genes. defined dataset for organism. 
    genes=get_family_name(gene_fam_fname)
    # get differential expression data 
    gids, exp_score=get_diff_exp(diff_exp_fname)
    # calculate q value(fdr) for each replicates.
    for rep_nb, rep_sc in sorted(exp_score.items()):
        rep_pv=sp.array(rep_sc)
        qval=estimate_q_values(rep_pv)
        fam_name=grading_genes(genes, gids, qval, rep_pv)
        # plot the gene family display distribution 
        plot_gene_class(fam_name, rep_nb)

def plot_gene_class(family_maps, prefix_out):

    maps=sorted(family_maps.keys())
    nmaps=len(maps)+1
    # TODO setting a legend for the three different class of genes.
    fig=plt.figure(figsize=(12,10))
    fig.subplots_adjust(top=0.94, bottom=0.06, left=0.2, right=0.9)
    # setting plotting structure 
    color_code=['blue', 'green', 'grey']
    for i,m in enumerate(maps):
        dis_nam=re.sub(r'_', r' ', m)
        gcnt=0 # n 
        for keys in ['dif_exp', 'mod_exp', 'not_exp']:
            try:cnt=family_maps[m][keys]
            except:cnt=0
            gcnt+=cnt
        xq, previous=[0], 0.0
        for keys in ['dif_exp', 'mod_exp', 'not_exp']:
            try:element=family_maps[m][keys]
            except:element=0
            xq.append(round(element/gcnt, 2)+previous)
            previous=xq[-1]
        ax=plt.subplot(nmaps, 1, i+1)
        plt.axis("off")
        cmap = mpl.colors.ListedColormap(color_code)
        cmap.set_over('0.15')
        cmap.set_under('0.75')
        norm = mpl.colors.BoundaryNorm(xq, cmap.N)
        cb2=mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal', spacing='proportional', boundaries=xq)
        pos=list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], dis_nam, fontsize=9, horizontalalignment='right')
        fig.text(0.91, pos[1], str(gcnt), fontsize=9, horizontalalignment='left')
        if m==maps[-1]:
            nmb=[]
            color_code=['white']
            cmap = mpl.colors.ListedColormap(color_code)
            for nb in range(0,11,1):
                nmb.append(nb/10)
            ax=plt.subplot(nmaps, 1, i+2)
            cb2=mpl.colorbar.ColorbarBase(ax, cmap=cmap, orientation='horizontal', boundaries=nmb)
    plt.suptitle('Genes expression/diff. expression')
    plt.savefig(str(prefix_out)+'_replicate_gene_distribution.svg')

def grading_genes(gene2family, gene_ids, qvalues, pvalues):
    """Classify genes according to the q value and check for which family it belongs.
    """
    family_class=dict()
    for i,fid in enumerate(gene_ids):
        if qvalues[i]<0.05: # 5% threshold value 
            if fid in gene2family:
                for family in gene2family[fid]:
                    if family in family_class:
                        if 'dif_exp' in family_class[family].keys():
                            family_class[family]['dif_exp']+=1
                        else:   
                            family_class[family]['dif_exp']=1
                    else:
                        family_class[family]={'dif_exp':1}
            else:
                if 'other' in family_class:
                    if 'dif_exp' in family_class['other'].keys():
                        family_class['other']['dif_exp']+=1
                    else:
                        family_class['other']['dif_exp']=1
                else:
                    family_class['other']={'dif_exp':1}
        elif qvalues[i]>0.05: # expressed or not-expressed
            if pvalues[i]<0.05: # based on p-value accept alternative hypothesis:- Ho:not expressed, Ha:expressed
                if fid in gene2family:
                    for family in gene2family[fid]:
                        if family in family_class:
                            if 'mod_exp' in family_class[family].keys():
                                family_class[family]['mod_exp']+=1
                            else:
                                family_class[family]['mod_exp']=1
                        else:
                            family_class[family]={'mod_exp':1}
                else:
                    if 'other' in family_class:
                        if 'mod_exp' in family_class['other'].keys():
                            family_class['other']['mod_exp']+=1
                        else:
                            family_class['other']['mod_exp']=1
                    else:
                        family_class['other']={'mod_exp':1}
            elif pvalues[i]>0.05: # not-expressed
                if fid in gene2family:
                    for family in gene2family[fid]:
                        if family in family_class:
                            if 'not_exp' in family_class[family].keys():
                                family_class[family]['not_exp']+=1
                            else:
                                family_class[family]['not_exp']=1
                        else:
                            family_class[family]={'not_exp':1}
                else:
                    if 'other' in family_class:
                        if 'not_exp' in family_class['other'].keys():
                            family_class['other']['not_exp']+=1
                        else:
                            family_class['other']['not_exp']=1
                    else:
                        family_class['other']={'not_exp':1}
    return family_class

def estimate_q_values(PV,m=None,pi=1):
    """estimate q vlaues from a list of Pvalues
    this algorithm is taken from Storey, significance testing for genomic ...
    m: number of tests, (if not len(PV)), pi: fraction of expected true null (1 is a conservative estimate)
    originally written by Oliver Stegel from MPI and edited by Vipin
    """
    if m is None:
        m = len(PV)
    lPV = len(PV)
    #1. sort pvalues
    PV = PV.squeeze()
    IPV = PV.argsort()
    PV  = PV[IPV]
    #2. estimate lambda
    if pi is None:
        lrange = sp.linspace(0.05,0.95,max(lPV/100,10))
        pil    = sp.double((PV[:,SP.newaxis]>lrange).sum(axis=0))/lPV
        pilr   = pil/(1-lrange)
        #ok, I think for SNPs this is pretty useless, pi is close to 1!
        pi =1
        #if there is something useful in there use the something close to 1
        if pilr[-1]<1:
            pi = pilr[-1]
    #3. initialise q values
    QV_ = pi * m/lPV* PV
    #4. update estimate
    for i in xrange(lPV-2,0,-1):
        QV_[i] = min(pi*m*PV[i]/(i+1),QV_[i+1])
    #5. inverst sorting
    QV = sp.zeros_like(PV)
    QV[IPV] = QV_
    return QV

def get_diff_exp(diff_fname):
    """Parse the output generated by rDiff/DESeq program. 
    """
    gid, exp_rate=[], dict()
    dfh=open(diff_fname, "rU")
    for line in dfh:
        line=line.strip('\n\r').split('\t')
        try:
            float(line[1])
        except:continue
        gid.append(line[0])
        for i, pscore in enumerate(line[1:]):
            if i in exp_rate:
                exp_rate[i].append(float(pscore))
            else:
                exp_rate[i]=[float(pscore)]
    dfh.close()
    return gid, exp_rate
    
def get_family_name(fname):
    """Parse gene family name.
    """
    gene_ids=dict()
    famh=open(fname, 'rU')
    for line in famh:
        gid, fam_name=line.strip("\n\r").split('\t')
        if gid in gene_ids:
            gene_ids[gid].append(fam_name)
        else:
            gene_ids[gid]=[fam_name]
    famh.close()
    return gene_ids

if __name__=="__main__":
	__main__()
