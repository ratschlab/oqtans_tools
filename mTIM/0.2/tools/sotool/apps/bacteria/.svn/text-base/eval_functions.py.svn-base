#!/usr/bin/env python2.6
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Written (W) 2011 Christian Widmer
# Copyright (C) 2011 Max-Planck-Society

"""
Created on 27.05.2011
@author: Christian Widmer
@summary: Evaluate predictions using three strategies

sensitivity = TP / (all pos) = TP / (FN + TP)
specificity = TP / (all pos predictions) = TP / (FP + TP)
f-score = 
"""

import sys


def eval_nucleotid_level(pred_seq, true_seq):
    """
    evaluate sensitivity and specificity on nucleotide level
    
    """
    
    assert len(pred_seq) == len(true_seq), "lengths don't match: %i!=%i" % (len(pred_seq), len(true_seq))
    
    tp = 1
    tn = 1
    fp = 1
    fn = 1
    
    # iterate over sequences
    for (idx, pred_n) in enumerate(pred_seq):
        
        true_n = true_seq[idx]
        
        # fill confusion matrix
        if pred_n == "1" and true_n == "1":
            tp += 1
        elif pred_n == "0" and true_n == "0":
            tn += 1
        elif pred_n == "1" and true_n == "0":
            fp += 1
        elif pred_n == "0" and true_n == "1":
            fn += 1
            
    sensitivity = compute_sensitivity(tp, fn)
    specificity = compute_sensitivity(tp, fn)
    
    fscore = compute_fscore(sensitivity, specificity)
    
    return fscore



def eval_start(pred_seq, true_seq):
    """
    evaluate sensitivity and specificity on nucleotide level
    
    """
    
    assert len(pred_seq) == len(true_seq), "lengths don't match: %i!=%i" % (len(pred_seq), len(true_seq))
    
    tp = 1
    tn = 1
    fp = 1
    fn = 1
    
    pred_prev_n = 0
    true_prev_n = 0
    
    
    # iterate over sequences
    for (idx, pred_n) in enumerate(pred_seq):
        
        true_n = true_seq[idx]
        
        # do we have a start?
        if true_prev_n == 0 and true_n == 1:
            true_start = 1
        else:
            true_start = 0
        
        # do we have a start?
        if pred_prev_n == 0 and pred_n == 1:
            pred_start = 1
        else:
            pred_start = 0
        
        # fill confusion matrix
        if pred_start == "1" and true_start == "1":
            tp += 1
        elif pred_start == "0" and true_start == "0":
            tn += 1
        elif pred_start == "1" and true_start == "0":
            fp += 1
        elif pred_start == "0" and true_start == "1":
            fn += 1
            
        # keep previous nucleotides around
        pred_prev_n = pred_n
        true_prev_n = true_n
        
        
    sensitivity = compute_sensitivity(tp, fn)
    specificity = compute_sensitivity(tp, fn)
    
    fscore = compute_fscore(sensitivity, specificity)
    
    return fscore


def eval_stop(pred_seq, true_seq):
    """
    evaluate sensitivity and specificity on nucleotide level
    
    """
    
    assert len(pred_seq) == len(true_seq), "lengths don't match: %i!=%i" % (len(pred_seq), len(true_seq))
    
    tp = 1
    tn = 1
    fp = 1
    fn = 1
    
    pred_prev_n = 0
    true_prev_n = 0
    
    
    # iterate over sequences
    for (idx, pred_n) in enumerate(pred_seq):
        
        true_n = true_seq[idx]
        
        # do we have a start?
        if true_prev_n == 0 and true_n == 1:
            true_start = 1
        else:
            true_start = 0
        
        # do we have a start?
        if pred_prev_n == 0 and pred_n == 1:
            pred_start = 1
        else:
            pred_start = 0
        
        # fill confusion matrix
        if pred_start == "1" and true_start == "1":
            tp += 1
        elif pred_start == "0" and true_start == "0":
            tn += 1
        elif pred_start == "1" and true_start == "0":
            fp += 1
        elif pred_start == "0" and true_start == "1":
            fn += 1
            
        # keep previous nucleotides around
        pred_prev_n = pred_n
        true_prev_n = true_n
        
        
    sensitivity = compute_sensitivity(tp, fn)
    specificity = compute_sensitivity(tp, fn)
    
    fscore = compute_fscore(sensitivity, specificity)
    
    return fscore



def eval_block(pred_seq, true_seq):
    """
    evaluate sensitivity and specificity on block level
    
    
    """
    
    #TODO: think about this one
    assert len(pred_seq) == len(true_seq), "lengths don't match: %i!=%i" % (len(pred_seq), len(true_seq))
    
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    
    # iterate over sequences
    for (idx, pred_n) in enumerate(pred_seq):
        
        true_n = true_seq[idx]
        
        # fill confusion matrix
        if pred_n == "1" and true_n == "1":
            tp += 1
        elif pred_n == "0" and true_n == "0":
            tn += 1
        elif pred_n == "1" and true_n == "0":
            fp += 1
        elif pred_n == "0" and true_n == "1":
            fn += 1
            
    sensitivity = compute_sensitivity(tp, fn)
    specificity = compute_sensitivity(tp, fn)
    
    fscore = compute_fscore(sensitivity, specificity)
    
    return fscore


def compute_fscore(sensitivity, specificity):
    """
    return fscore
    
    """

    #TODO: fix
    return float(sensitivity + specificity) / 2.0


def compute_sensitivity(tp, fn):
    """
    compute sensitivity
    sensitivity = TP / (all pos) = TP / (FN + TP)
    """
    
    return float(tp) / float(tp + fn) 

    
def compute_specificity(tp, fp):
    """
    compute specificity
    specificity = TP / (all pos predictions) = TP / (FP + TP)
    """
    
    return float(tp) / float(fp, tp)



if __name__ == '__main__':
    

    s1 = "00011100"
    s2 = "01111100"
    print eval_nucleotid_level(s1, s2)
    print eval_start(s1, s2)
    print eval_stop(s1, s2)
    
    
    # parse input args
    pred_seq = sys.argv[1]
    true_seq = sys.argv[2]
    mode = sys.argv[3]
    
    if mode == "nucleotide":
        ret = eval_nucleotid_level(pred_seq, true_seq)
    elif mode == "start":
        ret = eval_start(pred_seq, true_seq)
    elif mode == "stop":
        ret = eval_stop(pred_seq, true_seq)
    elif mode == "block":
        ret = eval_block(pred_seq, true_seq)
    
    # write back to std
    sys.stdout.write(str(ret))
    