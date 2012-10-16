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
@summary: Traverses given tree structure and outputs command line calls

"""


import taxonomy


def train(tax):
    """
    hierarchical training procedure 
    """

        
    root = tax

    print root
    
    grey_nodes = [root]
    
    
    #####################################################
    #         top-down processing of taxonomy           #
    #####################################################
    
    while len(grey_nodes)>0:
    
        node = grey_nodes.pop(0) # pop first item
    
        # enqueue children
        if node.children != None:
            grey_nodes.extend(node.children)
                
        # get data below current node
        data_keys = node.get_data_keys()
        
        # set up presvm
        if node.is_root():
            # no parent at root node
            parent_svm = "none" 
        
        else:
            # regularize against parent predictor
            parent_svm = node.parent.name

        #TODO: invoke system calls
        print node.name, data_keys, parent_svm


def main():
    
    tax = taxonomy.create_tax_four()
    
    train(tax)
    
    
if __name__ == '__main__':
    main()
    
    