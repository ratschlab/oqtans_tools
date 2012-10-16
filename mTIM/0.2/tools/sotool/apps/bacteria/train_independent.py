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

import os
import taxonomy
import pythongrid as pg
import system_call


def train(tax, tax_name, num_train):
    """
    hierarchical training procedure
    """

    datadir = "../../data"
    num_vald = num_train

    params_c = "[1 10 100 250]"
    params_b = "[0 0.25 0.5 0.75 1.0]"

    # get data below current node
    data_keys = tax.get_data_keys()
 

    commands = []    
       

    for org_name in data_keys:
    
                   
        # create org list in matlab format
        org_list = "{'%s'}" % org_name
        
       
        # determine current output directory
        current_out = "../out/%s/indep/%s" % (tax_name, org_name)

        # create dirs if not present
        if not os.path.exists(current_out):
            #Super-mkdir; create a leaf directory and all intermediate ones.
            print "creating dir %s" % (current_out)
            os.makedirs(current_out)


        # construct commands
        matlab_call_train = "matlab -nojvm -nosplash -nodesktop -r \"train_orgs('%s',%s,'%s',%i,%i,%s,%s,'%s'), exit;\" " % (datadir, org_list, current_out, num_train, num_vald, params_c, params_b, "none")
        matlab_call_eval = "matlab -nojvm -nosplash -nodesktop -r \"evaluate_result('%s'), exit;\" " % ("../" + current_out)

        command_tuple = matlab_call_train, matlab_call_eval
        commands.append(command_tuple)

        print matlab_call_train
        print matlab_call_eval
        # example call:
        # matlab -nojvm -nosplash -nodesktop -r train_orgs('../../data',{'Bacillus_anthracis_Ames_uid57909', 'Bacillus_subtilis_168_uid57675'},'../out/tax1/bacillus',10,50,[0.1 1 10],[0 0.1 0.5 0.8 0.9],'../../out/All_Tasks/evaluation.mat'), exit;

    result = pg.map(system_call.call2, commands, local=True, maxNumThreads=4, mem="16G")


def main():
    """
    choose taxonomy and output dir name and invoke training procedure
    """

    # TO USE DIFFERENT TAXONOMY CHANGE THIS    
    tax = taxonomy.create_tax_four_distant()

    num_train = 100

    # come up with better naming convention
    tax_name = "tax4distant_" + str(num_train)
    
    train(tax, tax_name, num_train)

    print "training done. exiting."
    
    
if __name__ == '__main__':
    main()
    
