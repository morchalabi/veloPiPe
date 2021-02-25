#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 13:18:50 2021

@author: zfd297
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 00:13:49 2021

@author: zfd297
"""
import sys
import os
import loompy
from scipy.io import mmread
import pandas as pd
import numpy as np
from argparse import ArgumentParser

def init_args():
    """
    init_args: parse args
    Args:
        None
    Returns:
        args
    Raises:
        None
    """
    arg_parser = ArgumentParser(description="export data from database")
 
    arg_parser.add_argument("-s", "--spliced", dest="spliced", required = True, \
                            help = "local or hadoop")
    arg_parser.add_argument("-u", "--unspliced", dest="unspliced", required = True, \
                            help = "input file path")
        
    arg_parser.add_argument("-a", "--ambi", dest="ambi", required = True, \
                            help = "input file path")
    
    arg_parser.add_argument("-r", "--gene.tsv", dest="gene", required = True, \
                            help = "input file path")
        
    arg_parser.add_argument("-c", "--barcode", dest="barcode", required = True, \
                            help = "input file path")

    arg_parser.add_argument("-o", "--output", dest="output", required = True, \
                            help = "input file path")
    args = None
    try:
        args = arg_parser.parse_args()
    except Exception as e:
        logging.fatal(str(e))
    return args

def create_loom(sp_mat,unsp_mat,amb_mat, gene, barcode, output):
    
    # path_1 = "./mtx_files/" + folder_name + "/" + "spliced.mtx"
    # path_2 = "./mtx_files/" + folder_name + "/" + "unspliced.mtx"
    # path_3 = "./mtx_files/" + folder_name + "/" + "genes.tsv"
    # bc_path = "./mtx_files/" + folder_name + "/" + "barcodes.tsv"

    isExists=os.path.exists(output)
    if not isExists:
        os.makedirs(output) 
    else:
        pass
    
    spliced = mmread(sp_mat).A
    unspliced = mmread(unsp_mat).A
    abi = mmread(amb_mat).A
    blank_mat = np.zeros((spliced.shape[0], spliced.shape[1]))
    barcode = pd.read_csv(barcode, header=None)
    cell_id = barcode[0]
    
    gene_id = pd.read_csv(gene, header=None)
    gene_id = np.array(gene_id)
    gene_id = [i[0] for i in  gene_id]
    
    gene_id = np.array(gene_id)
    cell_id = np.array(cell_id)
    gene_id = gene_id.astype(str)
    cell_id = cell_id.astype(str)
    
    dicts = {'':blank_mat,'spliced':spliced, 'unspliced':unspliced, 'ambiguous':abi}
#
#
    row_attrs = {'Gene':gene_id}
#
    col_attrs = {'CellID':cell_id}
#
    loompy.create(output+"out.loom",dicts, row_attrs, col_attrs)

args = init_args()
sp_mat = args.spliced
unsp_mat = args.unspliced
amb_mat = args.ambi
gene = args.gene
barcode = args.barcode
output = args.output
create_loom(sp_mat,unsp_mat,amb_mat, gene, barcode,output)