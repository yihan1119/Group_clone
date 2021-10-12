#!/usr/bin/env python3
import os
import re
import sys
print(sys.executable, file=sys.stderr)
print(sys.path, file=sys.stderr)
print(sys.version, file=sys.stderr)

import collections
import argparse
#import tables
import itertools
import matplotlib
import glob
import math


import scipy.io
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import scipy.sparse as sp_sparse
import scanpy as sc 
import scanpy.external as sce

from collections import defaultdict
from scipy import sparse, io

import scanpy.external as sce
import matplotlib

from scipy.sparse import csr_matrix
from multiprocessing import Pool
#from matplotlib_venn import venn2, venn2_circles
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

print('numpy', np.__version__, file=sys.stderr)
print('pandas', pd.__version__, file=sys.stderr)
print('scipy', scipy.__version__, file=sys.stderr)
print('matplotlib', matplotlib.__version__, file=sys.stderr)
print('scanpy', sc.__version__, file=sys.stderr)

def load_clone_tree(clone_dict_file):
    Clone_dict = {}
    with open(clone_dict_file) as f:
        first_line = f.readline()
        word = '}'
        for line in f:
            if not word in line:
                ID, clones = line.split(":")
                clone_IDs = clones.replace("'[", "").replace("]',", "").replace("'", "").strip(' \n')
                individual_clone_ID = clone_IDs.split(', ')
                Clone_dict.update({ID.strip("'") : individual_clone_ID})
    return Clone_dict


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-c', '--clone_tree', dest='clone_tree', required=True,
        type=str,
        help='specify the text file of clonal cells tree.'
    )

    parser.add_argument(
        '-s', '--sgRNA_df', dest='sgRNA_df', required=True,
        type=str,
        help='specify the pickle file of processed sgRNA df.'
    )

    parser.add_argument(
        '-n', '--number', dest='number', required=True,
        type=int,
        help='specify the number of random selection cells'
    )    

    parser.add_argument(
        '-g', '--group', dest='group', required=True,
        type=int,
        help='specify the group of cells to analyze: clone or non-clone'
    )

    parser.add_argument(
        '-o', '--output', dest='output', required=True,
        type=str,
        help='specify the output directory'
    )  
    args = parser.parse_args()
    SGRNA_DF = args.sgRNA_df
    Clone_dict_file = args.clone_tree
    GROUP = args.group
    NUMBER = args.number
    OUTPUT_FILE = args.output

    #check the normalization method
    if (GROUP != 'clone') and (GROUP != 'non-clone'):
        print("Incorrect groups of cells. Has to be either 'clone' or 'non-clone'.", file = sys.stderr, flush=True)
        sys.exit(0)


    #Load data     
    Clone_dict = load_clone_tree(Clone_dict_file)
    sgrna_df = pd.read_pickle(SGRNA_DF)
    sgrna_df_adj_bool = sgrna_df > 0
    sgrna_num, cell_num = sgrna_df_adj_bool.shape
    #remove cells without sgRNA 
    cells_with_sgRNA = sgrna_df_adj_bool.T[(np.sum(sgrna_df_adj_bool, axis=0) > 0).values].index

    #Flatten all the clones 
    All_cell_ID_list = []
    if (GROUP == 'clone' ):
        for i in Clone_dict.keys():
            if len(Clone_dict[i]) > 100:
                cell_ID_list = Clone_dict[i]
                for k in cell_ID_list:
                    All_cell_ID_list.append(k)
    
    elif (GROUP == 'non-clone'):
        for i in Clone_dict.keys():
            if len(Clone_dict[i]) < 2:
                cell_ID_list = Clone_dict[i]
                for k in cell_ID_list:
                    All_cell_ID_list.append(k)
    
    #randomly select cells 
    All_sgRNA_cell = set(All_cell_ID_list).intersection(set(cells_with_sgRNA))
    Random_cell_ID_list = np.random.choice(All_sgRNA_cell, size=NUMBER, replace=False)
    All_sgrna_overlap_df = pd.DataFrame(data=None, index=Random_cell_ID_list, columns=Random_cell_ID_list)
    print('The size of randomly selected cells overlap rate df is: ' + str(All_sgrna_overlap_df.shape), file=sys.stderr, flush=True)

    #find the sgRNAs in each randomly selected cells 
    All_clone_cells_sgRNA_dict = defaultdict(list)
    for i in Random_cell_ID_list:
        clonal_sgrna = list(sgrna_df_adj_bool.index[(sgrna_df_adj_bool[i] > 0).values])
        All_clone_cells_sgRNA_dict[i].append(clonal_sgrna)

    counter = 0
    for i in range(len(Random_cell_ID_list)):
        cell1_region = set(All_clone_cells_sgRNA_dict[All_sgrna_overlap_df.columns[i]][0])
        for j in range(len(Random_cell_ID_list)):
            cell2_region = set(All_clone_cells_sgRNA_dict[All_sgrna_overlap_df.columns[j]][0])
            if math.isnan(All_sgrna_overlap_df.iloc[j, i]) == True:
                pval = stats.hypergeom.sf(len(cell1_region.intersection(cell2_region))-1, sgrna_num, len(cell1_region), len(cell2_region))
                
                All_sgrna_overlap_df.iloc[j, i] = pval
                All_sgrna_overlap_df.iloc[i, j] = pval
                
                if counter % 10000 == 0:
                    print(counter, file=sys.stderr, flush=True)
                elif counter % 100 == 0:
                    print(".", end = '', file=sys.stderr, flush=True)
                counter += 1
            else:
                continue
        
    All_sgrna_overlap_df.to_pickle(OUTPUT_FILE)

if __name__ == '__main__':
    main()
