#!/usr/bin/env python3
import os
import re
import sys
print(sys.executable, file=sys.stderr, flush=True)
print(sys.path, file=sys.stderr, flush=True)
print(sys.version, file=sys.stderr, flush=True)

import collections
import argparse
import itertools
import matplotlib
import math

import scipy.io
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import scipy.sparse as sp_sparse

from collections import defaultdict
from scipy import sparse, io

print('numpy', np.__version__, file=sys.stderr, flush=True)
print('pandas', pd.__version__, file=sys.stderr, flush=True)
print('scipy', scipy.__version__, file=sys.stderr, flush=True)
print('matplotlib', matplotlib.__version__, file=sys.stderr, flush=True)
print('scanpy', sc.__version__, file=sys.stderr, flush=True)

print("start analyzing.", file=sys.stderr, flush=True)

def hypergeom_two_cell_sgrnas(array1, 
                              array2, 
                              sum_df, 
                              sgrna_num, 
                              cell_bc, 
                              k):
    
    x = np.dot(array1, array2) - 1
    M = sgrna_num
    n = sum_df[cell_bc]
    N = sum_df[k]
    p_val = stats.hypergeom.sf(x, M, n, N)
    
    return p_val


def group_clone_hypergeom(sgrna_df, 
                          cutoff):
    
    counter = 0
    clone_dict = defaultdict(list)
    clone_pval = []
    sgrna_num, cell_num = sgrna_df.shape
    sum_df = sgrna_df.sum(axis = 0)
    pval_threshold = (cutoff/((cell_num)*(cell_num)/2)) #bonferroni correction 
    
    for i in np.arange(0, cell_num):
        cell_bc = sgrna_df.columns[i]
        sgrna_pattern = sgrna_df[cell_bc].values.astype(int)
        if len(clone_dict) == 0:
            clone_dict[cell_bc].append(cell_bc)
            
        bool_found_overlap = False
        for k in clone_dict.keys():
            sgrna_pattern_uniq = sgrna_df[k].values.astype(int)
            
            pval = hypergeom_two_cell_sgrnas(sgrna_pattern,
                                             sgrna_pattern_uniq,
                                             sum_df,
                                             sgrna_num,
                                             cell_bc,
                                             k)
            
            if pval <= pval_threshold: #same clone
                bool_found_overlap = True
                clone_dict[k].append(cell_bc)
                clone_pval.append(pval)

        if (bool_found_overlap == False):
            clone_dict[cell_bc].append(cell_bc)
        
        if counter % 10000 == 0:
            print(counter, file=sys.stderr, flush=True)
        elif counter % 100 == 0:
            print(".", end = '', file=sys.stderr, flush=True)
        counter += 1
    
    return clone_dict, clone_pval



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-s', '--sgRNA_df', dest='sgRNA_df', required=True,
        type=str,
        help='specify the pickle file of processed sgRNA df.'
    )
    parser.add_argument(
        '-c', '--cutoff', dest='cutoff', required=True,
        type=float,
        help='specify the cutoff percentage (ex. 0.1 = 10%) for sgRNA overlap rate or hypergeom p-value (ex. 0.0001).'
    )
    parser.add_argument(
        '-o', '--output', dest='output', required=True,
        type=str,
        help='specify the output directory and file name.'
    )
    
    parser.add_argument(
        '-op', '--output_pval', dest='output_pval', required=True,
        type=str,
        help='specify the output directory and pvalue list file name.'
    )
    args = parser.parse_args()
    
    SGRNA_DF = args.sgRNA_df
    CUTOFF = args.cutoff
    OUTPUT_FILE = args.output
    OUTPUT_PVAL_FILE = args.output_pval
    
    
    print("start loading sgrna df.", file=sys.stderr, flush=True)
    sgrna_df = pd.read_pickle(SGRNA_DF)
    print("loaded sgrna df.", file=sys.stderr, flush=True)
    sgRNA_df_bool = sgrna_df > 0
    
    print("start calculte unique clones", file=sys.stderr, flush=True)

    clone_dict, clone_pval = group_clone_hypergeom(sgRNA_df_bool, 
                                                       CUTOFF)   
    
    #save file
    with open(OUTPUT_FILE,'w') as file:
        file.write("dictionary_name = { \n")
        for k in sorted (clone_dict.keys()):
            file.write("'%s':'%s', \n" % (k, clone_dict[k]))
        file.write("}")

    with open(OUTPUT_PVAL_FILE, "w") as f:
        for p_val in clone_pval:
            f.write(str(p_val) +"\n")

if __name__ == '__main__':
    main()
