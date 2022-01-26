#!/usr/bin/env python3
import os
import re
import sys
import collections
import argparse
import tables
import itertools
import matplotlib
import numba
#%matplotlib inline


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
from multiprocessing import Pool

np.random.seed(0)


def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))
def turn_point(sgRNA_name, df):
    sgRNA_count  = df.T.filter(items=[sgRNA_name]).sum(axis=1).sort_values(ascending=False)
    sgRNA_cumsum = sgRNA_count.cumsum()

    #get the total cell number of this sgRNA
    cell_num = np.argwhere(sgRNA_count > 0).size

    #calculate the turning point by using the max derivative
    turning_point = sgRNA_cumsum.loc[((sgRNA_cumsum.diff()) / sgRNA_count.sum() > (1/cell_num))].shape
    
    return(sgRNA_count.iloc[turning_point])

#load sgRNA-barcode data
def load_data(input_file):
    input_fh  = open(input_file, 'r')
    #define parameters
    cell_bc_list   = []
    num_sgRNA_list = np.array([])
    sgRNAs         = []
    umis           = []
    
    #initiate a 2D dictionary
    data_dict = nested_dict(2, list)
    for line in input_fh:
        cell_bc    = line.strip().split('\t')[0] + '-1'
        num_sgRNA  = line.strip().split('\t')[2]
        sgRNA_list = line.strip().split('\t')[3].split(';')
        umi_list   = line.strip().split('\t')[5].split(';')
        full_cell_bc = cell_bc
        for i in zip(sgRNA_list, umi_list):
            data_dict[full_cell_bc][i[0]] = i[1]
    
    #read the 2D dictionary into the pandas DataFrame
    df = pd.DataFrame(data_dict).fillna(0).astype(int)    
    return df

#filter the sgRNA UMI count based on the cutoff values
def filter_umi (df_input,replace=True):
    if replace:
        df = df_input.copy()
    else:
        df = df_input
    #calculate the cutoff value for each sgRNA in the dataset
    sgRNA_cutoff = [turn_point(i, df) for i in list(df.index)]
    for i in range(0, len(sgRNA_cutoff)):
        df.iloc[i, :].loc[df.iloc[i, :] <= sgRNA_cutoff[i]] = 0
    return df


def filter_umi_range (df_input, cutoff_values, replace=True):
    if replace:
        df = df_input.copy()
    else:
        df = df_input
    #calculate the cutoff value for each sgRNA in the dataset  
    num_sgRNA = df.shape[0]
    sgRNA_cutoff = [cutoff_values] * num_sgRNA
    for i in range(0, len(sgRNA_cutoff)):
        df.iloc[i, :].loc[df.iloc[i, :] <= sgRNA_cutoff[i]] = 0  
    return df

def rename_BC(df, number):
    BC_list = []
    for i in df.columns:
        i = i.strip('-1')
        new_BC = i + '-' + str(number)
        BC_list.append(new_BC)

    df.columns = BC_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input_dir', dest='input_dir', required=True,
        type=str,
        help='specify the input from regex output'
    )

    parser.add_argument(
        '-m', '--method', dest='method', required=True,
        type=str,
        help='specify the method for filter umi: curve or fix_range'
        )

    parser.add_argument(
        '-f', '--cutoff', dest='cutoff', required=False,
        type=int,
        help='specify the cutoff values for filter UMI fix range analysis, has to be integer.'
        )

    parser.add_argument(
        '-n', '--number', dest='number', required=True,
        type=int,
        help='specify the number of sgrna libraries to combine'
    )

    parser.add_argument(
        '-o', '--output_dir', dest='output_dir', required=True,
        type=str,
        help='specific the output directory and file name'
    )

    args = parser.parse_args()

    INPUT_FILE = args.input_dir
    NUMBER = args.number
    METHOD = args.method
    CUTOFF = args.cutoff
    OUTPUT_FILE = args.output_dir
    
    if (METHOD != 'curve') and (METHOD != 'fix_range'):
        print("Incorrect filter UMI method. Has to be either 'curve' or 'fix_range'.", file = sys.stderr)
        sys.exit(0)

    if METHOD == 'fix_range':
        assert CUTOFF is not None

    df = load_data(INPUT_FILE)
    if METHOD == 'curve':
        df_adj = filter_umi(df)
    elif METHOD == 'fix_range':
        df_adj = filter_umi_range(df, CUTOFF)
    
    # rename cell BC
    if NUMBER != 1:
        rename_BC(df_adj, NUMBER)
    
    df_adj.to_pickle(OUTPUT_FILE)

if __name__ == '__main__':
    main()
