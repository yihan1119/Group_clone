#!/usr/bin/env python3
'''
This version use the directional method to remove the UMIs.
'''
import sys
import argparse
import numpy as np
import editdistance
from collections import Counter, defaultdict

#get file handles based on the suffix
def get_handle(input_file):
    '''Get the file handle for the input file'''
    if input_file.endswith('gz'):
        import gzip
        handle = gzip.open(input_file, mode='rt')
    elif input_file.endswith('bz2'):
        import bz2
        handle = bz2.open(input_file, mode='rt')
    else:
        handle = open(input_file, 'r')
    return handle

#decide which cell barcode to choose for each read
def get_cell_bc(barcode_info):
    '''
    Examples: 
    ACGGAGACAACACCTA	0:16	1:0:0
    GTCTTCGAGATCGATA,GTCTTCGAGTCGATAA,TCTTCGGAGTCGATAA	0:15,0:16,1:16	0:0:1,0:0:0,0:0:1

    Rules for choosing Cell Barcodes: 
    (1) Check if match is unique. 
    (2) If not unique, check if there is just one match at the position 0:16.
    (3) If none of the criteria matches, discard the read. 
    '''
    seq_list   = barcode_info[0].split(',')
    span_list  = barcode_info[1].split(',')
    
    count = Counter(span_list)
    
    if (len(seq_list) == 1):
        return [seq_list[0], span_list[0]]
    elif (count["0:16"] == 1):
        for i in zip(seq_list, span_list):
            if (i[1] == '0:16'):
                return [i[0], i[1]]
    else:
        return ['NA', 'NA']

#decide which sgRNA to choose for each read    
def get_sgRNA(sgRNA_info):
    '''
    Example:
    GTAAACTATCGTGGCGTGG,TAAACTATCGTGGCGTGGT	18:37,19:38	0:0:0,0:0:0

    Rules for choosing sgRNAs:
    Some enhancers have sgRNAs that are staggered (off by one or two base pairs). Therefore, multiple sgRNAs
    could be returned, but the correct one should have to perfect match at the right postion(19:38). 
    
    (1) Check if match is unique. 
    (2) If not unique, check if there is just one match at the position 19:38.
    (3) If not unique match at postion 19:38, check if there is any perfect match.
    (4) If none of the criteria matches, discard the read. 
    '''
    seq_list   = sgRNA_info[0].split(',')
    span_list  = sgRNA_info[1].split(',')
    match_list = sgRNA_info[2].split(',')
    
    span_count = Counter(span_list)
    match_count = Counter(match_list)
    if (len(seq_list) == 1):
        return seq_list[0]
    elif (span_count["19:38"] == 1):
        for i in zip(seq_list, span_list):
            if (i[1] == "19:38"):
                return i[0]
    elif (match_count["0:0:0"] == 1):
        for i in zip(seq_list, span_list, match_list):
            if (i[2] == "0:0:0"):
                return i[0]
    else:
        return 'NA'
    
#Create an N-dimentional nested dictionary
def nested_dict(n, type):

    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))
    
#Summarize cell barcode and sgRNA into dictionary
def summarize_reads(input_file, umi_length = 10):

    data_dict = nested_dict(2, list)
    
    counter = 0
    print_counter = 0
    non_print_counter = 0

    #summarize the reads into a 2D dictionary as:
    # data_dict->cell_barcode->sgRNA_sequence = umi_list
    with get_handle(input_file) as f:
        for line in f:
            counter += 1
            i = line.strip().split('\t')
            barcode_info   = get_cell_bc(i[1:3])
            sgRNA_sequence = get_sgRNA(i[-3:])
            
            if ((barcode_info[0] == 'NA') or (sgRNA_sequence == 'NA')):
                non_print_counter += 1
            else:
                print_counter += 1
                read_sequence = i[0]
                cell_barcode  = barcode_info[0]
                barcode_span  = barcode_info[1].split(':')
                umi_start     = int(barcode_span[1])
                umi_end       = umi_start + umi_length
                umi           = read_sequence[umi_start : umi_end]

                data_dict[cell_barcode][sgRNA_sequence].append(umi)
    
    return data_dict

'''
Use Directional Adjacency method to define the UMIs of each sgRNA within every cell.
Modified from the algorithm used in UMI-tools
'''
#Function to calculate edit distance
def calc_edit_dist(first, second):
    return editdistance.eval(first, second)

#Define the UMI based on the Directional Adjacency (2 times less reads that the father nodes)
def get_adj_list_directional_adjacency(barcodes, barcode_counts, threshold=1):
    """barcode_counts: Counter object"""
    return {barcode: [barcode2 for barcode2 in barcodes
                  if 0 < calc_edit_dist(barcode, barcode2) <= threshold and
                  barcode_counts[barcode] >= (barcode_counts[barcode2] * 2) - 1]
                  for barcode in barcodes}

#Beadth First Search Tool
def breadth_first_search(node, graph):
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))

    while queue:
        node = (list(queue))[0]
        found.update(graph[node])
        queue.update(graph[node])
        searched.update((node,))
        queue.difference_update(searched)

    return found

#Get all the adjancecy nodes
def get_connected_components_adjacency(graph, barcode_counts):
    """barcode_counts: Counter object"""
    found = []
    components = []

    for node in sorted(graph, key=lambda x: barcode_counts[x], reverse=True):
        if node not in found:
            component = breadth_first_search(node, graph)
            found.extend(component)
            components.append(component)

    return components

#Dedup UMI main function
def dedup_umis(umi_list, threshold=1):
    umis = Counter(umi_list)
    dir_adj = get_adj_list_directional_adjacency(umis.keys(), umis, threshold)
    clusters = get_connected_components_adjacency(dir_adj, umis)

    return len(clusters)

#Count UMI main function
def count_umi(data_dict, output_fh, edit_distance_threshold):
    '''
    Count and output the sgRNA read count.
    Use directional UMI method.
    '''
    for i in data_dict.keys():
        #print cell barcode
        print(i, end ='\t', file=output_fh)

        output_sgRNA_list = np.array([])
        output_read_list  = np.array([])
        output_umi_list   = np.array([])
        for ii in data_dict[i].keys():
            umi_list = data_dict[i][ii]
                
            read_count = len(umi_list)
            umi_count  = dedup_umis(umi_list, edit_distance_threshold)
#            umi_count  = len(set(list))            
            
            output_sgRNA_list = np.append(output_sgRNA_list, ii)
            output_read_list  = np.append(output_read_list, read_count)
            output_umi_list   = np.append(output_umi_list, umi_count)
            
        #sort the output sgRNA based on the order of their umi numbers
        output_index = np.argsort(-output_umi_list)
        
        #print 
        print(np.sum(output_read_list).astype(int),end='\t', file=output_fh)
        print(len(output_read_list),end='\t',file=output_fh)
        print(';'.join(output_sgRNA_list[[output_index]]),end ='\t', file=output_fh)
        print(';'.join(map(str, output_read_list[[output_index]].astype(int))),end ='\t', file=output_fh)
        print(';'.join(map(str, output_umi_list[[output_index]].astype(int))), file=output_fh)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', required=True,
                        type=str,
                        help='specify an input file \
                        (output of evaluate_barcodes_from_reads.py)')
    parser.add_argument('-o', '--output', dest='output', required=False,
                        help='specify an output file')
    parser.add_argument('-m', '--edit_dist', required=False,
                        type=int,
                        default=1,
                        help='specify a edit distance threshold \
                        for UMI deduplication. \
                        The default is 1')


    args = parser.parse_args()
    input_file  = args.input
    output_file = args.output
    edit_distance_threshold = args.edit_dist
    
    output_fh = open(output_file, 'w')
    data_dict = summarize_reads(input_file)
    umi_adj_data_dict = count_umi(data_dict, output_fh, edit_distance_threshold)
    output_fh.close()

if __name__ == '__main__':
    main()
    
