#!/usr/bin/python3
#jweep_detection_phase1.py

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
from re import *
import allel
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict as dict 
from re import *
from scipy import cluster

# This function is taken from Nick's selection package
def find_clusters(dist, n, threshold=0.001, method='complete'):
        # build hierarchy
        clust = cluster.hierarchy.linkage(dist, method=method)
        # find clusters
        f = cluster.hierarchy.fcluster(clust, threshold, criterion='distance')
        # compute cluster sizes
        fsz = np.bincount(f)
        # sort largest first
        fsort = np.argsort(fsz)[::-1]
        # take largest n
        fsort = fsort[:n]
        # get haplotype indices for each cluster
        clusters = [set(np.nonzero(f == i)[0]) for i in fsort]
        return clusters



sample_info_fn = "/home/xavi/dades/Variation/phase1.AR3.1/samples/samples.meta.txt"
haplotypes_fn = "/home/xavi/dades/Variation/phase1.AR3.1/haplotypes/ag1000g.phase1.ar3.1.haplotypes.2R.h5" 
accessibility_fn = "/home/xavi/dades/Variation/phase1.AR3.1/accessibility/accessibility.h5"
chrom = '2R'
pop = 'UGS'

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')
stdout.flush()

sample_info = pd.read_csv(sample_info_fn, sep = '\t', index_col = 'ox_code')
# Get the population information and repeat each entry twice to have a record of haplotype populations
haplotype_pop = np.repeat(sample_info['population'], 2)
# Get the names of the UGgam samples
focal_haplotypes = np.array(haplotype_pop.index[np.where(haplotype_pop == pop)])
focal_samples = focal_haplotypes[1::2]

# Load the haplotypes file
print('Loading haplotypes.')
stdout.flush()
all_haplotypes = h5py.File(haplotypes_fn, 'r')
pos = allel.SortedIndex(all_haplotypes[chrom]['variants']['POS'])
# Get the indices of the samples of interest
sample_list = pd.Series([x.decode('utf-8') for x in all_haplotypes[chrom]['samples'].value])
genotype_indices_list = np.where(sample_list.isin(focal_samples))[0]

print('Loading P4 calls.')
stdout.flush()
# Get calls for P4 from the haplotype file. 
P4_SNP_pos = 28497967
P4_calls = all_haplotypes[chrom]['calldata']['genotype'][pos.locate_key(P4_SNP_pos)]
# Although the genotype calls include the crosses, and the sample_info object doesn't, we can still used the
# UGS indices to pull out the right samples because the crosses are all at the end of the table
P4_calls_UGgam = P4_calls[np.where(sample_info.population == pop)[0], :]
# Get the genotype
P4_UGgam = np.sum(P4_calls_UGgam, 1)

# Load the accessibility
print('Loading accessibility.')
accessibility = h5py.File(accessibility_fn, "r")

# Identify the swept clusters

print('Finding swept cluster.')
## For CYP6AA1 (this is the start position):
focus = 28480576

# Define the range (either side of the focus) over which we look for SNPs
SNP_range = 500

# Identify segregating, non-singleton positions
hap_matrix = allel.GenotypeArray(all_haplotypes[chrom]['calldata']['genotype'][:,genotype_indices_list,:]).to_haplotypes()
ac = hap_matrix.count_alleles()
non_singleton = ac.min(1)>1
non_singleton_pos = pos[non_singleton]

# Get the range of SNPs that are within the region of interest, and that are segregating non-singletons
pre_SNPs = non_singleton_pos[non_singleton_pos <= focus][(-SNP_range):]
post_SNPs = non_singleton_pos[non_singleton_pos > focus][:SNP_range]
SNP_range = np.concatenate([pre_SNPs, post_SNPs])
SNP_pos = np.where([x in SNP_range for x in pos])[0]

focal_hap_matrix = hap_matrix[SNP_pos, :]
dist = allel.pairwise_distance(focal_hap_matrix, metric = 'hamming')
is_accessible = accessibility[chrom]['is_accessible'][SNP_range[0]:SNP_range[-1]]
n_bases = np.count_nonzero(is_accessible)
dist_dxy = dist * focal_hap_matrix.n_variants / n_bases
focal_clusters = find_clusters(dist_dxy, n=3, threshold=0.001)
largest_focal_cluster = focal_clusters[0]
focal_cluster_members = focal_haplotypes[list(largest_focal_cluster)]

# Record the focal cluster as binary, with 1 showing haplotypes in the swept cluster, and 0 showing haplotypes 
# that are not
focal_cluster_calls = np.zeros(focal_hap_matrix.shape[1])
focal_cluster_calls[list(largest_focal_cluster)] = 1

# Get genotype-level calls for the swept cluster and the P4 SNP
focal_cluster_genotype = [int(sum(focal_cluster_calls[i:i+2])) for i in range(0, len(focal_cluster_calls), 2)]

# Load the Dup1 and ZZB genotype calls, making sure the sample order is the same as in the haplotype data
genotype_calls = pd.read_csv('/home/xavi/Documents/cyp6-AgUganda/data/p1_UGS_extramutations_CYP6.csv', sep = '\t', index_col = 0).loc[focal_samples, ]

# Create the output table joining these data together. For Dup1, the sample with coverage 3 is likely to be 
# homozygote. 
print('Writing output table to file.')
contingency_allsamples = pd.DataFrame({'Sweep':focal_cluster_genotype,
                                       'P4': P4_UGgam, 
                                       'Dup1': genotype_calls.dupaa1.replace(3,2), 
                                       'ZZB': genotype_calls.zzb})
contingency_allsamples.to_csv('zz_sweep_contingency_UGgam_phase1.csv', sep = '\t')

# And output the swept clusters
cluster_calls_table = pd.DataFrame({'Sweep_haplotype': [int(x) for x in focal_cluster_calls]})
cluster_calls_table.to_csv('zz_sweep_haplotypes.csv', sep = '\t', index = False)









