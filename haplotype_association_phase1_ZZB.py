#!/usr/bin/python3
#haplotype_association.py

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



phased_vcf_fn = 'mvncall_phase1/mvncall_200snps_lambda01_UGgam_fixed.vcf'
sample_info_fn = '/home/eric/Liverpool/AR3/samples/samples.meta.txt'
haplotypes_fn = '/media/eric/Ucalegon/haplotypes_phase1/hdf5/ag1000g.phase1.ar3.1.haplotypes.2R.h5'
accessibility_fn = "/media/eric/Ucalegon/accessibility/accessibility.h5"
chrom = '2R'
pop = 'UGS'

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')
stdout.flush()

# Load the haplotypes file
print('Loading haplotypes.')
stdout.flush()
all_haplotypes = h5py.File(haplotypes_fn, 'r')
pos = allel.SortedIndex(all_haplotypes[chrom]['variants']['POS'])

print('Loading phased vcf file.')
stdout.flush()
# Read the vcf file for the newly phased variants
callset = allel.read_vcf(phased_vcf_fn)
# Get the metadata
sample_info = pd.read_csv(sample_info_fn, sep = '\t', index_col = 'ox_code')
# Get the population information and repeat each entry twice to have a record of haplotype populations
haplotype_pop = np.repeat(sample_info['population'], 2)
# Get the haplpotype calls for the mutations of interest. These are initially produced as a 3-dimensional 
# array and are then flattened to 2 dimensions. 
mut_calls_3d = callset['calldata/GT']
mut_calls = mut_calls_3d.reshape(mut_calls_3d.shape[0], mut_calls_3d.shape[1]*2)

# Get calls for Dup1 and ZZB
mut_names = callset['variants/ID']
Dup1_UGgam = mut_calls[mut_names == 'Dup1'].flatten()
ZZB_UGgam = mut_calls[mut_names == 'ZZB'].flatten()

# Get calls for P4 from the haplotype file. 
P4_genotype = all_haplotypes[chrom]['calldata']['genotype'][pos.locate_key(28497967)]
# Although the genotype calls include the crosses, and the sample_info object doesn't, we can still used the
# UGS indices to pull out the right samples because the crosses are all at the end of the table
P4_genotype_UGgam = P4_genotype[np.where(sample_info.population == pop)[0], :]
# Get the haplotypes
P4_UGgam = P4_genotype_UGgam.flatten()

# Load the accessibility
accessibility = h5py.File(accessibility_fn, "r")

# Identify the swept clusters

## For CYP6AA1 (this is the start position):
focus = 28480576

# Find the index of the last SNP at or just before the focal point (I don't think I can use scikit-allele's 
# locate functions here because I don't know the exact position of the SNPs beforehand. 
focus_preloc = max(np.where(pos <= focus)[0])

# Get the names of the UGgam samples
focal_haplotypes = np.array(haplotype_pop.index[np.where(haplotype_pop == pop)])
focal_samples = np.unique(focal_haplotypes)

# Get the indices of the samples of interest
sample_list = pd.Series([x.decode('utf-8') for x in all_haplotypes[chrom]['samples'].value])
genotype_indices_list = np.where(sample_list.isin(focal_samples))[0]

# Define the starting window size (in SNPs)
starting_win_size = 1000

focal_start_index = focus_preloc - int(starting_win_size/2)
focal_end_index = focus_preloc + int(starting_win_size/2)
focal_hap_matrix = allel.GenotypeArray(all_haplotypes[chrom]['calldata']['genotype'][range(focal_start_index, focal_end_index)][:,genotype_indices_list,:]).to_haplotypes()
dist = allel.stats.pairwise_distance(focal_hap_matrix, metric = 'hamming')
is_accessible = accessibility[chrom]['is_accessible'][pos[focal_start_index]:pos[focal_end_index]]
n_bases = np.count_nonzero(is_accessible)
dist_dxy = dist * focal_hap_matrix.n_variants / n_bases
focal_clusters = find_clusters(dist_dxy, n=3, threshold=0.001)
largest_focal_cluster = focal_clusters[0]
focal_cluster_members = focal_haplotypes[list(largest_focal_cluster)]
focal_cluster2_members = focal_haplotypes[list(focal_clusters[1])]
focal_cluster3_members = focal_haplotypes[list(focal_clusters[2])]

# Record the focal cluster as binary, with 1 showing haplotypes in the swept cluster, and 0 showing haplotypes that are not
focal_cluster_genotype = np.zeros(focal_hap_matrix.shape[1])
focal_cluster_genotype[list(largest_focal_cluster)] = 1

# Let's associate the sweep with the mutations
contingency_allsamples = pd.DataFrame({'Sweep':focal_cluster_genotype, 'P4': P4_UGgam, 'Dup1': Dup1_UGgam, 'ZZB': ZZB_UGgam})
contingency_allsamples.to_csv('sweep_contingency_UGgam_phase1.csv', sep = '\t')










