# Cyp6 cluster UGS

Haplotype analyses of the CYP6 cluster in *Anopheles gambiae* samples from Uganda (UGS), using Ag1000G Phase 1 data.

## Output

Outputs are in these folders:

* `haplotype_analysis_25set19_CYP6P4_I236M.py`: python script to carry out analyses
* `haplotype_analysis_output_CYP6P4_I236M`: all figures and tables related to the haplotype cluster and positive selection analyses.
* `haplotype_phylogeny`: phylogenies of haplotypes using variants in the cluster (using IQTREE). Results in newick (`.treefile`) & pdf.

Linkage disequilibrium analyses are carried out across the entire CYP6 cluster, but haplotype clustering and downstream analyses are specific to the CYP6P4 gene (`I236M` allele +/- 2000bp ).

Default cluster names are not very informative -- here's what they actually are:

* `cluster_0`: a cluster of UGS haplotypes carrying duplications, ZZB, P4 mutations etc.
* `cluster_1`: a wt cluster of UGS haplotypes, that doesn't have any interesting mutation
* `cluster_no_wt`: all other wt haplotypes with less than 1% frequency in UGS are put in this group
* `cluster_no_alt`: all other non-wt haplotypes  with less than 1% frequency in UGS are put in this group

All other folders and files are used for input/metadata/etc.

## Input

Download Ag1000G Phase 1 data from here

```
ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3
```
