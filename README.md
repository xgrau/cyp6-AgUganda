# Cyp6 cluster UGS

Haplotype analyses of the CYP6 cluster in *Anopheles gambiae* samples from Uganda (UGS), using Ag1000G Phase 1 data.

## Input

Download Ag1000G Phase 1 data from here

```bash
ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3
```

## Contents

Notebooks:

* `2020-02-27_hap_analysis_UGS_sweep.ipynb`: notebook with analysis of haplotype similarity and selection signals in UGS. The cluster of swept haplotypes is defined using genetic distances of genotypes located +/- 1000 variants around the *Cypaa1* start coordinate (results taken from `sweep_haplotypes.csv`). Selection statistics are estimated for the Cyp6 custer (EHH, Garud H, hap diversity). Results and plots go to folder `results_sweep`.

* `results_phylogeny`: phylogenies of haplotypes using variants in the cluster (using IQTREE). Results in newick (`.treefile`) & pdf.

* in folder `other/`, there is another notebook I used to do a quick analysis of haplotype networks built around another variant of interest, `2R:28491424`, which was in LD with both the duplication and the ZZB insertion in UGS. I wanted to see if we can find it elsewhere in the Ag1k samples, but it's not there, so this is not useful. Results go to `results_TAG`.

## Results

### Main cluster

Phylogeny & haplotype cluster shows a group of swept haplotypes with strong signals of selection. It contains the following genotypes:

* all sequences have the *I236M* mutation
* all sequences have the *zzb* transposon OR come from genomes that have it (het)
* all sequences have the *indel* OR come from genomes that have it (het)
* it contains sequences from specimens with and without *duplication* (including duplicated homozygotes).
* all of them have the tagging variant that is in LD with the *duplication* and *zzb* (`tag1` in sequence name)

Therefore, the selective sweep occurred **after** the *I236M*, *zzb* insertion and indels; but **before** the duplication occurred. This is consistent with the late emergence of the duplication (2009) relative to other mutations (*zzb* in 2004, *I236M* in 2005).

### Secondary cluster

There is a second cluster of highly similar haplotypes that do not have any resistance mutations. This cluster has signals of selection too. Not explored right now.
