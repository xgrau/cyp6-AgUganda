# Cyp6 cluster UGS

Haplotype analyses of the CYP6 cluster in *Anopheles gambiae* samples from Uganda (UGS), using Ag1000G Phase 1 data.

## Input

Download Ag1000G Phase 1 data from here

```bash
ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3
```

## Contents

Notebooks:

* `2020-02-27_hap_analysis_CYP6P4_236M.ipynb`: notebook with analysis of haplotype similarity and selection signals in UGS. Haplotype networks are built from variants around CYP6P4 236M allele, and selection statistics are estimated for the Cyp6 custer (EHH, Garud H). Results and plots go to folder `results_236M_UGS`.

* `results_phylogeny`: phylogenies of haplotypes using variants in the cluster (using IQTREE). Results in newick (`.treefile`) & pdf.

* `2020-02-27_hap_analysis_ZZBD_tag.ipynb`: quick analysis of haplotype networks built around another variant of interest, `2R:28491424`, which is in LD with both the duplication and the ZZB insertion in UGS. I wanted to see if we can find it elsewhere in the Ag1k samples, but it's not there, so this is not useful. Results go to `results_TAG`.

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

There is a second cluster of highly similar haplotypes that do not have the *I236M* mutation (`cluster_1`). This cluster has signals of selection, and the following genotypes:

* no *I236M* mutation
* mostly from specimens that are heterozygous for the *indel* and *zzb*, with some coming from specimens that don't have indels or zzb.
* mostly from specimens *without duplication*, with some coming from het specimens.
* yet, none of them has the tagging variant (`tag0`), so we can conclude that they're 100% wt.
