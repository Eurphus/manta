Manta Workflow Demo
-------------------

This directory contains two small datasets which can be used to verify
correct installation and demonstrate basic elements of the
workflow. To run the demonstration, run the demo script found in the
installation bin directory:

```
python3 ${MANTA_INSTALL_PATH}/bin/runMantaWorkflowDemo.py
```

This script creates a `MantaDemoAnalysis` directory under the current
working directory, runs Manta on each demo dataset in a subdirectory
of its own, and compares the candidate and germline structural variant
output of each analysis to an expected result.

The demo verifies that an installation reproduces the output of the
build it was made from. It is not a measurement of accuracy.

### SimulatedGermline

A simulated diploid sample at 30x depth over 80 kb of chromosome 21
(GRCh38 chr21:33200001-33280000; positions in the demo files count from
the start of this region). The reads are not real sequencing data and
were not mapped by a read mapper: each read is written with the
alignment an ideal mapper would give it. The sample was created by
`scratch/simulation/simulateSVReads.py` in the source tree, with the
options recorded in the header of the alignment file.

The simulated variants are listed in
`data/SimulatedGermline.chr21_region.truth.tsv`. The expected result
in `diploidSV.vcf.gz` is:

Simulated variant | Simulated genotype | Expected call
--- | --- | ---
300 base insertion | heterozygous | `0/1`, PASS
1404 base deletion | homozygous | `1/1`, PASS, depth ratio (`RDR`) 0.06
1381 base inversion | heterozygous | 4 breakend records, `0/1`, PASS
215 base deletion | heterozygous | `0/1`, PASS
5720 base tandem duplication | heterozygous | `1/1` with GQ 3 and sample filter `MinGQ`, record PASS, `RDR` 1.42
44 base insertion | heterozygous | none, it is below the minimum scored size and is in `candidateSmallIndels.vcf.gz`
1950 base deletion in 10% of reads | not a variant of the sample | none

The genotype of the duplication is wrong. The call shows how to read
such a record: the variant passes, the sample filter and the genotype
quality say that its zygosity is not known, and the depth ratio of 1.42
is the value expected for a heterozygous duplication (1.5). The depth
ratio is reported only and is not used for the genotype.

### HCC1954

Reads from the cell lines HCC1954 and HCC1954BL mapped in the vicinity
of translocation breakends corresponding to COSMIC variant
[COST16011][1]. The sequencing data is extracted from TCGA Benchmark 4.
This is real sequencing data, of a tumor cell line and its matched
normal. The two alignment files are analyzed as two samples of a joint
germline analysis.

The expected result has three candidate junctions (six breakend
records in `candidateSV.vcf.gz`) and no records in `diploidSV.vcf.gz`.
The junctions are somatic rearrangements in a highly amplified region
of the tumor cell line. They are found as candidates, but are
supported by only 1% to 4% of the reads at their breakends, so the
germline model genotypes both samples as homozygous reference. This
analysis checks that candidate discovery on real data is reproduced
and that these junctions are not called as germline variants.

[1]:http://grch37-cancer.sanger.ac.uk/cosmic/rearrangement/overview?id=16011
