# Changes in this fork

This is a modified version of Manta. It starts from commit `75b5c38` of
[Illumina/manta](https://github.com/Illumina/manta) (17 July 2019, three commits after the v1.6.0 release),
which is the last revision published under the GPLv3. Nothing from later upstream revisions is included. The
fork is distributed under the same license, see [LICENSE.txt](../LICENSE.txt).

**This branch calls germline variants only.** Somatic (tumor/normal), tumor-only and RNA-Seq calling have been
removed, see [germline only](#germline-only). The branch `gpl-main` has the same method changes with those
modes still present.

Every change to the calling methods has a switch, and the default selects the changed method. With all
switches set to their v1.6.0 values, the candidate and diploid VCF records and headers of the benchmark below
are identical to those of a build without any method change.

**Validation status: simulation only.** None of the accuracy figures here come from real sequencing data, and
the simulation does not include a read mapper. Read the [validation](#validation) section before relying on
any of them.

## Table of Contents

* [Summary](#summary)
* [Germline only](#germline-only)
* [Build and runtime](#build-and-runtime)
* [Method changes](#method-changes)
* [Failures fixed](#failures-fixed)
* [Reproducing v1.6.0 results](#reproducing-v160-results)
* [Validation](#validation)
* [Known gaps](#known-gaps)

## Summary

Change | Default | Switch (`configManta.py.ini`)
--- | --- | ---
More than one allele per small variant locus | 2 alleles | `maxLocalHaplotypeCount = 1`
Join left and right contigs of a large insertion | on | `enableInsertionContigMerge = 0`
Contig trimming which tolerates basecall errors | on | `enableContigConservativeRangeRefinement = 0`
Tandem duplications scored as duplications | on | `enableTandemDupModel = 0`
Uncertain zygosity no longer filters a variant | on | `filterUncertainZygosity = 1`
Depth ratio for deletions and duplications (`FORMAT/RDR`) | on | `reportDepthRatio = 0`
Near-identical duplicate records merged | 5% mismatches | `maxDuplicateAltMismatchFraction = 0`
Python 3, htslib/samtools 1.25, current CMake and compilers | | none
Basecall qualities above 70 accepted | | none
Reads which cannot be split-read aligned are skipped | | none
csi index accepted for `--callRegions` | | none
Somatic, tumor-only and RNA-Seq calling removed | | none (use branch `gpl-main`)

## Germline only

Every input sample is a germline sample, scored under the diploid model. What was removed:

* `configManta.py` options `--tumorBam`, `--rna` and `--unstrandedRNA`. Using one of them is an error.
  `--bam` and its synonym `--normalBam` are unchanged.
* Outputs `somaticSV.vcf.gz`, `tumorSV.vcf.gz` and `rnaSV.vcf.gz`, and with them the INFO keys `SOMATIC`,
  `SOMATICSCORE` and `JUNCTION_SOMATICSCORE` and the filter `MinSomaticScore`.
* Configuration file options `minSomaticScore`, `minPassSomaticScore`, `rnaMinCandidateVariantSize` and
  `enableRemoteReadRetrievalForInsertionsInCancerCallingModes`. A configuration file which still sets them is
  accepted and the values are ignored without a warning.
* The somatic and tumor-only scoring models, the RNA scoring model and the intron-aware contig aligner, the
  transcript strand bookkeeping, and the tests of all of these. About 5500 lines in total.

Germline results do not change. `diploidSV.vcf.gz`, `candidateSV.vcf.gz` and `candidateSmallIndels.vcf.gz`
were compared before and after the removal for 14 runs, and all 42 files have identical records and headers
(leaving out the header lines for date, command line, version and reference path). The runs are 7 inputs, each with the
default settings and with the [v1.6.0 settings](#reproducing-v160-results): the three simulated benchmark
samples, a joint analysis of two simulated samples (502 diploid records), and the real reads of the demo
analyzed as one sample with `--exome`, as two samples with `--exome`, and as two samples without.

Not covered by that comparison: CRAM input, `--callRegions`, more than two samples, and whole genome sized
input.

A tumor sample given to this branch as `--bam` is genotyped as a diploid germline sample. A variant in a
small fraction of its reads is then called `0/0` and is not reported in `diploidSV.vcf.gz`.

Some internals which only the removed modes needed are still present, because taking them out would change
shared code paths for no change in output: storage for reads below the minimum mapping quality in the
candidate data, the second (lower) mapping quality tier of evidence counting, and the read scanner option
`--ignore-anom-proper-pair`, which nothing sets any more.

## Build and runtime

* The workflow and all scripts run on Python 3 (3.6 or later). Python 2 is no longer supported. The bundled
  pyflow is a Python 3 port of pyflow 1.1.20.
* htslib and samtools are updated from 1.9 to 1.25. htslib 1.9 fails an internal assertion when built against
  a current zlib.
* CMake 3.5 or later is required. The build is tested with GCC 13 and CMake 3.28 on Ubuntu 24.04.
* The workflow version is found from `git describe`, also when building from a git worktree. Output from this
  fork reports a version such as `1.6.0-25-g<hash>` in VCF headers, which is how to tell it from Illumina's
  1.6.0.

## Method changes

### More than one allele per locus

v1.6.0 aligns every contig of a local assembly to the reference and then uses only the best one. A site where
the two haplotypes carry different variants yields one of them at most. Now up to two contigs are used: after
the best contig, the next best one adds any candidate which is not already present (breakends or insertion
length differ by more than 3 bases), if at least two reads support it.

Each allele is still scored against the reference on its own. There is no joint genotype, so a site with two
different variant alleles gives two records and never a `1/2` genotype. Where neither haplotype has the
reference allele, for instance two different insertions at one position, both records tend to be genotyped
`1/1`.

Records in this situation are marked with the `OVERLAPPING_ALT` flag in INFO: another allele from a different
contig of the same assembly overlaps this one. Treat the genotype of a flagged record as 'carries the
variant, zygosity unknown'. In the benchmark, 75% to 83% of flagged records have the right zygosity, against
91% to 93% of records with neither this flag nor a `MinGQ` sample filter.

A second contig can also repeat a variant which a neighboring locus reports. Records at one position whose
assembled alleles have the same length and differ at no more than 5% of bases are now merged into the best
one. Previously only records with identical alleles were merged.

### Large insertions

Two changes increase how often a large insertion is fully assembled:

* The part of a contig supported by at least two reads (the 'conservative range', to which a contig is
  trimmed) was found from exact word matches during assembly. One basecall error in a read breaks every word
  across it, so contigs were trimmed harder as error rates went up, sometimes into their reference anchor.
  The range is now recomputed by aligning the reads to the finished contig with a mismatch allowance.
* When no contig spans an insertion, the left and right contigs are now joined if the end of one overlaps the
  start of the other by at least 20 bases of inserted sequence. The join is refused if more than one overlap
  fits or the overlap is low complexity. Without a join, the insertion is reported as before, with
  `LEFT_SVINSSEQ` and `RIGHT_SVINSSEQ` and no length.

### Tandem duplications

A haplotype with a tandem duplication keeps both reference junctions next to the one new junction. v1.6.0
scores a duplication as if the reference junctions were gone, like any other variant. For a duplication
called as a single junction:

* **Zygosity.** A heterozygous duplication is supported by 1 in 5 informative reads and a homozygous one by 1
  in 3, not 1 in 2 and 99 in 100. v1.6.0 cannot genotype a duplication `1/1`. The choice between `0/1` and
  `1/1` now uses the fractions expected for a duplication. Whether the sample has the duplication at all is
  still decided with the standard fractions, so QUAL and the set of duplications called are not changed by
  this part.
* **Fragments across two junctions.** For a duplication smaller than a fragment, one read often crosses the
  duplication junction while its mate crosses a reference junction. v1.6.0 then compares the two reads against
  each other, which cancels the evidence or turns it into reference support. Both alleles are now scored on
  the read which crosses the duplication junction.
* **Paired-read evidence** is phased out for duplications smaller than the fragment size, as it always has
  been for small deletions and insertions.

Expect lower GQ for duplications than before. Telling 1 in 5 from 1 in 3 at 30x is harder than the old model
made it look, and the lower GQ is the honest one.

### Uncertain zygosity

GQ is low whenever two genotypes are similarly likely. If the two are `0/1` and `1/1`, the variant itself is
not in doubt, but v1.6.0 marks the sample `MinGQ` and then the record `SampleFT`, which removes it from the
passing calls of a single sample analysis. Now a sample which fails only `MinGQ`, has a non-reference
genotype, and has a phred-scaled reference genotype probability of at least 15 counts as passing when
`SampleFT` is decided. `FORMAT/FT` still says `MinGQ` and GQ is unchanged.

**This changes which records are `PASS`.** Records with `FILTER=PASS` and `FT=MinGQ` in every sample did not
exist before. A pipeline which needs a confident genotype should test `FORMAT/FT` as well as `FILTER`.

If `filterUncertainZygosity = 1` is set, also set `enableTandemDupModel = 0`. The duplication model lowers GQ
for duplications, and under the old filter rule many real duplications are then filtered.

### Depth ratio

For deletions and tandem duplications of at least 500 bases, the diploid VCF reports `FORMAT/RDR`: mean read
depth within the variant divided by mean depth of the 1000 bases on each side. Expect about 0.5 and 0 for a
heterozygous and homozygous deletion, 1.5 and 2 for duplications, and 1 when the junction is not a copy number
change. Variants over 20 kb are sampled in 20 windows of 1 kb.

RDR is reported only. It is not used in QUAL, GT or any filter. It is not corrected for GC content or
mappability, so it is least reliable for small variants and in repeats.

## Failures fixed

* A read with a basecall quality above 70 used in scoring ended the run with an exception. Read merging and
  consensus tools write such values. They are now treated as 70. An input with a quality offset error is no
  longer rejected as a side effect.
* A read longer than the allele sequence built for a candidate, or an alignment record with only a few
  sequenced bases (such as a mostly hard-clipped supplementary alignment), ended the run with an exception in
  the split read aligner. These reads are now skipped.
* `--callRegions` required a `.tbi` index. A `.csi` index is now accepted.
* The `enableEvidenceSignalFilter` configuration option had no effect. Setting it to 0 now turns the filter
  off. The default (on) is unchanged.
* The serialization headers of boost 1.74 and 1.75 fail to compile unless a definition is included first. The
  include is added. **Not tested**: only the bundled boost 1.58 was available.

## Reproducing v1.6.0 results

Copy `bin/configManta.py.ini`, set the following, and pass the copy with `--config`:

```
maxLocalHaplotypeCount = 1
enableInsertionContigMerge = 0
enableContigConservativeRangeRefinement = 0
enableTandemDupModel = 0
filterUncertainZygosity = 1
reportDepthRatio = 0
maxDuplicateAltMismatchFraction = 0
```

The same switches exist as `GenerateSVCandidates` options (`--help`).

## Validation

What has been checked:

* All unit tests pass, including tests added for each change.
* On branch `gpl-main`, the bundled demo (HCC1954 tumor/normal region) reproduces its expected somatic calls
  with the default settings. On this branch the demo is a different test, see
  [the demo README](../src/demo/README.md).
* Germline output is identical before and after the removal of the other calling modes, see
  [germline only](#germline-only).
* With the seven switches above, benchmark output is identical to a build of the unmodified methods.
* The branch builds from a fresh clone, and every commit which changes C++ code was built and unit tested
  on its own.
* A synthetic benchmark, in [scratch/simulation](../scratch/simulation). `runSimBenchmark.sh` reproduces it.

The benchmark simulates a diploid sample at 30x over 3 Mb of GRCh38 chr21, with variants of 50 bases to 8 kb,
chimeric pairs, reads with junk tails, and a set of 'artifact' junctions at 10% of sample depth which should
not be called.

Mixed variant set used during development (214 alleles) | v1.6.0 methods | this fork
--- | --- | ---
Recall | 0.776 | 0.967
Precision | 0.971 | 0.976
Genotype concordance | 0.867 | 0.889
Compound heterozygous sites, recall | 27/52 | 50/52
Large insertions (500 to 3000 bases), recall | 4/10 | 9/10

Held-out mixed variant set (216 alleles) | v1.6.0 methods | this fork
--- | --- | ---
Recall | 0.843 | 0.972
Precision | 0.948 | 0.955
Genotype concordance | 0.885 | 0.862
Compound heterozygous sites, recall | 33/54 | 50/54
Large insertions (500 to 3000 bases), recall | 6/8 | 8/8

The held-out set uses a different region and random seed, and was first run after the method changes were
finished. No method was changed in response to it. It did expose the duplicate record problem described
above, which was fixed in the VCF merge step.

Genotype concordance over all passing calls does not improve, and falls on the held-out set. The calls which
are gained are mostly at compound heterozygous sites and at variants with uncertain zygosity, and those are
the calls whose genotypes are least reliable. Both groups are marked (`OVERLAPPING_ALT`, `FT=MinGQ`). For
records with neither mark, concordance goes from 0.88 to 0.91 on the held-out set and from 0.86 to 0.93 on
the development set.

Tandem duplications (246 alleles) | v1.6.0 methods | this fork
--- | --- | ---
Recall | 0.919 | 0.972
Precision | 0.970 | 0.972
Genotype concordance | 0.611 | 0.908

Both versions report the same false positives on all three sets. All of them are at simulated artifact
junctions, apart from one 101 base duplication in the duplication set which is reported as an insertion of
the same length and so is not matched by the evaluation.

Why these figures are an upper bound and not an estimate:

* **There is no read mapper in the simulation.** Each read is placed where a mapper following bwa-mem
  conventions would ideally place it. Mismapping, which causes most false SV calls in real data, is absent.
* The simulated region is unique sequence on one chromosome arm. Real duplications and insertions concentrate
  in repeats and segmental duplications.
* The 10% artifact haplotype is an arbitrary stand-in for low level noise. It also puts reference reads at
  every homozygous variant, which inflates the benefit of the zygosity filter change.
* Accuracy was measured for one sample at a time. Joint calling of several samples was run (two simulated
  samples, and the two demo samples) and its output is stable across the removal, but it was not compared to
  a truth set.

Before using this fork for anything that matters, run it and v1.6.0 on a sample with a truth set (for
instance HG002 with the GIAB SV benchmark), per variant type and size, and compare with each switch.

## Known gaps

Planned and not done:

* Mobile element annotation of insertions.
* Mapping quality aware read likelihoods. The scoring model still uses constant mismapping probabilities.
  A related weakness: a read which matches the variant allele with two or more mismatches contributes almost
  no evidence, because its likelihood falls below a fixed noise floor.
* A joint genotype for two variant alleles at one locus (see above).
* Using read depth in the genotype of duplications and deletions. RDR separates heterozygous from homozygous
  duplications well in simulation, and is not used.

Open upstream issues which were looked at and not resolved: failure to estimate fragment statistics on
DRAGMAP alignments (#299), input restricted to FR read pairs (#309), no single-record representation of
inversions (#303).
