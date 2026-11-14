Manta Structural Variant Caller
===============================

**This is a modified version of Manta.** It is based on the last revision of
[Illumina/manta](https://github.com/Illumina/manta) published under the GPLv3
(commit 75b5c38, July 2019, three commits after release v1.6.0), with changes
made in November 2026. It is not an Illumina release and is not supported by
Illumina. The changes, their defaults, how to turn each of them off and what
has and has not been validated are described in
[docs/forkChanges.md](docs/forkChanges.md). The accuracy of the changed
methods has been measured on simulated data only.

The remainder of this file describes Manta as released by Illumina, updated
where this version differs.

Manta calls structural variants (SVs) and indels from mapped
paired-end sequencing reads. It is optimized for analysis of germline
variation in small sets of individuals and somatic variation in
tumor/normal sample pairs. Manta discovers, assembles and scores
large-scale SVs, medium-sized indels and large insertions within a
single efficient workflow. The method is designed for rapid analysis
on standard compute hardware: NA12878 at 50x genomic coverage is
analyzed in less than 20 minutes on a 20 core server, and most WGS
tumor/normal analyses can be completed within 2 hours. Manta combines
paired and split-read evidence during SV discovery and scoring to
improve accuracy, but does not require split-reads or successful
breakpoint assemblies to report a variant in cases where there is
strong evidence otherwise. It provides scoring models for germline
variants in small sets of diploid samples and somatic variants in
matched tumor/normal sample pairs. There is experimental support for
analysis of unmatched tumor samples as well. Manta accepts input read
mappings from BAM or CRAM files and reports all SV and indel inferences
in VCF 4.1 format. See the [user guide][UserGuide] for a full
description of capabilities and limitations.

[UserGuide]:docs/userGuide/README.md

Methods and benchmarking details are described in:

Chen, X. *et al.* (2016) Manta: rapid detection of structural variants and
indels for germline and cancer sequencing applications. *Bioinformatics*,
32, 1220-1222. [doi:10.1093/bioinformatics/btv710][bpaper]

...and the corresponding [open-access pre-print][preprint].

[bpaper]:https://doi.org/10.1093/bioinformatics/btv710
[preprint]:https://doi.org/10.1101/024232


License
-------

Manta source code is provided under the [GPLv3 license](LICENSE.txt).
Manta includes several third party packages provided under other
open source licenses, please see [COPYRIGHT.txt](COPYRIGHT.txt)
for additional details.


Getting Started
---------------

This version must be built from source, please see the
[installation instructions](docs/userGuide/installation.md) for full
build and installation details. After installation the build can be
tested by [running a small demo](docs/userGuide/installation.md#demo).
The binary distributions on the [Illumina Manta releases page][releases]
are builds of the unmodified v1.6.0 and earlier releases, they do not
contain any of the changes made here.

[releases]:https://github.com/Illumina/manta/releases


Data Analysis and Interpretation
--------------------------------

After completing installation, see the [Manta user guide][UserGuide]
for instructions on how to run Manta, interpret results and estimate
hardware requirements/compute cost, in addition to a high-level methods
overview.


Manta Code Development
----------------------

For manta code development and debugging details, see the
[Manta developer guide][DeveloperGuide]. This includes details
on Manta's developement protocols, special build instructions,
recommended workflows for investigating
calls, and internal documentation details.

[DeveloperGuide]:docs/developerGuide/README.md
