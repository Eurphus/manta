Manta Workflow Demo
-------------------

This directory contains a small dateset which can be used to verify
correct installation and demonstrate basic elements of the
workflow. To run the demonstration, run the demo script found in the
installation bin directory:

```
python ${MANTA_INSTALL_PATH}/bin/runMantaWorkflowDemo.py
```

This script creates a `MantaDemoAnalysis` directory under the current
working directory, runs Manta on a small demo dataset, and compares
the candidate and germline structural variant output to an expected
result.

The demo data contain reads from the cell lines HCC1954 and HCC1954BL
mapped in the vicinity of translocation breakends corresponding to
COSMIC variant [COST16011][1]. The demo sequencing data is extracted
from TCGA Benchmark 4. The two alignment files are analyzed as two
samples of a joint germline analysis.

The expected result has three candidate junctions (six breakend
records in `candidateSV.vcf.gz`) and no records in `diploidSV.vcf.gz`.
The junctions are somatic rearrangements in a highly amplified region
of the HCC1954 tumor cell line. They are found as candidates, but are
supported by only 1% to 4% of the reads at their breakends, so the
germline model genotypes both samples as homozygous reference. The
demo therefore checks that the workflow runs and that candidate
discovery on real data is reproduced, it does not demonstrate a
germline call. See `scratch/simulation` in the source tree for a
benchmark of germline calling.

[1]:http://grch37-cancer.sanger.ac.uk/cosmic/rearrangement/overview?id=16011
