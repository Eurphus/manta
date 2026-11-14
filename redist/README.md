# External packages

This directory contains external packages for building manta.

## Package modification notes

### boost

To reduce size, certain files have been removed from the
full boost source distribution acccording to:

    ${ROOT_PATH}/scratch/make_boost_subset.bash

### htslib/samtools

htslib and samtools were updated from 1.9 to 1.25 in this fork. htslib 1.9 aborts at
startup on current zlib releases (`bgzf_hopen: Assertion 'compressBound(BGZF_BLOCK_SIZE)
< BGZF_MAX_BLOCK_SIZE' failed`), cannot read CRAM 3.1, and predates a series of CRAM
decoder security fixes.

To reduce size, both packages have been modified from the upstream release tarballs
(<https://github.com/samtools/htslib/releases>, <https://github.com/samtools/samtools/releases>)
by removing files which are not needed to build the static library and the tools used by the
workflow. No source files were edited:

- htslib: removed `test/`, `htscodecs/tests/` and `samples/`
- samtools: removed `test/`, `examples/` and the bundled `htslib-1.25/` copy

Only the `samtools` program itself is built (`make samtools`) because the default `all`
target also builds the removed test programs.

### pyflow

pyflow 1.1.20 is included as plain source in `pyflow/` instead of the original release
tarball, because it carries the python3 port required by this fork. The demo and doc
directories of the pyflow release are not included.

### cmake-modules

cmake-modules-c99fd3 modified to show git describe --dirty

