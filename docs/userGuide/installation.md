Manta User Guide - Installation
===============================

[User Guide Home](README.md)

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Prerequisites to build from source](#prerequisites-to-build-from-source)
* [Runtime prerequisites](#runtime-prerequisites)
* [Operating System Guidelines](#operating-system-guidelines)
    * [Linux](#linux)
    * [OS X](#os-x)
    * [Windows](#windows)
* [Linux Package Additions](#linux-package-additions)
    * [Ubuntu 14.04 and 16.04](#ubuntu-1404-and-1604)
    * [Ubuntu 12.04](#ubuntu-1204)
    * [CentOS 7](#centos-7)
    * [CentOS 6](#centos-6)
* [Build procedure](#build-procedure)
    * [Workflow relocation](#workflow-relocation)
* [Demo](#demo)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


It is recommended to start from one of the [binary distributions on
the releases page][releases] if a suitable binary is available
(note that the CentOS 6 binary is expected to support a
large variety of linux systems). If building from source, then start from
the release distributions of the source code, also provided on the
[releases page][releases]. Cloning/archiving the source
directly from git could result in missing version number entries,
undesirably stringent build requirements or an unstable development
version between releases. Additional build notes for methods developers can
be found in the [developer guide][developerGuide].

[releases]:https://github.com/Illumina/manta/releases
[DeveloperGuide]:../developerGuide/README.md


### Prerequisites to build from source

A C++11 compiler is required. These are the current minimum compiler versions enforced
by the build system:

* python 3.6+
* gcc 4.8+ OR clang 3.2+ (OR Visual Studio 2013+, see windows note below)
* cmake 3.5+
* libz (including headers)
* (optional) libbz2 and liblzma (including headers). When found, these are used to enable
  support for CRAM files written with the bzip2/lzma block codecs.

This fork is routinely built with gcc 13 and cmake 3.28.

### Runtime prerequisites

* python 3.6+

Note that this is a change from the original manta distribution, which required python 2.
The workflow and all helper scripts were ported to python 3 in this fork, and no longer
run under python 2.

### Operating System Guidelines

##### Linux

Manta is known to build and run on the following linux distributions
(with additional packages as described below):

- Ubuntu 12.04, 14.04, 16.04
- CentOS 6, 7

##### OS X

Manta builds and passes basic tests on OS X 10.9, but full WGS analyses
are not tested for this platform.

##### Windows

Manta does not build or run on windows. Library-level compilation is
possible for Visual Studio users. See the [developer guide][DeveloperGuide] for details.

### Linux Package Additions

##### Ubuntu 22.04 and 24.04

    apt-get update -qq
    apt-get install -qq bzip2 cmake gcc g++ make python3 zlib1g-dev libbz2-dev liblzma-dev

##### Rocky/Alma Linux 8 and 9

    dnf install -y tar bzip2 cmake make gcc gcc-c++ python3 zlib-devel bzip2-devel xz-devel

### Build procedure

After acquiring a release distribution of the source code, the build
procedure is:

* Unpack source code
* Create and move to a separate `build` directory (out-of-source build is required.)
* Configure build
* Compile & Install

Example (building on 4 cores):

    wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2
    tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2
    mkdir build && cd build
    # Ensure that CC and CXX are updated to target compiler if needed, e.g.:
    #     export CC=/path/to/cc
    #     export CXX=/path/to/c++
    ../manta-${MANTA_VERSION}.release_src/configure --jobs=4 --prefix=/path/to/install
    make -j4 install

Note that there is one other dependency, boost, which can optionally be provided by the user. It will
automatically be built from source if not detected. The minimum required version for users planning to provide
it to the build process is

* boost 1.58.0+ (must include static libraries)

...the build process will find an existing boost installation indicated by the environment variable `BOOST_ROOT`.

If not detected, then boost will be built from source and installed to a temporary location under
the build directory automatically. This step can make installation more convenient, but does increase the time
required for configuration. To accelerate this process the configuration step can be parallelized over multiple cores,
as demonstrated in the example above with the`--jobs=4` argument to configure.

htslib and samtools (version 1.25) are always built from the copies bundled in the `redist` directory.

To see more configure options, run:

    ${MANTA_SRC_PATH}/configure --help

##### Workflow relocation

After Manta is built the installation directory can be relocated to
another directory.  All internal paths used in the workflow are
relative.

### Demo

To help verify a successful installation, Manta includes a small demo
data set and test script. After completing the installation steps
above, the demo can be run as follows:

    python3 ${MANTA_INSTALL_PATH}/bin/runMantaWorkflowDemo.py

This script creates a `MantaDemoAnalysis` directory under the current
working directory, runs Manta on a small demo dataset, and compares the
somatic structural variant output to an expected result.

See [the demo README](../../src/demo/README.md) for additional information
on the test script and data.

