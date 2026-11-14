#!/usr/bin/env bash
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2019 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

################################################################################
##
## Installation script for cmake
##
## author Come Raczy
##
################################################################################

set -o nounset
set -o pipefail

REDIST_DIR=$1
INSTALL_DIR=$2
if [[ $# -ge 3 ]] ; then PARALLEL=$3 ; else PARALLEL=1 ; fi

script_dir="$(dirname "$0")"
source $script_dir/common.bash

BUILD_DIR=${INSTALL_DIR}/build
BIN_DIR=${INSTALL_DIR}/bin
INCLUDE_DIR=${INSTALL_DIR}/include

# cmake 3.5 is the minimum version supported by the project CMakeLists.txt
CMAKE_MAJOR=3
CMAKE_MINOR=5
CMAKE_PATCH_MIN=0
CMAKE_REQUIRED="$CMAKE_MAJOR.$CMAKE_MINOR.$CMAKE_PATCH_MIN"

SCRIPT=`basename "$0"`

AVAILABLE_CMAKE_VERSION=`cmake --version 2> /dev/null`
if [[ "${AVAILABLE_CMAKE_VERSION}" =~ ^cmake\ version\ ([0-9]+)\.([0-9]+)\.([0-9]+) ]] ; then
    MAJOR=${BASH_REMATCH[1]}
    MINOR=${BASH_REMATCH[2]}
    PATCH=${BASH_REMATCH[3]}
    if [[ "$MAJOR" -gt "$CMAKE_MAJOR" || ( "$MAJOR" -eq "$CMAKE_MAJOR" && ( "$MINOR" -gt "$CMAKE_MINOR" || ( "$MINOR" -eq "$CMAKE_MINOR" && "$PATCH" -ge "$CMAKE_PATCH_MIN"  ) ) ) ]] ; then
        ilog "${BASH_REMATCH[0]} (>= $CMAKE_REQUIRED) is already installed"
        echo "cmake"
        exit 0
    fi
fi

# This fork no longer bundles a cmake source distribution to bootstrap from (the previously
# bundled cmake 2.8.12 is older than the minimum required version), so a system cmake is required.
ilog "ERROR: cmake $CMAKE_REQUIRED or newer is required, but was not found on PATH."
ilog "\tPlease install cmake (version $CMAKE_REQUIRED+) and rerun configure."
exit 1
