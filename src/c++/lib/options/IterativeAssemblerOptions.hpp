//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#pragma once

#include <string>

/// Input parameters for IterativeAssembler
///
struct IterativeAssemblerOptions {
  IterativeAssemblerOptions() {}

  /// the symbol set used during assembly
  std::string alphabet = "ACGT";

  /// minimum basecall quality for assembly input
  int minQval = 5;

  /// initial word (kmer) length
  unsigned minWordLength = 41;

  unsigned maxWordLength   = 76;
  unsigned wordStepSize    = 5;
  unsigned minContigLength = 15;

  /// min. coverage required for contig extension
  unsigned minCoverage = 1;

  /// coverage required for conservative contig sub-range
  unsigned minConservativeCoverage = 2;

  /// max error rates allowed during contig extension
  double maxError = 0.35;

  /// min. number of unused reads to enable search for more contigs
  unsigned minUnusedReads = 3;

  /// min. number of reads required to start assembly
  unsigned minSupportReads = 2;

  /// Max. number of assembly returned for a given set of reads
  unsigned maxAssemblyCount = 10;

  /// If true, the conservative sub-range of each returned contig is recomputed after assembly by
  /// placing the input reads back onto the finished contig while tolerating a low mismatch rate.
  ///
  /// During contig extension the conservative range is ended by the first k-mer which is observed in
  /// fewer than minConservativeCoverage reads, and it is never re-opened. Because k-mers have to match
  /// exactly, a single basecall error in one of a small number of reads is enough to do this, which
  /// regularly excludes most of a contig, including sequence which every read agrees on.
  bool isRefineConservativeRange = true;

  /// Length of the exact match used to seed the placement of a read on a contig
  unsigned refineSeedLength = 20;

  /// Minimum overlap between a read and a contig for the read to count towards contig coverage
  unsigned refineMinOverlap = 40;

  /// Maximum fraction of mismatching bases in the overlap of a read and a contig for the read to count
  /// towards contig coverage
  double refineMaxMismatchFrac = 0.04;
};
