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

#include "alignment/AlignmentScores.hpp"
#include "options/IterativeAssemblerOptions.hpp"

typedef IterativeAssemblerOptions AssemblerOptions;

namespace ALIGNERPARAM {
static const int largeGapOpenScore(-100);
}

/// Options for the SV refiner step
///
/// Note that we have two categories of options for assembly and alignment,
/// one for small events, and one for large events
///
struct SVRefinerOptions {
  /// match, mismatch, open score ratios taken from bwa defaults (but not extend!) :
  ///
  SVRefinerOptions()
    : largeSVAlignScores(2, -8, -24, -1, -1),
      largeInsertEdgeAlignScores(2, -8, -18, -1, -1),
      largeInsertCompleteAlignScores(2, -8, ALIGNERPARAM::largeGapOpenScore, 0, -1),
      spanningAlignScores(2, -8, -12, -1, -1),
      largeGapOpenScore(ALIGNERPARAM::largeGapOpenScore),
      jumpScore(-100),
      contigFilterScores(2, -8, -18, 0, -1)
  {
    // For breakend-spanning assemblies we require a larger contig than for small-variant assemblies
    spanningAssembleOpt.minContigLength = 75;
  }

  // parameters for small SV assembly/alignment:
  AlignmentScores<int> largeSVAlignScores;  // large SV but at a single assembly locus
  AlignmentScores<int> largeInsertEdgeAlignScores;
  AlignmentScores<int> largeInsertCompleteAlignScores;
  AssemblerOptions     smallSVAssembleOpt;

  // parameters for large SV assembly/alignment:
  AlignmentScores<int> spanningAlignScores;
  const int            largeGapOpenScore;
  const int            jumpScore;
  AssemblerOptions     spanningAssembleOpt;

  /// Scores used to check for low quality assembly contig alignments
  AlignmentScores<int> contigFilterScores;

  /// Maximum number of contigs from one local (single locus) assembly which can contribute SV candidates.
  ///
  /// A value of 1 reproduces the v1.6.0 behavior, where only the single best contig alignment is used, so
  /// that at most one alternate haplotype can be found at a locus. Higher values allow additional contigs to
  /// add candidates which differ from those already found, such as the second allele of a compound
  /// heterozygous site.
  unsigned maxSmallSVHaplotypeCount = 2;

  /// Minimum number of supporting reads for a contig to add candidates beyond those of the best contig
  unsigned minSecondaryHaplotypeSupportReads = 2;

  /// Two candidates from the same local assembly are treated as the same allele if their breakend positions
  /// and insertion lengths all agree to within this many bases
  int sameAlleleTolerance = 3;

  /// If true, a left-anchored and a right-anchored large insertion contig are joined into one contig when the
  /// insertion sequence at the end of one overlaps that at the start of the other. This yields a fully
  /// assembled insertion where the assembler could not join the two sides itself because the overlap is
  /// shorter than its minimum word length.
  bool isMergeLargeInsertionContigs = true;

  /// Minimum overlap required to join two large insertion contigs
  unsigned minInsertionContigOverlap = 20;

  /// Maximum mismatch fraction in the overlap used to join two large insertion contigs
  double maxInsertionContigOverlapMismatchFrac = 0.05;
};
