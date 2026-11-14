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

#include "boost/program_options.hpp"

struct CallOptionsDiploid {
  float indelPrior = 1e-5f;

  /// breakpoints where the depth is greater than the chromosome average x this factor are filtered
  /// out:
  float       maxDepthFactor      = 3.0f;
  std::string maxDepthFilterLabel = "MaxDepth";

  /// minimum QUAL score to print out a diploid variant
  unsigned minOutputAltScore = 10;

  /// minimum QUAL score to PASS a diploid variant
  unsigned    minPassAltScore   = 20;
  std::string minAltFilterLabel = "MinQUAL";

  // control filtration based on MQ0 fraction:
  float       maxMQ0Frac      = 0.4f;
  std::string maxMQ0FracLabel = "MaxMQ0Frac";

  /// filter for large SVs with no pair support
  std::string noPairSupportLabel = "NoPairSupport";

  /// no sample passes all sample-specific filters
  std::string failedSampleFTLabel = "SampleFT";

  /// below this GQ value, the SAMPLE filter is marked in the VCF
  unsigned    minPassGTScore   = 15;
  std::string minGTFilterLabel = "MinGQ";

  /// the SAMPLE filter for home-ref calls
  std::string homRefLabel = "HomRef";

  /// If true, the evidence for a tandem duplication is scored in a way which accounts for the structure of a
  /// haplotype carrying it. This applies to germline scoring of a duplication called as a single junction:
  ///
  /// 1. The choice between a heterozygous and homozygous genotype is made using the allele fractions expected
  /// for a duplication instead of those expected for other variant types. The sequence across both reference
  /// breakends of a tandem duplication is still present on a haplotype carrying the duplication, so reads
  /// supporting the reference allele are observed even for a homozygous duplication, see
  /// getTandemDupAlleleFractions(). The probability that the sample carries the duplication at all does not
  /// depend on these fractions.
  ///
  /// 2. If one read of a fragment crosses the duplication junction and the other read crosses a reference
  /// breakend, the fragment is scored on the read crossing the duplication junction. Both are present on one
  /// molecule from a haplotype carrying the duplication, so the read at the reference breakend is not evidence
  /// against the duplication, see getRefAltFragLnLhood().
  ///
  /// 3. Paired-read evidence is phased out for small duplications, in the same way as it always has been for
  /// small deletions and insertions. Read pairs supporting a tandem duplication are only found through their
  /// anomalous orientation, which requires a duplication at least as large as the sequenced fragments. For a
  /// smaller duplication paired-read evidence can only ever be found for the reference allele.
  ///
  /// If false, tandem duplications are scored exactly as all other variant types (v1.6.0 behavior)
  bool isTandemDupModel = true;

  /// Set by the command line to turn off isTandemDupModel
  bool isDisableTandemDupModel = false;

  /// If true, the read depth within deletions and tandem duplications is reported relative to the depth of
  /// their flanking regions, for each sample (FORMAT/RDR).
  ///
  /// The split and paired read evidence used to call a variant says nothing about the amount of sequence
  /// between its breakends. A real deletion or duplication also changes the read depth between them, so this
  /// value is an independent check on the call. It is reported only, and has no influence on any score or
  /// filter.
  bool isReportDepthRatio = true;

  /// Set by the command line to turn off isReportDepthRatio
  bool isDisableDepthRatio = false;

  /// Smallest deletion or tandem duplication for which the depth ratio is reported
  unsigned minDepthRatioSVSize = 500;

  /// If false, a variant is not filtered at the record level when its only confident sample fails the minimum
  /// genotype quality filter because the sample's zygosity is uncertain.
  ///
  /// The genotype quality is low whenever two genotypes are similarly likely. If those two genotypes are
  /// heterozygous and homozygous, there is no doubt that the sample carries the variant. Such a sample is
  /// still marked with the genotype quality filter at the sample level, but it is sufficient to keep the
  /// record level sample filter from being applied, if the probability of the reference genotype is low
  /// enough to pass the genotype quality threshold.
  ///
  /// If true, any sample failing the minimum genotype quality filter counts as a failed sample (v1.6.0
  /// behavior).
  bool isFilterUncertainZygosity = false;
};

boost::program_options::options_description getOptionsDescription(CallOptionsDiploid& opt);
