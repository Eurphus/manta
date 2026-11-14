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

#include "options/CallOptionsDiploid.hpp"

boost::program_options::options_description getOptionsDescription(CallOptionsDiploid& opt)
{
  namespace po = boost::program_options;
  po::options_description desc("germline-variant-calling");
  // clang-format off
  desc.add_options()
  ("diploid-max-depth-factor", po::value(&opt.maxDepthFactor)->default_value(opt.maxDepthFactor),
   "Variants where the depth around the breakpoint is greater than this factor x the chromosomal mean will be filtered out")
  ("min-qual-score", po::value(&opt.minOutputAltScore)->default_value(opt.minOutputAltScore),
   "minimum QUAL score for variants included in the germline output vcf")
  ("min-pass-qual-score", po::value(&opt.minPassAltScore)->default_value(opt.minPassAltScore),
   "minimum QUAL score for variants to PASS in germline output vcf")
  ("min-pass-gt-score", po::value(&opt.minPassGTScore)->default_value(opt.minPassGTScore),
   "minimum genotype quality score below which samples are filtered for a variant in the germline output vcf")
  ("no-tandem-dup-model", po::value(&opt.isDisableTandemDupModel)->zero_tokens(),
   "Score the evidence for tandem duplications exactly as for all other variant types (v1.6.0 behavior). By default the germline scoring of a tandem duplication accounts for the reference junctions which remain on a haplotype carrying it")
  ("no-depth-ratio", po::value(&opt.isDisableDepthRatio)->zero_tokens(),
   "Do not report the read depth within deletions and tandem duplications relative to their flanking regions (FORMAT/RDR)")
  ("min-depth-ratio-sv-size", po::value(&opt.minDepthRatioSVSize)->default_value(opt.minDepthRatioSVSize),
   "Smallest deletion or tandem duplication for which the depth ratio is reported")
  ("filter-uncertain-zygosity", po::value(&opt.isFilterUncertainZygosity)->zero_tokens(),
   "Apply the record level sample filter to a variant if its only non-reference samples have low genotype quality, even if the low quality is due to uncertain zygosity and not to an uncertain presence of the variant (v1.6.0 behavior)")
  ;
  // clang-format on

  return desc;
}
