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

/// \file
/// \author Chris Saunders and Xiaoyu Chen
///

#include "SVScorer.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "SVScorePairAltProcessor.hpp"
#include "SplitReadAlignment.hpp"
#include "blt_util/LinearScaler.hpp"
#include "blt_util/log.hpp"
#include "blt_util/math_util.hpp"
#include "blt_util/prob_util.hpp"
#include "blt_util/qscore.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/align_path_bam_util.hpp"
#include "manta/BamStreamerUtils.hpp"
#include "manta/ReadGroupStatsSet.hpp"
#include "manta/SVCandidateUtil.hpp"

//#define DEBUG_SCORE
//#define DEBUG_SOMATIC_SCORE

#if defined(DEBUG_SCORE) || defined(DEBUG_SOMATIC_SCORE)
#define ANY_DEBUG_SCORE
#endif

SVScorer::SVScorer(const GSCOptions& opt, const SVLocusScanner& readScanner, const bam_header_info& header)
  : _callOpt(opt.callOpt),
    _callDopt(_callOpt),
    _diploidOpt(opt.diploidOpt),
    _diploidDopt(_diploidOpt),
    _scanOpt(opt.scanOpt),
    _refineOpt(opt.refineOpt),
    _dFilterDiploid(opt.chromDepthFilename, _diploidOpt.maxDepthFactor, header),
    _readScanner(readScanner),
    _header(header)
{
  openBamStreams(opt.referenceFilename, opt.alignFileOpt.alignmentFilenames, _bamStreams);

  _sampleCount = opt.alignFileOpt.alignmentFilenames.size();
}

/// add bam alignment to simple short-range vector depth estimate
///
/// \param[in] beginPos this is the begin position of the range covered by the depth array
///
static void addReadToDepthEst(const bam_record& bamRead, const pos_t beginPos, std::vector<unsigned>& depth)
{
  using namespace ALIGNPATH;

  const pos_t endPos(beginPos + depth.size());

  // get cigar:
  path_t apath;
  bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);

  pos_t refPos(bamRead.pos() - 1);
  for (const path_segment& ps : apath) {
    if (refPos >= endPos) return;

    if (is_segment_align_match(ps.type)) {
      for (pos_t pos(refPos); pos < (refPos + static_cast<pos_t>(ps.length)); ++pos) {
        if (pos >= beginPos) {
          if (pos >= endPos) return;
          depth[pos - beginPos]++;
        }
      }
    }
    if (is_segment_type_ref_length(ps.type)) refPos += ps.length;
  }
}

void SVScorer::getBreakendMaxMappedDepthAndMQ0(
    const bool        isMaxDepth,
    const double      cutoffDepth,
    const SVBreakend& bp,
    unsigned&         maxDepth,
    float&            MQ0Frac)
{
  /// define a new interval -/+ 50 bases around the center pos
  /// of the breakpoint
  static const pos_t regionSize(50);

  maxDepth = 0;
  MQ0Frac  = 0;

  unsigned totalReads(0);
  unsigned totalMQ0Reads(0);

  const pos_t            centerPos(bp.interval.range.center_pos());
  const known_pos_range2 searchRange(std::max((centerPos - regionSize), 0), (centerPos + regionSize));

  if (searchRange.size() == 0) return;

  std::vector<unsigned> depth(searchRange.size(), 0);

  bool isCutoff(false);
  bool isBamFound(false);

  const unsigned bamCount(_bamStreams.size());
  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    isBamFound = true;

    bam_streamer& bamStream(*_bamStreams[bamIndex]);

    // set bam stream to new search interval:
    bamStream.resetRegion(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

    while (bamStream.next()) {
      const bam_record& bamRead(*(bamStream.get_record_ptr()));

      const pos_t refPos(bamRead.pos() - 1);
      if (refPos >= searchRange.end_pos()) break;

      if (isReadUnmappedOrFilteredCore(bamRead)) continue;

      addReadToDepthEst(bamRead, searchRange.begin_pos(), depth);

      totalReads++;
      if (0 == bamRead.map_qual()) totalMQ0Reads++;

      if (isMaxDepth) {
        const pos_t depthOffset(refPos - searchRange.begin_pos());
        if (depthOffset >= 0) {
          if (depth[depthOffset] > cutoffDepth) {
            isCutoff = true;
            break;
          }
        }
      }
    }

    if (isCutoff) break;
  }

  assert(isBamFound);

  maxDepth = *(std::max_element(depth.begin(), depth.end()));
  if (totalReads >= 10) {
    MQ0Frac = static_cast<float>(totalMQ0Reads) / static_cast<float>(totalReads);
  }
}

/// \return The number of reference bases in \p range which are covered by an aligned (match) segment of a
/// read alignment starting at \p readBeginPos
static unsigned getAlignmentRangeOverlap(
    const ALIGNPATH::path_t& apath, const pos_t readBeginPos, const known_pos_range2& range)
{
  using namespace ALIGNPATH;

  pos_t overlap(0);
  pos_t refPos(readBeginPos);
  for (const path_segment& ps : apath) {
    if (refPos >= range.end_pos()) break;
    const pos_t segmentEndPos(refPos + static_cast<pos_t>(ps.length));
    if (is_segment_align_match(ps.type)) {
      const pos_t overlapBegin(std::max(refPos, range.begin_pos()));
      const pos_t overlapEnd(std::min(segmentEndPos, range.end_pos()));
      if (overlapEnd > overlapBegin) overlap += (overlapEnd - overlapBegin);
    }
    if (is_segment_type_ref_length(ps.type)) refPos = segmentEndPos;
  }
  return overlap;
}

/// Regions used to compare the read depth within a variant to the depth next to it
struct SVDepthRegions {
  std::vector<known_pos_range2> variant;
  std::vector<known_pos_range2> flank;
};

/// Choose the regions used to estimate the read depth within and next to a deletion or duplication
///
/// The variant depth is taken from the whole interval between the breakends, leaving out a margin at each end
/// where reads are clipped or split at the breakend. For a large variant a fixed number of evenly spaced
/// windows is used instead of the whole interval, which bounds the cost of the estimate for any variant size.
/// The flanking depth is taken from one region on each side of the variant, separated from it by a gap.
///
/// \param[in] svBeginPos First reference position within the variant
///
/// \param[in] svEndPos Position following the last reference position within the variant
///
static SVDepthRegions getSVDepthRegions(const pos_t svBeginPos, const pos_t svEndPos, const pos_t chromSize)
{
  static const pos_t flankSize(1000);
  static const pos_t minFlankSize(200);
  static const pos_t flankGap(100);
  static const pos_t maxEdgeTrim(100);
  static const pos_t maxFullScanSize(20000);
  static const pos_t sampleWindowSize(1000);
  static const pos_t sampleWindowCount(20);

  SVDepthRegions regions;

  const pos_t svSize(svEndPos - svBeginPos);
  if ((svSize <= 0) || (svBeginPos < 0) || (svEndPos > chromSize)) return regions;

  const pos_t edgeTrim(std::min(maxEdgeTrim, svSize / 10));
  const pos_t innerBeginPos(svBeginPos + edgeTrim);
  const pos_t innerEndPos(svEndPos - edgeTrim);
  const pos_t innerSize(innerEndPos - innerBeginPos);
  if (innerSize <= 0) return regions;

  if (innerSize <= maxFullScanSize) {
    regions.variant.emplace_back(innerBeginPos, innerEndPos);
  } else {
    const double windowStep(
        static_cast<double>(innerSize - sampleWindowSize) / static_cast<double>(sampleWindowCount - 1));
    for (pos_t windowIndex(0); windowIndex < sampleWindowCount; ++windowIndex) {
      const pos_t windowBeginPos(innerBeginPos + static_cast<pos_t>(windowIndex * windowStep));
      regions.variant.emplace_back(windowBeginPos, windowBeginPos + sampleWindowSize);
    }
  }

  {
    const pos_t flankEndPos(svBeginPos - flankGap);
    const pos_t flankBeginPos(std::max(0, flankEndPos - flankSize));
    if ((flankEndPos - flankBeginPos) >= minFlankSize) regions.flank.emplace_back(flankBeginPos, flankEndPos);
  }
  {
    const pos_t flankBeginPos(svEndPos + flankGap);
    const pos_t flankEndPos(std::min(chromSize, flankBeginPos + flankSize));
    if ((flankEndPos - flankBeginPos) >= minFlankSize) regions.flank.emplace_back(flankBeginPos, flankEndPos);
  }

  return regions;
}

double SVScorer::getMeanRegionDepth(
    const unsigned bamIndex, const int32_t tid, const std::vector<known_pos_range2>& regions)
{
  double alignedBaseCount(0);
  double regionBaseCount(0);

  bam_streamer&     bamStream(*_bamStreams[bamIndex]);
  ALIGNPATH::path_t apath;
  for (const known_pos_range2& region : regions) {
    if (region.end_pos() <= region.begin_pos()) continue;
    regionBaseCount += (region.end_pos() - region.begin_pos());

    bamStream.resetRegion(tid, region.begin_pos(), region.end_pos());
    while (bamStream.next()) {
      const bam_record& bamRead(*(bamStream.get_record_ptr()));

      const pos_t refPos(bamRead.pos() - 1);
      if (refPos >= region.end_pos()) break;

      if (isReadUnmappedOrFilteredCore(bamRead)) continue;

      bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);
      alignedBaseCount += getAlignmentRangeOverlap(apath, refPos, region);
    }
  }

  if (regionBaseCount <= 0) return -1;
  return (alignedBaseCount / regionBaseCount);
}

void SVScorer::getSVDepthSupport(const SVCandidate& sv, SVScoreInfo& baseInfo)
{
  if (sv.bp1.interval.tid != sv.bp2.interval.tid) return;

  const EXTENDED_SV_TYPE::index_t svType(getExtendedSVType(sv));
  if (!((svType == EXTENDED_SV_TYPE::DELETE) || (svType == EXTENDED_SV_TYPE::TANDUP))) return;

  if (isSVBelowMinSize(sv, _diploidOpt.minDepthRatioSVSize)) return;

  // A region at many times the expected depth is not a place where a depth ratio means much, and is expensive
  // to scan. Evidence enumeration is cut off for the same reason.
  const int32_t tid(sv.bp1.interval.tid);
  if (_dFilterDiploid.isMaxDepthFilter()) {
    const double maxDepth(_dFilterDiploid.maxDepth(tid));
    if ((baseInfo.bp1MaxDepth > maxDepth) || (baseInfo.bp2MaxDepth > maxDepth)) return;
  }

  // the variant covers the reference positions following the first breakend, up to and including that of the
  // second breakend:
  const pos_t bp1Pos(sv.bp1.interval.range.center_pos());
  const pos_t bp2Pos(sv.bp2.interval.range.center_pos());
  const pos_t svBeginPos(std::min(bp1Pos, bp2Pos) + 1);
  const pos_t svEndPos(std::max(bp1Pos, bp2Pos) + 1);

  assert(tid < static_cast<int32_t>(_header.chrom_data.size()));
  const pos_t          chromSize(_header.chrom_data[tid].length);
  const SVDepthRegions regions(getSVDepthRegions(svBeginPos, svEndPos, chromSize));
  if (regions.variant.empty() || regions.flank.empty()) return;

  const unsigned bamCount(_bamStreams.size());
  assert(baseInfo.samples.size() == bamCount);
  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    const double flankDepth(getMeanRegionDepth(bamIndex, tid, regions.flank));
    if (flankDepth <= 0) continue;
    const double variantDepth(getMeanRegionDepth(bamIndex, tid, regions.variant));
    if (variantDepth < 0) continue;

    SVSampleInfo& sampleInfo(baseInfo.samples[bamIndex]);
    sampleInfo.isDepthEstimated = true;
    sampleInfo.variantDepth     = static_cast<float>(variantDepth);
    sampleInfo.flankDepth       = static_cast<float>(flankDepth);
  }
}

/// Convert log likelihoods from a 2 state space to probabilities:
static void lnToProb(float& lower, float& higher)
{
  lower  = std::exp(lower - higher);
  higher = 1 / (lower + 1);
  lower  = lower / (lower + 1);
}

/// Get ref and alt allele likelihoods for the given read
///
/// \param isForcedSupport If true, proceed with the likelihood calculation even if there is no split support
/// indicated for this read for any breakend of any allele.
///
/// \return False if likelihoods were not computed
static bool getSampleSplitReadLnLhood(
    const SVFragmentEvidence& fragev,
    const bool                isRead1,
    float&                    refLnLhood,
    float&                    altLnLhood,
    const bool                isForcedSupport = false)
{
  refLnLhood = 1.;
  altLnLhood = 1.;

  const std::pair<bool, bool> isBpSupport(fragev.isAnySplitReadSupport(isRead1));
  if (!isForcedSupport) {
    // Abort the process if the read doesn't support any allele breakends
    if (!(isBpSupport.first || isBpSupport.second)) return false;
  }

  bool isUseBp1Score(isBpSupport.first);

  if (isForcedSupport || (isBpSupport.first == isBpSupport.second)) {
    isUseBp1Score =
        (fragev.alt.bp1.getRead(isRead1).splitLnLhood >= fragev.alt.bp2.getRead(isRead1).splitLnLhood);
  }

  altLnLhood =
      (isUseBp1Score ? fragev.alt.bp1.getRead(isRead1).splitLnLhood
                     : fragev.alt.bp2.getRead(isRead1).splitLnLhood);

  if (isBpSupport.first && isBpSupport.second) {
    isUseBp1Score =
        (fragev.ref.bp1.getRead(isRead1).splitLnLhood >= fragev.ref.bp2.getRead(isRead1).splitLnLhood);
  }

  refLnLhood =
      (isUseBp1Score ? fragev.ref.bp1.getRead(isRead1).splitLnLhood
                     : fragev.ref.bp2.getRead(isRead1).splitLnLhood);

  return true;
}

static void addConservativeSplitReadSupport(
    const SVFragmentEvidence& fragev, const bool isRead1, SVSampleInfo& sampleBaseInfo)
{
  // A read must support a single allele with at least this probability to be enumerated in the VCF output's
  // split read count:
  static const float splitSupportProb(0.999f);

  // only consider reads where at least one allele and one breakend is confident
  //
  // ...note this is done in the absence of having a noise state in the model
  //
  float refLnLhood;
  float altLnLhood;
  if (!getSampleSplitReadLnLhood(fragev, isRead1, refLnLhood, altLnLhood)) return;

  // convert to normalized prob:
  if (altLnLhood > refLnLhood) {
    lnToProb(refLnLhood, altLnLhood);
    if (altLnLhood > splitSupportProb) sampleBaseInfo.alt.confidentSplitReadCount++;

  } else {
    lnToProb(altLnLhood, refLnLhood);
    if (refLnLhood > splitSupportProb) {
      sampleBaseInfo.ref.confidentSplitReadCount++;
      if (fragev.ref.bp1.getRead(isRead1).isSplitSupport)
        sampleBaseInfo.ref.confidentSplitReadAndPairCountRefBp1++;
      if (fragev.ref.bp2.getRead(isRead1).isSplitSupport)
        sampleBaseInfo.ref.confidentSplitReadAndPairCountRefBp2++;
    }
  }
}

static float getSpanningPairAlleleLhood(const SVFragmentEvidenceAllele& allele)
{
  float fragProb(0);
  if (allele.bp1.isFragmentSupport) {
    fragProb = allele.bp1.fragLengthProb;
  }

  if (allele.bp2.isFragmentSupport) {
    fragProb = std::max(fragProb, allele.bp2.fragLengthProb);
  }

  return fragProb;
}

static void addSpanningPairSupport(const SVFragmentEvidence& fragev, SVSampleInfo& sampleBaseInfo)
{
  if (fragev.alt.bp1.isFragmentSupport || fragev.alt.bp2.isFragmentSupport) {
    sampleBaseInfo.alt.spanningPairCount++;
  }
  if (fragev.ref.bp1.isFragmentSupport || fragev.ref.bp2.isFragmentSupport) {
    sampleBaseInfo.ref.spanningPairCount++;
  }
}

static void addConservativeSpanningPairSupport(const SVFragmentEvidence& fragev, SVSampleInfo& sampleBaseInfo)
{
  static const float pairSupportProb(0.9f);

  if (!fragev.isAnySpanningPairSupport()) return;

  float altLhood(getSpanningPairAlleleLhood(fragev.alt));
  float refLhood(getSpanningPairAlleleLhood(fragev.ref));

  assert(altLhood >= 0);
  assert(refLhood >= 0);
  if ((altLhood <= 0) && (refLhood <= 0)) {
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "Spanning likelihood is zero for all alleles. Fragment: " << fragev;
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
  }

  static const bool isTier2(false);
  const bool isFullyMapped(fragev.read1.isObservedAnchor(isTier2) && fragev.read2.isObservedAnchor(isTier2));

  // convert to normalized prob:
  const float sum(altLhood + refLhood);
  if (altLhood > refLhood) {
    if ((altLhood / sum) > pairSupportProb) {
#ifdef DEBUG_SCORE
      log_os << __FUNCTION__ << ": semi-mapped alt pair support\n";
#endif
      sampleBaseInfo.alt.confidentSemiMappedSpanningPairCount++;
      if (isFullyMapped) {
#ifdef DEBUG_SCORE
        log_os << __FUNCTION__ << ": fully-mapped alt pair support\n";
#endif
        sampleBaseInfo.alt.confidentSpanningPairCount++;
      }
    }

  } else {
    if ((refLhood / sum) > pairSupportProb) {
#ifdef DEBUG_SCORE
      log_os << __FUNCTION__ << ": semi-mapped ref pair support\n";
#endif
      sampleBaseInfo.ref.confidentSemiMappedSpanningPairCount++;
      if (isFullyMapped) {
#ifdef DEBUG_SCORE
        log_os << __FUNCTION__ << ": fully-mapped ref pair support\n";
#endif
        sampleBaseInfo.ref.confidentSpanningPairCount++;
        if (fragev.ref.bp1.isFragmentSupport) sampleBaseInfo.ref.confidentSplitReadAndPairCountRefBp1++;
        if (fragev.ref.bp2.isFragmentSupport) sampleBaseInfo.ref.confidentSplitReadAndPairCountRefBp2++;
      }
    }
  }
}

static void getSampleCounts(const SVEvidence::evidenceTrack_t& sampleEvidence, SVSampleInfo& sampleBaseInfo)
{
  for (const SVEvidence::evidenceTrack_t::value_type& val : sampleEvidence) {
    const SVFragmentEvidence& fragev(val.second);
#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": Counting read: " << val.first << "\n";
#endif
    // evaluate read1 and read2 from this fragment
    //
    addConservativeSplitReadSupport(fragev, true, sampleBaseInfo);
    addConservativeSplitReadSupport(fragev, false, sampleBaseInfo);
    addSpanningPairSupport(fragev, sampleBaseInfo);
    addConservativeSpanningPairSupport(fragev, sampleBaseInfo);
  }
}

/// get conservative count of reads which support only one allele, ie. P ( allele | read ) is high
///
static void getSVSupportSummary(const SVEvidence& evidence, SVScoreInfo& baseInfo)
{
  const unsigned sampleCount(baseInfo.samples.size());
  assert(sampleCount == evidence.samples.size());

  for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex) {
    getSampleCounts(evidence.getSampleEvidence(sampleIndex), baseInfo.samples[sampleIndex]);
  }
}

/// Check for and fix conflicts in pair and split read support for one sample
///
/// \param isFindAltPairConflict If true, resolve conflicts in support for the alt allele (ref allele
/// conflicts are always resolved)
static void resolvePairSplitConflictsSample(
    const bool isFindAltPairConflict, SVEvidence::evidenceTrack_t& sampleEvidence)
{
  for (SVEvidence::evidenceTrack_t::value_type& val : sampleEvidence) {
#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": conflict check for read name " << val.first << "\n";
#endif
    SVFragmentEvidence& fragev(val.second);

    // This filtration scheme is only relevant if there's pair and split support for the same fragment:
    if (!fragev.isAnySpanningPairSupport()) continue;
    //    if (! (fragev.isAnySplitReadSupport(true) || fragev.isAnySplitReadSupport(false))) continue;

    // If there's a difference in fragment support for one allele, then
    // there must also be either neutral split support or split support in favor of the same allele

    const float refPairLhood(getSpanningPairAlleleLhood(fragev.ref));
    const float altPairLhood(getSpanningPairAlleleLhood(fragev.alt));

    static const bool isForcedSupport(true);

    float      refSplitLnLhoodRead1;
    float      altSplitLnLhoodRead1;
    const bool isRead1Split(
        getSampleSplitReadLnLhood(fragev, true, refSplitLnLhoodRead1, altSplitLnLhoodRead1, isForcedSupport));

    float      refSplitLnLhoodRead2;
    float      altSplitLnLhoodRead2;
    const bool isRead2Split(getSampleSplitReadLnLhood(
        fragev, false, refSplitLnLhoodRead2, altSplitLnLhoodRead2, isForcedSupport));

#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": fragev " << fragev << "\n";
    log_os << __FUNCTION__ << ": ref/alt pair " << refPairLhood << " " << altPairLhood << "\n";
    log_os << __FUNCTION__ << ": r1 issplit/refLhood/altLhood " << isRead1Split << " " << refSplitLnLhoodRead1
           << " " << altSplitLnLhoodRead1 << "\n";
    log_os << __FUNCTION__ << ": r2 issplit/refLhood/altLhood " << isRead2Split << " " << refSplitLnLhoodRead2
           << " " << altSplitLnLhoodRead2 << "\n";
#endif

    const bool isRefPair(refPairLhood > altPairLhood);
    const bool isAltPair(altPairLhood > refPairLhood);

    if (isAltPair) {
      if (!isFindAltPairConflict) continue;
    }

    if (isRead1Split) {
      if (altSplitLnLhoodRead1 > refSplitLnLhoodRead1) {
        if (isRefPair) {
#ifdef DEBUG_SCORE
          log_os << __FUNCTION__ << ": clearing alt1/ref\n";
#endif
          fragev.clearPairSupport();
        }
      }
      if (refSplitLnLhoodRead1 > altSplitLnLhoodRead1) {
        if (isAltPair) {
#ifdef DEBUG_SCORE
          log_os << __FUNCTION__ << ": clearing ref1/alt\n";
#endif
          fragev.clearPairSupport();
        }
      }
    }

    if (isRead2Split) {
      if (altSplitLnLhoodRead2 > refSplitLnLhoodRead2) {
        if (isRefPair) {
#ifdef DEBUG_SCORE
          log_os << __FUNCTION__ << ": clearing alt2/ref\n";
#endif
          fragev.clearPairSupport();
        }
      }
      if (refSplitLnLhoodRead2 > altSplitLnLhoodRead2) {
        if (isAltPair) {
#ifdef DEBUG_SCORE
          log_os << __FUNCTION__ << ": clearing ref2/alt\n";
#endif
          fragev.clearPairSupport();
        }
      }
    }
  }
}

/// Check for and fix conflicts in pair and split read support in all samples
///
/// Check for cases where pair support was added in error, the fragment does span the breakpoint, but the
/// alignment past the breakpoint is poor, and is better in the other allele.
///
/// A schematic of the intended filtration scenario is as follows:
///
/// >>X>>XXXX>XXXXX>>>read1>-------<read2<<<<<<<<<<<<<<<
///                |
///                Breakpoint Location
///
/// Here 'X' shows mismatch locations in the read, and thus poor split read support for this breakend, even
/// though it is possible that read pair analysis showed this to be supporting.
///
/// note this might be done more naturally during the pair computation, but all the info we need is added
/// during the split routine, so it's at least checked here as well:
static void resolvePairSplitConflicts(const SVCandidate& sv, SVEvidence& evidence)
{
  if (sv.isImprecise()) return;

  // Only find conflicts in pairs supporting the ALT allele for SVs below this size:
  static const pos_t maxAltPairConflictSearch(1000);
  const bool         isFindAltPairConflict(sv.centerSize() <= maxAltPairConflictSearch);

  const unsigned sampleCount(evidence.size());

  for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex) {
    resolvePairSplitConflictsSample(isFindAltPairConflict, evidence.getSampleEvidence(sampleIndex));
  }
}

void SVScorer::getSVSupportingEvidence(
    const SVCandidateSetData&      svData,
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate&             sv,
    const SVId&                    svId,
    SVScoreInfo&                   baseInfo,
    SVEvidence&                    evidence,
    SVEvidenceWriterData&          svSupports)
{
  // at what factor above the maxDepth FILTER criteria do we stop enumerating scoring components?
  static const unsigned cutoffDepthFactor(2);

  const bool isMaxDepth(_dFilterDiploid.isMaxDepthFilter());

  double bp1CutoffDepth(0);
  double bp2CutoffDepth(0);
  if (isMaxDepth) {
    const double bp1MaxMaxDepth(_dFilterDiploid.maxDepth(sv.bp1.interval.tid));
    const double bp2MaxMaxDepth(_dFilterDiploid.maxDepth(sv.bp2.interval.tid));

    bp1CutoffDepth = cutoffDepthFactor * bp1MaxMaxDepth;
    bp2CutoffDepth = cutoffDepthFactor * bp2MaxMaxDepth;
  }

  // get breakend center_pos depth estimate:
  getBreakendMaxMappedDepthAndMQ0(
      isMaxDepth, bp1CutoffDepth, sv.bp1, baseInfo.bp1MaxDepth, baseInfo.bp1MQ0Frac);
  const bool isBp1OverDepth(baseInfo.bp1MaxDepth > bp1CutoffDepth);
  if (!(isMaxDepth && isBp1OverDepth)) {
    getBreakendMaxMappedDepthAndMQ0(
        isMaxDepth, bp2CutoffDepth, sv.bp2, baseInfo.bp2MaxDepth, baseInfo.bp2MQ0Frac);
  }
  const bool isBp2OverDepth(baseInfo.bp2MaxDepth > bp2CutoffDepth);
  const bool isOverDepth(isBp1OverDepth || isBp2OverDepth);
  const bool isSkipEvidenceSearch(isMaxDepth && isOverDepth);

  if (!isSkipEvidenceSearch) {
    // count the paired-read fragments supporting the ref and alt alleles in each sample:
    //
    getSVPairSupport(svData, assemblyData, sv, svId, evidence, svSupports);

    // count the split reads supporting the ref and alt alleles in each sample
    //
    getSVSplitReadSupport(assemblyData, sv, svId, baseInfo, evidence, svSupports);

    // fix erroneous pair support based on split evidence:
    resolvePairSplitConflicts(sv, evidence);
  }

  // compute allele likelihoods, and any other summary metric shared between all models:
  //
  getSVSupportSummary(evidence, baseInfo);
}

/// record a set of convenient companion values for any probability
///
struct ProbSet {
  explicit ProbSet(const double initProb)
    : prob(initProb), comp(1 - prob), lnProb(std::log(prob)), lnComp(std::log(comp))
  {
  }

  double prob;
  double comp;
  double lnProb;
  double lnComp;
};

static void incrementSpanningPairAlleleLnLhood(
    const ProbSet&                  selfChimeraProb,
    const ProbSet&                  otherChimeraProb,
    const SVFragmentEvidenceAllele& allele,
    const double                    power,
    double&                         bpLnLhood)
{
  const float fragProb(getSpanningPairAlleleLhood(allele));
  bpLnLhood += std::log(selfChimeraProb.comp * fragProb + otherChimeraProb.prob) * power;
}

static double incrementAlleleSplitReadLhood(
    const ProbSet&                  selfMapProb,
    const ProbSet&                  otherMapProb,
    const SVFragmentEvidenceAllele& allele,
    const double /*readLnPrior*/,
    const std::pair<bool, bool>& isSupported,
    const bool                   isRead1,
    bool&                        isReadEvaluated)
{
  if (!(allele.bp1.getRead(isRead1).isSplitEvaluated && allele.bp2.getRead(isRead1).isSplitEvaluated)) {
    isReadEvaluated = false;
  }

  const double alignBp1LnLhood(allele.bp1.getRead(isRead1).splitLnLhood);
  const double alignBp2LnLhood(allele.bp2.getRead(isRead1).splitLnLhood);

  bool isUseBp1Lhood(isSupported.first);
  if (isSupported.first && isSupported.second) {
    isUseBp1Lhood = (alignBp1LnLhood >= alignBp2LnLhood);
  }

  const double alignLnLhood(isUseBp1Lhood ? alignBp1LnLhood : alignBp2LnLhood);

  const double fragLnLhood =
      log_sum((selfMapProb.lnComp + alignLnLhood), (otherMapProb.lnProb));  //+readLnPrior));

#ifdef DEBUG_SCORE
  static const std::string logtag("incrementAlleleSplitReadLhood: ");
  log_os << logtag  //<< "readPrior: " << readLnPrior
         << " isRead1?: " << isRead1 << "\n";
  log_os << logtag << "isEval " << isReadEvaluated << "\n";
  log_os << logtag << "alignBp1LnLhood " << alignBp1LnLhood << "\n";
  log_os << logtag << "alignBp2LnLhood " << alignBp2LnLhood << "\n";
  log_os << logtag << "selfMap " << selfMapProb.lnProb << "\n";
  log_os << logtag << "otherMap " << otherMapProb.lnProb << "\n";
  log_os << logtag << "increment " << fragLnLhood << "\n";
#endif

  return fragLnLhood;
}

static void incrementSplitReadLhood(
    const std::string& /*fragLabel*/,
    const SVFragmentEvidence& fragev,
    const ProbSet&            refMapProb,
    const ProbSet&            altMapProb,
    const bool                isPermissive,
    const bool                isRead1,
    double&                   refSplitLnLhood,
    double&                   altSplitLnLhood,
    bool&                     isReadEvaluated)
{
  static const double baseLnPrior(std::log(0.25));

#ifdef DEBUG_SCORE
  log_os << __FUNCTION__ << ": pre-support\n";
#endif

  std::pair<bool, bool> isSupported;
  if (isPermissive) {
    isSupported = fragev.isAnyTier2SplitReadSupport(isRead1);
  } else {
    isSupported = fragev.isAnySplitReadSupport(isRead1);
  }

  if (!(isSupported.first || isSupported.second)) {
    isReadEvaluated = false;
    return;
  }

#ifdef DEBUG_SCORE
  log_os << __FUNCTION__ << ": post-support\n";
#endif

  const unsigned readSize(fragev.getRead(isRead1).size);
  const double   readLnPrior(baseLnPrior * readSize);

#ifdef DEBUG_SCORE
  log_os << __FUNCTION__ << ": starting ref\n";
#endif
  const double refSplit = incrementAlleleSplitReadLhood(
      refMapProb, altMapProb, fragev.ref, readLnPrior, isSupported, isRead1, isReadEvaluated);
#ifdef DEBUG_SCORE
  log_os << __FUNCTION__ << ": starting alt\n";
#endif
  const double altSplit = incrementAlleleSplitReadLhood(
      altMapProb, refMapProb, fragev.alt, readLnPrior, isSupported, isRead1, isReadEvaluated);

  // filter out split read evidence with a poor alignment to both alleles:
  /// TODO: fraction of these noise reads could be informative for filtration
  static const float pseudoLnProb(0.5);
  if ((refSplit < (altMapProb.lnProb + pseudoLnProb)) && (altSplit < (refMapProb.lnProb + pseudoLnProb)))
    return;

  refSplitLnLhood += refSplit;
  altSplitLnLhood += altSplit;
}

struct AlleleLnLhood {
  double fragPair   = 0.;
  double read1Split = 0.;
  double read2Split = 0.;
};

static double getFragLnLhood(
    const AlleleLnLhood& al, const bool isRead1Evaluated, const bool isRead2Evaluated)
{
#ifdef DEBUG_SCORE
  log_os << "getFragLnLhood: frag/read1/read2 " << al.fragPair << " " << al.read1Split << " " << al.read2Split
         << "\n";
  log_os << "getFragLnLhood: isread1/isread2 " << isRead1Evaluated << " " << isRead2Evaluated << "\n";
#endif

  double ret(al.fragPair);

  // limit split read evidence to only one read, because it's only possible for one section
  // of the molecule to independently cross the breakend:
  if (isRead1Evaluated) {
    if (isRead2Evaluated) {
      ret += std::max(al.read1Split, al.read2Split);
    } else {
      ret += al.read1Split;
    }
  } else if (isRead2Evaluated) {
    ret += al.read2Split;
  }

  return ret;
}

/// Combine the paired and split read likelihoods of one fragment, for the ref and the alt allele
///
/// \param[in] isUseBestAltRead If false, the best split read likelihood of the fragment is chosen separately
/// for each allele, see getFragLnLhood(). If true and both reads have been evaluated, the split read
/// likelihoods of both alleles are taken from the same read: the read which supports the alt allele most
/// strongly, relative to the reference allele.
///
/// The second case is designed for a tandem duplication. A molecule from a haplotype carrying the duplication
/// can cross both the duplication junction and one of the reference breakends, which are still intact on that
/// haplotype. For a duplication smaller than the fragment size it is common for one read to cross each. The
/// read at the reference breakend then matches the reference allele at least as well as the other read matches
/// the alt allele, so that choosing the best read per allele makes the fragment uninformative at best, and
/// evidence for the reference allele whenever the read on the duplication junction has more basecall errors.
/// In terms of where the molecule came from, the likelihood ratio between 'covers the duplication junction'
/// and 'covers reference sequence only' is the largest such ratio of any one read, because the other read is
/// reference sequence under both hypotheses.
///
static void getRefAltFragLnLhood(
    const AlleleLnLhood& refLnLhoodSet,
    const AlleleLnLhood& altLnLhoodSet,
    bool                 isRead1Evaluated,
    bool                 isRead2Evaluated,
    const bool           isUseBestAltRead,
    double&              refLnFragLhood,
    double&              altLnFragLhood)
{
  if (isUseBestAltRead && isRead1Evaluated && isRead2Evaluated) {
    const double read1AltSupport(altLnLhoodSet.read1Split - refLnLhoodSet.read1Split);
    const double read2AltSupport(altLnLhoodSet.read2Split - refLnLhoodSet.read2Split);
    if (read1AltSupport >= read2AltSupport) {
      isRead2Evaluated = false;
    } else {
      isRead1Evaluated = false;
    }
  }
  refLnFragLhood = getFragLnLhood(refLnLhoodSet, isRead1Evaluated, isRead2Evaluated);
  altLnFragLhood = getFragLnLhood(altLnLhoodSet, isRead1Evaluated, isRead2Evaluated);
}

/// when an sv is treated as 'small', we skip all paired-read evidence and rely on split reads only:
///
/// with further model improvements we can add pairs back into the small variant calls:
///
/// this function returns 1 for a variant which is "fully large" and 0 for a variant which is "fully small",
/// with intermediate values for sizes in between
///
/// \param[in] tandemDupFullPairSize The tandem duplication size from which paired-read evidence is given full
/// weight. The weight is zero up to half of this size. Read pairs supporting a tandem duplication are only
/// found through their anomalous orientation, so this should be set to the size of a large fragment in the
/// sequenced library. If zero, paired-read evidence is given full weight for tandem duplications of any size.
///
static float getSpanningPairWeight(const SVCandidate& sv, const int tandemDupFullPairSize = 0)
{
  const auto svType(getExtendedSVType(sv));
  if (svType == EXTENDED_SV_TYPE::TANDUP) {
    if (tandemDupFullPairSize <= 1) return 1.f;
    const LinearScaler<int> dupSizeRamp(tandemDupFullPairSize / 2, tandemDupFullPairSize);
    return dupSizeRamp.getScale(sv.centerSize());
  }
  if (!((svType == EXTENDED_SV_TYPE::INSERT) || (svType == EXTENDED_SV_TYPE::DELETE))) return 1.f;

  if ((svType == EXTENDED_SV_TYPE::INSERT) && SVScorePairAltProcessor::isLargeInsertSV(sv)) {
    static const int               minInsertSmallSize(100);
    static const int               maxInsertSmallSize(150);
    static const LinearScaler<int> insertSizeRamp(minInsertSmallSize, maxInsertSmallSize);

    return insertSizeRamp.getScale(sv.insertSeq.size());
  } else {
    /// TODO set these numbers from insert size:
    static const int               minSmallSize(300);
    static const int               maxSmallSize(500);
    static const LinearScaler<int> svSizeRamp(minSmallSize, maxSmallSize);

    return svSizeRamp.getScale(sv.centerSize());
  }
}
/// return true if any evidence exists for fragment:
///
/// \param semiMappedPower multiply out semi-mapped reads (in log space) by this value
///
static bool getRefAltFromFrag(
    const float               spanningPairWeight,
    const double              semiMappedPower,
    const ProbSet&            refChimeraProb,
    const ProbSet&            altChimeraProb,
    const ProbSet&            refSplitMapProb,
    const ProbSet&            altSplitMapProb,
    const bool                isPermissive,
    const std::string&        fragLabel,
    const SVFragmentEvidence& fragev,
    AlleleLnLhood&            refLnLhoodSet,
    AlleleLnLhood&            altLnLhoodSet,
    bool&                     isRead1Evaluated,
    bool&                     isRead2Evaluated)
{
#ifdef DEBUG_SCORE
  log_os << __FUNCTION__ << ": qname: " << fragLabel << " fragev: " << fragev << "\n";
#endif

  /// TODO: add read pairs with one shadow read to the alt read pool

  bool isFragEvaluated(false);

  // high-quality spanning support relies on read1 and read2 mapping well:
  bool isPairUsable;
  if (isPermissive) {
    isPairUsable =
        (fragev.read1.isObservedAnchor(isPermissive) || fragev.read2.isObservedAnchor(isPermissive));
  } else {
    isPairUsable =
        ((fragev.read1.isScanned && fragev.read2.isScanned) &&
         (fragev.read1.isAnchored(isPermissive) || fragev.read2.isAnchored(isPermissive)));
  }

  if (isPairUsable) {
    /// only add to the likelihood if the fragment supports at least one allele:
    if (fragev.isAnySpanningPairSupport()) {
      // reduce the impact of spanning reads to zero as svs become small, this is because of complex
      // signal/noise which the scoring models haven't (yet) been designed to handle.
      const bool isSemiMapped(
          !(fragev.read1.isAnchored(isPermissive) && fragev.read2.isAnchored(isPermissive)));
      double spanPower(spanningPairWeight);

      if (isSemiMapped) {
        // only count semi-mapped reads for the alt allele
        if (getSpanningPairAlleleLhood(fragev.alt) > getSpanningPairAlleleLhood(fragev.ref)) {
          spanPower *= semiMappedPower;
        } else {
          spanPower = 0.;
        }
      }

      incrementSpanningPairAlleleLnLhood(
          refChimeraProb, altChimeraProb, fragev.ref, spanPower, refLnLhoodSet.fragPair);
      incrementSpanningPairAlleleLnLhood(
          altChimeraProb, refChimeraProb, fragev.alt, spanPower, altLnLhoodSet.fragPair);
      isFragEvaluated = true;
    }
  }

  /// split support is less dependent on mapping quality of the individual read, because
  /// we're potentially relying on shadow reads recovered from the unmapped state
  isRead1Evaluated = true;
  isRead2Evaluated = true;
#ifdef DEBUG_SCORE
  log_os << __FUNCTION__ << ": starting read1 split\n";
#endif
  incrementSplitReadLhood(
      fragLabel,
      fragev,
      refSplitMapProb,
      altSplitMapProb,
      isPermissive,
      true,
      refLnLhoodSet.read1Split,
      altLnLhoodSet.read1Split,
      isRead1Evaluated);
#ifdef DEBUG_SCORE
  log_os << __FUNCTION__ << ": starting read2 split\n";
#endif
  incrementSplitReadLhood(
      fragLabel,
      fragev,
      refSplitMapProb,
      altSplitMapProb,
      isPermissive,
      false,
      refLnLhoodSet.read2Split,
      altLnLhoodSet.read2Split,
      isRead2Evaluated);

#ifdef DEBUG_SCORE
  log_os << __FUNCTION__ << ": iseval frag/read1/read2: " << isFragEvaluated << " " << isRead1Evaluated << " "
         << isRead1Evaluated << "\n";
#endif
  return (isFragEvaluated || isRead1Evaluated || isRead2Evaluated);
}

/// Expected fraction of the fragments informative for a variant which support its alt allele, in log space,
/// for each diploid genotype
struct DiploidAlleleFractions {
  std::array<double, DIPLOID_GT::SIZE> altLnFraction;
  std::array<double, DIPLOID_GT::SIZE> altLnCompFraction;
};

/// Allele fractions used for all variant types other than tandem duplications
static DiploidAlleleFractions getDefaultDiploidAlleleFractions()
{
  DiploidAlleleFractions fractions;
  for (unsigned gt(0); gt < DIPLOID_GT::SIZE; ++gt) {
    const DIPLOID_GT::index_t gtid(static_cast<DIPLOID_GT::index_t>(gt));
    fractions.altLnFraction[gt]     = DIPLOID_GT::altLnFraction(gtid);
    fractions.altLnCompFraction[gt] = DIPLOID_GT::altLnCompFraction(gtid);
  }
  return fractions;
}

/// Allele fractions expected for a tandem duplication
///
/// For most variant types the reference allele junction is destroyed on a haplotype carrying the variant, so
/// that reads supporting the reference allele can only come from a reference haplotype. A tandem duplication is
/// different: a haplotype carrying it still contains the unchanged reference sequence across the start of the
/// first copy and across the end of the second copy, in addition to the one new junction between the copies.
/// Reads across those two positions are indistinguishable from reference reads and are (correctly) counted as
/// reference allele support.
///
/// Write d for the duplication size and w for the number of positions at which a read can be placed such that
/// it crosses one given junction with enough sequence on both sides to be counted as evidence. The observations
/// which are informative for the duplication are the reads crossing at least one junction:
///
/// - On a duplication haplotype A-D-D-B, w placements cross the new junction between the two copies and
///   support the alt allele. The placements which cross A-D or D-B without also crossing the new junction
///   support the reference allele, there are min(w, d) of these on each side.
/// - On a reference haplotype A-D-B, w placements cross the first reference breakend and a further min(w, d)
///   cross the second breakend only. All of these support the reference allele.
///
/// With r = min(1, d/w), the expected fraction of alt-supporting observations is therefore:
///
///   heterozygous: 1 / (2 + 3r)
///   homozygous:   1 / (1 + 2r)
///
/// So for any duplication larger than a read, 1 in 5 informative observations support a heterozygous
/// duplication and 1 in 3 support a homozygous duplication, compared to 1 in 2 and (nearly) all for other
/// variant types. Genotyping a duplication with the fractions used for other variant types can therefore
/// never yield a homozygous genotype. As the duplication size goes to zero, the expected fractions converge
/// on those of a simple insertion, which is what a very small tandem duplication is.
///
/// These fractions are used to choose between the heterozygous and homozygous genotype only, see
/// applyTandemDupZygosity().
///
/// The fractions are scaled such that they are equal to the default values of all other variant types at r = 0.
///
/// \param[in] dupSize duplication size (d)
///
/// \param[in] junctionWindow number of read placements crossing one junction (w)
///
static DiploidAlleleFractions getTandemDupAlleleFractions(const double dupSize, const double junctionWindow)
{
  assert(junctionWindow > 0);

  const double r(std::min(1., std::max(0., dupSize) / junctionWindow));

  std::array<double, DIPLOID_GT::SIZE> altFraction;
  altFraction[DIPLOID_GT::REF] = 0.;
  altFraction[DIPLOID_GT::HET] = DIPLOID_GT::altFraction(DIPLOID_GT::HET) / (1. + 1.5 * r);
  altFraction[DIPLOID_GT::HOM] = DIPLOID_GT::altFraction(DIPLOID_GT::HOM) / (1. + 2. * r);

  DiploidAlleleFractions fractions;
  for (unsigned gt(0); gt < DIPLOID_GT::SIZE; ++gt) {
    fractions.altLnFraction[gt]     = std::log(altFraction[gt]);
    fractions.altLnCompFraction[gt] = std::log(1. - altFraction[gt]);
  }
  return fractions;
}

/// Replace the heterozygous and homozygous genotype likelihoods of a tandem duplication with values which
/// split the non-reference probability according to the tandem duplication allele fraction model
///
/// The question of whether the sample carries the duplication at all is left to the standard model used for
/// every variant type, so that neither the set of duplications called nor their QUAL depends on the
/// duplication model. This is the conservative choice: a heterozygous duplication is supported by only 1 in 5
/// informative reads, so a model expecting this fraction also finds a junction supported by a few percent of
/// reads far more plausible as a variant than the standard model does, and such low fraction junctions are
/// more likely to be artifacts than germline variants. The duplication model is used only for the question the
/// standard model cannot answer, which is how many copies of the duplication are present.
///
/// \param[in] logPrior Genotype prior
///
/// \param[in] dupLoglhood Genotype likelihoods from the tandem duplication allele fractions
///
/// \param[in,out] loglhood Genotype likelihoods from the standard allele fractions. On return the total
/// posterior probability of the two non-reference genotypes is unchanged, and the ratio between them is that
/// found from \p dupLoglhood
///
template <typename PriorType>
static void applyTandemDupZygosity(
    const PriorType&                            logPrior,
    const std::array<double, DIPLOID_GT::SIZE>& dupLoglhood,
    std::array<double, DIPLOID_GT::SIZE>&       loglhood)
{
  using namespace DIPLOID_GT;

  const double nonRefLnPost(log_sum(loglhood[HET] + logPrior[HET], loglhood[HOM] + logPrior[HOM]));
  const double dupNonRefLnPost(log_sum(dupLoglhood[HET] + logPrior[HET], dupLoglhood[HOM] + logPrior[HOM]));

  // nothing to divide up if the duplication model has no finite likelihood for either genotype
  if (!(std::isfinite(nonRefLnPost) && std::isfinite(dupNonRefLnPost))) return;

  for (const index_t gt : {HET, HOM}) {
    loglhood[gt] = nonRefLnPost + (dupLoglhood[gt] - dupNonRefLnPost);
  }
}

/// \return The mean size of all reads scanned as evidence in one sample, or zero if there are none
static double getMeanEvidenceReadSize(const SVEvidence::evidenceTrack_t& sampleEvidence)
{
  double   sizeSum(0);
  unsigned readCount(0);
  for (const SVEvidence::evidenceTrack_t::value_type& val : sampleEvidence) {
    for (const bool isRead1 : {true, false}) {
      const SVFragmentEvidenceRead& evRead(val.second.getRead(isRead1));
      if ((!evRead.isScanned) || (evRead.size == 0)) continue;
      sizeSum += evRead.size;
      readCount++;
    }
  }
  if (readCount == 0) return 0;
  return (sizeSum / readCount);
}

/// \return true if the sample has a confident non-reference genotype with uncertain zygosity
///
/// This is the case for a sample which fails the minimum genotype quality filter and no other sample filter,
/// while the probability of the homozygous reference genotype is as low as the genotype quality filter would
/// require of the probability of all genotypes other than the called one.
static bool isConfidentNonRefSample(
    const CallOptionsDiploid& diploidOpt, const SVScoreInfoDiploidSample& diploidSampleInfo)
{
  if (diploidSampleInfo.gt == DIPLOID_GT::REF) return false;
  if (diploidSampleInfo.filters.size() != 1) return false;
  if (diploidSampleInfo.filters.count(diploidOpt.minGTFilterLabel) != 1) return false;

  const int nonRefScore(error_prob_to_qphred(diploidSampleInfo.pprob[DIPLOID_GT::REF]));
  return (nonRefScore >= static_cast<int>(diploidOpt.minPassGTScore));
}

/// score diploid germline specific components:
///
/// \param[in] isUseBestAltRead See getRefAltFragLnLhood()
///
static void addDiploidLoglhood(
    const float                           spanningPairWeight,
    const DiploidAlleleFractions&         alleleFractions,
    const bool                            isUseBestAltRead,
    const SVEvidence::evidenceTrack_t&    sampleEvidence,
    std::array<double, DIPLOID_GT::SIZE>& loglhood)
{
  for (const SVEvidence::evidenceTrack_t::value_type& val : sampleEvidence) {
    const std::string&        fragLabel(val.first);
    const SVFragmentEvidence& fragev(val.second);

    AlleleLnLhood refLnLhoodSet, altLnLhoodSet;
    bool          isRead1Evaluated(true);
    bool          isRead2Evaluated(true);

    /// TODO: set this value from error rates observed in input data:
    //
    // put some more thought into this -- is this P (spurious | any old read) or P( spurious | chimera ) ??
    // it seems like it should be the latter in the usages that really matter.
    //
    static const ProbSet chimeraProb(1e-3);

    // use a constant mapping prob for now just to get the zero-th order concept into the model
    // that "reads are mismapped at a non-trivial rate"
    /// TODO: experiment with per-read mapq values
    static const ProbSet refSplitMapProb(1e-6);
    static const ProbSet altSplitMapProb(1e-5);

    // don't use semi-mapped reads for germline calling:
    static const double semiMappedPower(0.);

    static const bool isPermissive(false);

    if (!getRefAltFromFrag(
            spanningPairWeight,
            semiMappedPower,
            chimeraProb,
            chimeraProb,
            refSplitMapProb,
            altSplitMapProb,
            isPermissive,
            fragLabel,
            fragev,
            refLnLhoodSet,
            altLnLhoodSet,
            isRead1Evaluated,
            isRead2Evaluated)) {
      // continue if this fragment was not evaluated for pair or split support for either allele:
      continue;
    }

    double refLnFragLhood(0);
    double altLnFragLhood(0);
    getRefAltFragLnLhood(
        refLnLhoodSet,
        altLnLhoodSet,
        isRead1Evaluated,
        isRead2Evaluated,
        isUseBestAltRead,
        refLnFragLhood,
        altLnFragLhood);

#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": refLnFragLhood: " << refLnFragLhood << "\n";
    log_os << __FUNCTION__ << ": altLnFragLhood: " << altLnFragLhood << "\n";
#endif

    for (unsigned gt(0); gt < DIPLOID_GT::SIZE; ++gt) {
      using namespace DIPLOID_GT;

#ifdef DEBUG_SCORE
      log_os << __FUNCTION__ << ": starting gt: " << gt << " " << label(gt) << "\n";
#endif

      const double refLnLhood(refLnFragLhood + alleleFractions.altLnCompFraction[gt]);
      const double altLnLhood(altLnFragLhood + alleleFractions.altLnFraction[gt]);
      loglhood[gt] += log_sum(refLnLhood, altLnLhood);

#ifdef DEBUG_SCORE
      log_os << __FUNCTION__ << ": gt/fragref/ref/fragalt/alt/loglhood: " << label(gt) << " "
             << refLnFragLhood << " " << refLnLhood << " " << altLnFragLhood << " " << altLnLhood << " "
             << loglhood[gt] << "\n";
#endif
    }
  }
}

/// score diploid germline specific components:
static void scoreDiploidSV(
    const CallOptionsDiploid&            diploidOpt,
    const SVLocusScanner&                readScanner,
    const CallOptionsDiploidDeriv&       diploidDopt,
    const ChromDepthFilterUtil&          dFilter,
    const std::vector<JunctionCallInfo>& junctionData,
    SVScoreInfoDiploid&                  diploidInfo)
{
  //
  // compute qualities
  //
  static const int maxQ(999);

  assert(!junctionData.empty());

  double jointRefProb(1.);

  // The tandem duplication model only applies to a duplication called as a single junction, not to a
  // duplication-like junction of a larger event.
  const bool isTandemDupModel(
      diploidOpt.isTandemDupModel && (junctionData.size() == 1) &&
      (getSVType(junctionData[0].getSV()) == SV_TYPE::TANDUP));

  static const DiploidAlleleFractions defaultAlleleFractions(getDefaultDiploidAlleleFractions());

  const unsigned diploidSampleCount(diploidInfo.samples.size());
  for (unsigned diploidSampleIndex(0); diploidSampleIndex < diploidSampleCount; ++diploidSampleIndex) {
    SVScoreInfoDiploidSample& diploidSampleInfo(diploidInfo.samples[diploidSampleIndex]);

    std::array<double, DIPLOID_GT::SIZE> loglhood;
    std::fill(loglhood.begin(), loglhood.end(), 0);
    for (unsigned junctionIndex(0); junctionIndex < junctionData.size(); ++junctionIndex) {
      const JunctionCallInfo&            junction(junctionData[junctionIndex]);
      const SVEvidence::evidenceTrack_t& etrack(junction.getEvidence().samples[diploidSampleIndex]);

      if (!isTandemDupModel) {
        addDiploidLoglhood(junction.getSpanningWeight(), defaultAlleleFractions, false, etrack, loglhood);
        continue;
      }

      // paired-read support for a duplication is found from anomalous read pair orientation only, which is
      // not reliably observed until the duplication is as large as a large fragment of the library
      const float spanningPairWeight(getSpanningPairWeight(
          junction.getSV(), static_cast<int>(readScanner.getExtremeFifthRange().max)));

      addDiploidLoglhood(spanningPairWeight, defaultAlleleFractions, true, etrack, loglhood);

      // the number of placements at which a read crosses a junction with the minimum flank size on both
      // sides. Without any reads there is no evidence to score.
      const double meanReadSize(getMeanEvidenceReadSize(etrack));
      const double junctionWindow(meanReadSize - (2. * getMinSplitReadFlankSize()));
      if (junctionWindow >= 1.) {
        const DiploidAlleleFractions dupAlleleFractions(
            getTandemDupAlleleFractions(junction.getSV().centerSize(), junctionWindow));

        std::array<double, DIPLOID_GT::SIZE> dupLoglhood;
        std::fill(dupLoglhood.begin(), dupLoglhood.end(), 0);
        addDiploidLoglhood(spanningPairWeight, dupAlleleFractions, true, etrack, dupLoglhood);
        applyTandemDupZygosity(diploidDopt.logPrior, dupLoglhood, loglhood);
      }
    }
    std::array<double, DIPLOID_GT::SIZE> pprob;
    for (unsigned gt(0); gt < DIPLOID_GT::SIZE; ++gt) {
      pprob[gt] = loglhood[gt] + diploidDopt.logPrior[gt];
    }

    unsigned maxGt(0);
    normalizeLogDistro(pprob.begin(), pprob.end(), maxGt);

#ifdef DEBUG_SCORE
    for (unsigned gt(0); gt < DIPLOID_GT::SIZE; ++gt) {
      log_os << __FUNCTION__ << ": gt/lhood/prior/pprob: " << DIPLOID_GT::label(gt) << " " << loglhood[gt]
             << " " << diploidDopt.prior[gt] << " " << pprob[gt] << "\n";
    }
#endif

    diploidSampleInfo.gt = static_cast<DIPLOID_GT::index_t>(maxGt);
    diploidSampleInfo.gtScore =
        std::min(maxQ, error_prob_to_qphred(prob_comp(pprob.begin(), pprob.end(), diploidSampleInfo.gt)));

    // set phredLoghood:
    {
      unsigned maxIndex(0);
      for (unsigned gt(1); gt < DIPLOID_GT::SIZE; ++gt) {
        if (loglhood[gt] > loglhood[maxIndex]) maxIndex = gt;
      }
      for (unsigned gt(0); gt < DIPLOID_GT::SIZE; ++gt) {
        diploidSampleInfo.pprob[gt] = pprob[gt];
        diploidSampleInfo.phredLoghood[gt] =
            std::min(maxQ, ln_error_prob_to_qphred(loglhood[gt] - loglhood[maxIndex]));
      }
    }

    jointRefProb *= pprob[DIPLOID_GT::REF];
  }
  diploidInfo.altScore = std::min(maxQ, error_prob_to_qphred(jointRefProb));

  //
  // apply filters
  //
  {
    if (diploidInfo.altScore < diploidOpt.minPassAltScore) {
      diploidInfo.filters.insert(diploidOpt.minAltFilterLabel);
    }

    // add sample specific filters
    bool isAllSampleFiltered(true);
    for (unsigned sampleIndex(0); sampleIndex < diploidSampleCount; ++sampleIndex) {
      SVScoreInfoDiploidSample& diploidSampleInfo(diploidInfo.samples[sampleIndex]);
      if (diploidSampleInfo.gt == DIPLOID_GT::REF) diploidSampleInfo.filters.insert(diploidOpt.homRefLabel);

      if (diploidSampleInfo.gtScore < diploidOpt.minPassGTScore)
        diploidSampleInfo.filters.insert(diploidOpt.minGTFilterLabel);

      if (diploidSampleInfo.filters.empty()) {
        isAllSampleFiltered = false;
      } else if ((!diploidOpt.isFilterUncertainZygosity) && isConfidentNonRefSample(diploidOpt, diploidSampleInfo)) {
        // The sample keeps its genotype quality filter, but the variant is not filtered at the record level
        // because of it: the uncertainty is in the number of copies of the variant, not in its presence.
        isAllSampleFiltered = false;
      }
    }

    // If no sample passes all sample-specific filters, apply sample FT filter at the record level
    if (isAllSampleFiltered) diploidInfo.filters.insert(diploidOpt.failedSampleFTLabel);

    const unsigned junctionCount(junctionData.size());

    // apply high depth filter:
    if (dFilter.isMaxDepthFilter()) {
      unsigned filteredJunctionCount(0);
      for (const JunctionCallInfo& junction : junctionData) {
        const SVScoreInfo& baseInfo(junction.getBaseInfo());
        const SVCandidate& sv(junction.getSV());

        // apply maxdepth filter if either of the breakpoints exceeds the maximum depth:
        if ((baseInfo.bp1MaxDepth > dFilter.maxDepth(sv.bp1.interval.tid)) ||
            (baseInfo.bp2MaxDepth > dFilter.maxDepth(sv.bp2.interval.tid))) {
          filteredJunctionCount++;
        }
      }

      if ((filteredJunctionCount * 2) > junctionCount) {
        diploidInfo.filters.insert(diploidOpt.maxDepthFilterLabel);
      }
    }

    // apply MQ0 filter
    {
      unsigned filteredJunctionCount(0);
      for (const JunctionCallInfo& junction : junctionData) {
        const SVScoreInfo& baseInfo(junction.getBaseInfo());
        const SVCandidate& sv(junction.getSV());

        const bool isMQ0FilterSize(isSVBelowMinSize(sv, 1000));
        if (isMQ0FilterSize) {
          // apply MQ0 filter for one junction if either breakend meets the filter criteria:
          if ((baseInfo.bp1MQ0Frac > diploidOpt.maxMQ0Frac) ||
              (baseInfo.bp2MQ0Frac > diploidOpt.maxMQ0Frac)) {
            filteredJunctionCount++;
          }
        }
      }

      if ((filteredJunctionCount * 2) > junctionCount) {
        diploidInfo.filters.insert(diploidOpt.maxMQ0FracLabel);
      }
    }

    // apply zero pair filter
    {
      // this size represents the outer edge of variant size above which we expect pair
      // discovery to suffer no dropouts due to normal pair distro sizes
      static const double insertSizeFactor(1);
      const unsigned      maxClosePairSize(readScanner.getExtremeFifthRange().max * insertSizeFactor);

      unsigned filteredJunctionCount(0);
      for (const JunctionCallInfo& junction : junctionData) {
        const SVScoreInfo& baseInfo(junction.getBaseInfo());
        const SVCandidate& sv(junction.getSV());

        // only apply the zero pair filter to variants that should definitely have located supporting pairs:
        const EXTENDED_SV_TYPE::index_t svType(getExtendedSVType(sv));
        const bool                      isZeroPairFilterSize(
            (svType != EXTENDED_SV_TYPE::INSERT) && (!isSVBelowMinSize(sv, maxClosePairSize)));

        if (isZeroPairFilterSize) {
          unsigned totalDiploidSpanningPairCount(0);
          for (unsigned diploidSampleIndex(0); diploidSampleIndex < diploidSampleCount;
               ++diploidSampleIndex) {
            totalDiploidSpanningPairCount +=
                baseInfo.samples[diploidSampleIndex].alt.confidentSpanningPairCount;
          }

          if (totalDiploidSpanningPairCount == 0) {
            filteredJunctionCount++;
          }
        }
      }

      if ((filteredJunctionCount * 2) > junctionCount) {
        diploidInfo.filters.insert(diploidOpt.noPairSupportLabel);
      }
    }
  }
}

void SVScorer::computeAllScoreModels(
    const std::vector<JunctionCallInfo>& junctionData, SVModelScoreInfo& modelScoreInfo)
{
  scoreDiploidSV(
      _diploidOpt, _readScanner, _diploidDopt, _dFilterDiploid, junctionData, modelScoreInfo.diploid);
}

void SVScorer::scoreSV(
    const SVCandidateSetData&                   svData,
    const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
    const SVMultiJunctionCandidate&             mjSV,
    const std::vector<SVId>&                    mjSVId,
    const std::vector<bool>&                    isJunctionFiltered,
    std::vector<SVModelScoreInfo>&              mjModelScoreInfo,
    SVModelScoreInfo&                           mjJointModelScoreInfo,
    bool&                                       isMJEvent,
    SVEvidenceWriterData&                       svEvidenceWriterData)
{
  // scoring is roughly divided into two parts -- treating individual dna-junctions
  // independently (the simpler call mechanism used the great majority of the time) and
  // joint junction analysis for larger scale events
  //
  const unsigned junctionCount(mjSV.junction.size());
  mjModelScoreInfo.resize(junctionCount);
  std::vector<SVEvidence> junctionEvidence(junctionCount);

  // set the degree to which we rely on spanning pair evidence for each junction, see getSpanningPairWeight()
  std::vector<float> junctionSpanningPairWeight(junctionCount);

  for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
    mjModelScoreInfo[junctionIndex].setSampleCount(_sampleCount);
    junctionEvidence[junctionIndex].samples.resize(_sampleCount);
  }

  mjJointModelScoreInfo.setSampleCount(_sampleCount);
  mjJointModelScoreInfo.clear();

  unsigned unfilteredJunctionCount(0);

  std::vector<JunctionCallInfo> junctionData;

  for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
    if (isJunctionFiltered[junctionIndex]) continue;

#ifdef ANY_DEBUG_SCORE
    log_os << __FUNCTION__ << ": Scoring single junction " << junctionIndex << "/" << junctionCount << "\n";
#endif

    unfilteredJunctionCount++;

    const SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
    const SVCandidate&             sv(mjSV.junction[junctionIndex]);
    const SVId&                    svId(mjSVId[junctionIndex]);
    SVModelScoreInfo&              modelScoreInfo(mjModelScoreInfo[junctionIndex]);

    modelScoreInfo.clear();

    try {
      // accumulate model-agnostic evidence for each candidate (or its corresponding reference allele)
      SVEvidence& evidence(junctionEvidence[junctionIndex]);
      getSVSupportingEvidence(
          svData, assemblyData, sv, svId, modelScoreInfo.base, evidence, svEvidenceWriterData);

      if (_diploidOpt.isReportDepthRatio) {
        getSVDepthSupport(sv, modelScoreInfo.base);
      }

      // score components specific to diploid-germline model:
      float& spanningPairWeight(junctionSpanningPairWeight[junctionIndex]);
      ;
      spanningPairWeight = (getSpanningPairWeight(sv));

      junctionData.resize(1);
      junctionData[0].init(sv, evidence, modelScoreInfo.base, spanningPairWeight);

      computeAllScoreModels(junctionData, modelScoreInfo);
    } catch (illumina::common::ExceptionData& e) {
      std::ostringstream oss;
      oss << "Exception caught while attempting to score " << sv;
      e << boost::error_info<struct scoring_candidate_info, std::string>(oss.str());
      throw;
    }

    catch (...) {
      log_os << "Exception caught while attempting to score " << sv << "\n";
      throw;
    }
  }

  //
  // handle multi-junction case:
  //
  if (unfilteredJunctionCount == 1) {
    isMJEvent = false;
  } else if (unfilteredJunctionCount == 2) {
#ifdef ANY_DEBUG_SCORE
    log_os << __FUNCTION__ << ": Scoring multi-junction " << junctionCount << "\n";
#endif
    isMJEvent = true;

    junctionData.resize(unfilteredJunctionCount);
    for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
      if (isJunctionFiltered[junctionIndex]) continue;

      junctionData[junctionIndex].init(
          mjSV.junction[junctionIndex],
          junctionEvidence[junctionIndex],
          mjModelScoreInfo[junctionIndex].base,
          junctionSpanningPairWeight[junctionIndex]);
    }

    computeAllScoreModels(junctionData, mjJointModelScoreInfo);
  } else {
    using namespace illumina::common;
    std::ostringstream oss;
    oss << "Unexpected junction count: " << unfilteredJunctionCount;
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
  }
}
