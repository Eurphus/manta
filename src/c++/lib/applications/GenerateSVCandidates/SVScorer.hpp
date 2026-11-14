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

#pragma once

#include <string>
#include <vector>

#include "GSCOptions.hpp"
#include "JunctionCallInfo.hpp"
#include "SVEvidence.hpp"
#include "SVScorePairProcessor.hpp"
#include "SplitReadAlignment.hpp"
#include "assembly/AssembledContig.hpp"
#include "blt_util/qscore_snp.hpp"
#include "htsapi/bam_header_info.hpp"
#include "htsapi/bam_streamer.hpp"
#include "manta/ChromDepthFilterUtil.hpp"
#include "manta/SVCandidateAssemblyData.hpp"
#include "manta/SVCandidateSetData.hpp"
#include "manta/SVLocusScanner.hpp"
#include "manta/SVModelScoreInfo.hpp"
#include "manta/SVMultiJunctionCandidate.hpp"

struct CallOptionsSharedDeriv {
  explicit CallOptionsSharedDeriv(const CallOptionsShared& opt) : refQ(opt.snpPrior), altQ(0) {}

  const qscore_snp refQ;
  const qscore_snp altQ;
};

struct CallOptionsDiploidDeriv : private boost::noncopyable {
  explicit CallOptionsDiploidDeriv(const CallOptionsDiploid& opt)
  {
    using namespace DIPLOID_GT;

    assert(opt.indelPrior < 0.5);

    prior[HET] = opt.indelPrior;
    prior[HOM] = opt.indelPrior / 2;
    prior[REF] = 1 - prior[HET] - prior[HOM];

    for (unsigned i(0); i < SIZE; ++i) {
      logPrior[i] = std::log(prior[i]);
    }
  }

  std::array<float, DIPLOID_GT::SIZE> prior;
  std::array<float, DIPLOID_GT::SIZE> logPrior;
};

/// Implements SV scoring/genotyping process
///
struct SVScorer {
  SVScorer(const GSCOptions& opt, const SVLocusScanner& readScanner, const bam_header_info& header);

  /// Gather supporting read evidence and generate:
  /// diploid quality score and genotype for SV candidate
  ///
  /// \param mjModelScoreInfo standard (junction-independent) scoring information for each junction
  ///
  /// \param mjJointModelScoreInfo  joint junction scoring information used for cases where a multi-junction
  /// event is detected
  void scoreSV(
      const SVCandidateSetData&                   svData,
      const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
      const SVMultiJunctionCandidate&             mjSV,
      const std::vector<SVId>&                    mjSVId,
      const std::vector<bool>&                    isJunctionFiltered,
      std::vector<SVModelScoreInfo>&              mjModelScoreInfo,
      SVModelScoreInfo&                           mjJointModelScoreInfo,
      bool&                                       isMJEvent,
      SVEvidenceWriterData&                       svEvidenceWriterData);

  typedef std::shared_ptr<SVScorePairProcessor> pairProcPtr;
  typedef std::shared_ptr<bam_streamer>         streamPtr;

#if 0
    unsigned
    sampleCount() const
    {
        return _sampleCount;
    }

    const std::vector<std::string>&
    sampleNames() const
    {
        return _sampleNames;
    }
#endif

  /// TestSVScorer is a friend structure of SV scorer. So that it can access private
  /// methods of SVScorer. As SVScorer has many private methods, for unit test writing
  /// this friend structure has been created.
  friend struct TestSVScorer;

private:
  void processExistingAltPairInfo(
      const PairOptions&        pairOpt,
      const SVCandidateSetData& svData,
      const SVCandidate&        sv,
      const SVId&               svId,
      SVEvidence&               evidence,
      SVEvidenceWriterData&     svSupports);

  /// Estimate pair support for an sv candidate
  ///
  /// This is restricted to simple indel-style SVs
  void getSVAltPairSupport(
      const PairOptions&             pairOpt,
      const SVCandidateAssemblyData& assemblyData,
      const SVCandidate&             sv,
      SVEvidence&                    evidence,
      std::vector<pairProcPtr>&      pairProcList);

  /// Find spanning read support for the reference allele for sv candidate
  void getSVRefPairSupport(
      const PairOptions&        pairOpt,
      const SVCandidate&        sv,
      SVEvidence&               evidence,
      std::vector<pairProcPtr>& pairProcList);

  /// Find paired read support for ref and alt alleles
  void getSVPairSupport(
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& assemblyData,
      const SVCandidate&             sv,
      const SVId&                    svId,
      SVEvidence&                    evidence,
      SVEvidenceWriterData&          svSupports);

  /// \brief Find split read support for ref and alt alleles
  ///
  /// \param assemblyData Assembly results for the given candidate SV
  ///
  /// \param sv Candidate SV for which split read support is being enumerated
  ///
  /// \param svId Structural variant ID tag used to annotate BAM output used for debugging
  ///
  /// \param ssInfo Model-agnostic summary information for the candidate SV
  ///
  /// \param evidence The structure used to return all split read information for the variant, used for
  /// generating scores under specific variant models.
  ///
  /// \param svSupports Structure into which supporting evidence reads can be added for BAM output used in
  /// debugging
  void getSVSplitReadSupport(
      const SVCandidateAssemblyData& assemblyData,
      const SVCandidate&             sv,
      const SVId&                    svId,
      SVScoreInfo&                   ssInfo,
      SVEvidence&                    evidence,
      SVEvidenceWriterData&          svSupports);

  /// Determine maximum depth and MQ0 frac in region around breakend, over all samples
  void getBreakendMaxMappedDepthAndMQ0(
      const bool        isMaxDepth,
      const double      cutoffDepth,
      const SVBreakend& bp,
      unsigned&         maxDepth,
      float&            MQ0Frac);

  /// \return Mean read depth of one alignment file over a set of regions on one chromosome, or a negative value
  /// if the regions are empty
  double getMeanRegionDepth(
      const unsigned bamIndex, const int32_t tid, const std::vector<known_pos_range2>& regions);

  /// For a deletion or tandem duplication large enough for read depth to be informative, estimate the depth
  /// within the variant and in its flanking regions for each sample
  void getSVDepthSupport(const SVCandidate& sv, SVScoreInfo& baseInfo);

  /// Apply all scoring models relevant to this event:
  ///
  /// \param junctionData one element describing each junction of an event, for normal (single-junction)
  /// candidates, the vector size should be one
  void computeAllScoreModels(
      const std::vector<JunctionCallInfo>& junctionData,
      SVModelScoreInfo&                    modelScoreInfo);

  /// Accumulate (model-agnostic) evidence for the SV alt/ref alleles
  ///
  /// This step gathers information (such as counts and allele likelihoods) to be used downstream by more
  /// specific scoring models.
  ///
  void getSVSupportingEvidence(
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& assemblyData,
      const SVCandidate&             sv,
      const SVId&                    svId,
      SVScoreInfo&                   ssInfo,
      SVEvidence&                    evidence,
      SVEvidenceWriterData&          svSupports);

  const CallOptionsShared       _callOpt;
  const CallOptionsSharedDeriv  _callDopt;
  const CallOptionsDiploid      _diploidOpt;
  const CallOptionsDiploidDeriv _diploidDopt;
  const ReadScannerOptions      _scanOpt;
  const SVRefinerOptions        _refineOpt;
  const ChromDepthFilterUtil    _dFilterDiploid;
  const SVLocusScanner&         _readScanner;

  const bam_header_info& _header;

  std::vector<streamPtr> _bamStreams;

  unsigned _sampleCount;
};
