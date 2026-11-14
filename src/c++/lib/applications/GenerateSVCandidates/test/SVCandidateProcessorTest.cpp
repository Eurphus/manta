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

#include <fstream>
#include "boost/make_unique.hpp"
#include "boost/test/unit_test.hpp"
#include "htsapi/vcf_streamer.hpp"
#include "manta/BamStreamerUtils.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"
#include "test/testSVLocusScanner.hpp"
#include "test/testSVLocusUtil.hpp"

#include "SVCandidateProcessor.cpp"  // cpp include is required to test 'checkJunctionsToFilter"
#include "SVCandidateProcessor.hpp"
#include "options/SVLocusSetOptions.hpp"
#include "svgraph/SVLocusSet.hpp"

/// The whole purpose of this test file is given an SV candidate, whether it is correctly writing
/// SV information in the corresponding vcf file.

BOOST_AUTO_TEST_SUITE(SVCandidateProcessor_test_suite)

// A candidate sv needs to be checked after assembly whether it is good candidate sv or not.
// A candidate sv is not good if either of the following three cases is  satisfied:
// 1. Spanning evidence count of a candidate is less than min candidate spanning count threshold (default 3).
// 2. Non-spanning low-res candidate went into assembly but did not produce a successful contig alignment
//    that means SV is imprecise.
// 3. Variant size is smaller than minCandidateVariantSize (default is 10)
// A successful case is designed in test_Diploid.
BOOST_AUTO_TEST_CASE(test_JunctionFilter)
{
  GSCOptions options;
  // Creating SV candidate
  SVCandidate candidate1;
  candidate1.setPrecise();
  candidate1.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate1.bp1.interval = GenomeInterval(0, 40, 50);
  candidate1.bp1.lowresEvidence.add(0, 1);
  candidate1.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate1.bp2.interval = GenomeInterval(1, 65, 75);
  // Adding a pair evidence to this candidate
  candidate1.bp2.lowresEvidence.add(0, 1);
  SVMultiJunctionCandidate junctionCandidate;
  junctionCandidate.junction = {candidate1};
  std::vector<bool> isInputJunctionFiltered1;
  isInputJunctionFiltered1.resize(1);
  SVCandidateAssemblyData candidateAssemblyData;
  candidateAssemblyData.isCandidateSpanning = true;
  std::vector<SVCandidateAssemblyData> assemblyData;
  assemblyData.push_back(candidateAssemblyData);
  // Case-1 is designed where number of spanning evidence observations required is 3.
  // But here spanning observation count is 1.
  checkJunctionsToFilter(junctionCandidate, assemblyData, isInputJunctionFiltered1, options);
  BOOST_REQUIRE_EQUAL(isInputJunctionFiltered1[0], true);

  // Case-2 is designed
  // In this case a non-spanning low-res candidate went into assembly but
  // did not produce a successful contig alignment.
  SVCandidate candidate2;
  candidate2.bp1.state    = SVBreakendState::COMPLEX;
  candidate2.bp1.interval = GenomeInterval(0, 40, 50);
  candidate2.bp1.lowresEvidence.add(0, 1);
  candidate2.bp2.state    = SVBreakendState::UNKNOWN;
  candidate2.bp2.interval = GenomeInterval(0, 65, 75);
  candidate2.bp2.lowresEvidence.add(0, 1);
  junctionCandidate.junction.clear();
  junctionCandidate.junction.push_back(candidate2);
  candidateAssemblyData.isCandidateSpanning = false;
  assemblyData.clear();
  assemblyData.push_back(candidateAssemblyData);
  std::vector<bool> isInputJunctionFiltered2;
  isInputJunctionFiltered2.resize(1);
  checkJunctionsToFilter(junctionCandidate, assemblyData, isInputJunctionFiltered2, options);
  BOOST_REQUIRE_EQUAL(isInputJunctionFiltered2[0], true);

  // Case-3 is designed.
  // variant size is less than min variant size(40)
  options.minCandidateSpanningCount = 1;
  SVCandidate candidate3;
  candidate3.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate3.bp1.interval = GenomeInterval(0, 40, 50);
  candidate3.bp1.lowresEvidence.add(0, 1);
  candidate3.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate3.bp2.interval = GenomeInterval(0, 65, 75);
  candidate3.bp2.lowresEvidence.add(0, 1);
  junctionCandidate.junction.clear();
  junctionCandidate.junction.push_back(candidate2);
  candidateAssemblyData.isCandidateSpanning = true;
  assemblyData.clear();
  assemblyData.push_back(candidateAssemblyData);
  options.scanOpt.minCandidateVariantSize = 40;
  std::vector<bool> isInputJunctionFiltered3;
  isInputJunctionFiltered3.resize(1);
  // Variant size = 65 - 40 -1 = 24 which is less than 40.
  checkJunctionsToFilter(junctionCandidate, assemblyData, isInputJunctionFiltered3, options);
  BOOST_REQUIRE_EQUAL(isInputJunctionFiltered3[0], true);
}

// Create Temporary bam streams of a bam file which contains
// three inter-chromosomal read pairs.
struct BamStream {
  BamStream()
  {
    const bam_header_info bamHeader(buildTestBamHeader());

    std::string querySeq1 = "GTCTATCACCCTATTAACCACTCACGGGAGAAAAA";
    std::string querySeq2 = "AAAAATGCTCATCAGTTGATGATACGCCCGAGCAGATGCCAACACAGCAGCCATTCAAGAGTCA";
    bam_record  bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 8, 1, 69, 35, 15, "35M", querySeq1, 0);
    bamRecord1.set_qname("Read-1");
    bamRecord1.toggle_is_first();
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 8, 1, 69, 35, 15, "35M", querySeq1, 0);
    bamRecord2.set_qname("Read-2");
    bamRecord2.toggle_is_first();
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 8, 1, 69, 35, 15, "35M", querySeq1, 0);
    bamRecord3.set_qname("Read-3");
    bamRecord3.toggle_is_first();
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 1, 69, 0, 8, 64, 50, "64M", querySeq2, 0);
    bamRecord4.toggle_is_mate_fwd_strand();
    bamRecord4.toggle_is_fwd_strand();
    bamRecord4.set_qname("Read-1");
    bamRecord4.toggle_is_second();
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5, 1, 69, 0, 8, 64, 50, "64M", querySeq2, 0);
    bamRecord5.set_qname("Read-2");
    bamRecord5.toggle_is_mate_fwd_strand();
    bamRecord5.toggle_is_fwd_strand();
    bamRecord5.toggle_is_second();
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 1, 69, 0, 8, 64, 50, "64M", querySeq2, 0);
    bamRecord6.toggle_is_mate_fwd_strand();
    bamRecord6.toggle_is_fwd_strand();
    bamRecord6.set_qname("Read-3");
    bamRecord6.toggle_is_second();
    readsToAdd.push_back(bamRecord1);
    readsToAdd.push_back(bamRecord2);
    readsToAdd.push_back(bamRecord3);
    readsToAdd.push_back(bamRecord4);
    readsToAdd.push_back(bamRecord5);
    readsToAdd.push_back(bamRecord6);
    bamFileName = _bamFilename();
    buildTestBamFile(bamHeader, readsToAdd, bamFileName);

    const std::string                          referenceFilename = getTestReferenceFilename();
    std::vector<std::string>                   bamFilenames      = {bamFileName};
    std::vector<std::shared_ptr<bam_streamer>> bamStreams;
    openBamStreams(referenceFilename, bamFilenames, bamStreams);
    bamStream = bamStreams[0];
  }

  std::shared_ptr<bam_streamer> bamStream;
  std::vector<bam_record>       readsToAdd;
  std::string                   bamFileName;

private:
  const std::string&     _bamFilename() const { return _bamFilenameMaker.getFilename(); }
  const BamFilenameMaker _bamFilenameMaker;
};

BOOST_FIXTURE_TEST_SUITE(SVCandidateProcessor_test_suite, BamStream)

// For diploid caller, VCF records are written in diploid file.
// Verify those VCF records in diploid file.
BOOST_AUTO_TEST_CASE(test_Diploid)
{
  const bam_header_info bamHeader(buildTestBamHeader());
  TestFilenameMaker     filenameMaker1;
  TestFilenameMaker     filenameMaker3;
  TestFilenameMaker     filenameMaker4;
  GSCOptions            options;

  // Assembly options
  options.refineOpt.spanningAssembleOpt.minWordLength = 3;
  options.refineOpt.spanningAssembleOpt.maxWordLength = 9;
  options.refineOpt.spanningAssembleOpt.wordStepSize  = 3;
  options.refineOpt.smallSVAssembleOpt.minWordLength  = 3;
  options.refineOpt.smallSVAssembleOpt.maxWordLength  = 9;
  options.refineOpt.smallSVAssembleOpt.wordStepSize   = 3;
  options.alignFileOpt.alignmentFilenames             = {bamFileName};
  options.referenceFilename                           = getTestReferenceFilename();
  options.edgeRuntimeFilename                         = filenameMaker1.getFilename();
  // VCF records should be written in the diploid file
  options.diploidOutputFilename        = filenameMaker3.getFilename();
  options.candidateOutputFilename      = filenameMaker4.getFilename();
  options.minCandidateSpanningCount    = 1;
  options.minScoredVariantSize         = 20;
  options.diploidOpt.minOutputAltScore = 0;

  TestStatsFileMaker statsFileMaker;
  options.statsFilename = statsFileMaker.getFilename();
  auto                edgeTrackerPtr(std::make_shared<EdgeRuntimeTracker>(options.edgeRuntimeFilename));
  GSCEdgeStatsManager edgeStatMan;

  // Creating SV candidate
  SVCandidate candidate1;
  candidate1.setPrecise();
  candidate1.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate1.bp1.interval = GenomeInterval(0, 40, 50);
  candidate1.bp1.lowresEvidence.add(0, 1);
  candidate1.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate1.bp2.interval = GenomeInterval(1, 65, 75);
  candidate1.bp2.lowresEvidence.add(0, 1);
  SVCandidateAssemblyData         candidateAssemblyData1;
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  std::string                     programName = "Manta";
  std::string                     version     = "Test:0:1";
  EdgeInfo                        edgeInfo;
  edgeInfo.nodeIndex1 = 0;
  edgeInfo.nodeIndex2 = 1;
  SVMultiJunctionCandidate junctionCandidate;
  junctionCandidate.junction = {candidate1};
  std::vector<SVMultiJunctionCandidate> mjSvs;
  mjSvs.push_back(junctionCandidate);
  // Adding fragment information
  SVCandidateSetData                         svData;
  SVCandidateSetSequenceFragmentSampleGroup& group = svData.getDataGroup(0);
  group.add(bamHeader, readsToAdd[0], false, true, true);
  group.add(bamHeader, readsToAdd[3], false, true, true);
  SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
  group.begin()->svLink.push_back(association);

  // Generate SV locus graph file
  // Two nodes are created.
  SVLocus locus1;
  locusAddPair(locus1, 0, 40, 50, 1, 65, 75);
  SVLocus locus2;
  locusAddPair(locus2, 1, 65, 75, 0, 40, 50);
  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 10;
  SVLocusSet set1(sopt, bamHeader, {bamFileName});
  set1.merge(locus1);
  set1.merge(locus2);
  set1.checkState(true, true);
  TestFilenameMaker testFilenameMaker;
  std::string       graphFilename = testFilenameMaker.getFilename();
  const char*       testFilenamePtr(graphFilename.c_str());
  // serialize
  set1.save(testFilenamePtr);

  // block-scope svWriter to force file flush at end of scope:
  {
    const SVWriter       svWriter(options, bamHeader, programName.c_str(), version.c_str());
    auto                 svEvidenceWriterSharedData(std::make_shared<SVEvidenceWriterSharedData>(options));
    SVCandidateProcessor candidateProcessor(
        options,
        scanner.operator*(),
        set1,
        svWriter,
        svEvidenceWriterSharedData,
        edgeTrackerPtr,
        edgeStatMan);
    candidateProcessor.evaluateCandidates(edgeInfo, mjSvs, svData);
  }

  // Check output vcf file
  std::ifstream diploidFile(options.diploidOutputFilename);
  std::string   line;
  int           count(0);
  // Record-1:    chrFoo	41	MantaBND:0:0:1:0:0:0:0	C	C[chrBar:66[	.	.
  //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:1;CIPOS=0,9;HOMLEN=9;HOMSEQ=TCCATGCAT;BND_PAIR_COUNT=0;PAIR_COUNT=1
  // Record-2:    chrBar	66	MantaBND:0:0:1:0:0:0:1	T	]chrFoo:41]T	.	.
  //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:0;CIPOS=0,9;HOMLEN=9;HOMSEQ=ATCCTTCTT;BND_PAIR_COUNT=0;PAIR_COUNT=1
  if (diploidFile.is_open()) {
    while (std::getline(diploidFile, line)) {
      if (line.find("#") == std::string::npos) {
        count++;
      }
      if (count == 1)
        BOOST_REQUIRE(line.find("C[chrBar:66[") != std::string::npos);
      else if (count == 2)
        BOOST_REQUIRE(line.find("]chrFoo:41]T") != std::string::npos);
    }
    diploidFile.close();
  }
  BOOST_REQUIRE_EQUAL(count, 2);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
