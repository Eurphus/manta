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
/// \author Xiaoyu Chen
///

#include "boost/test/unit_test.hpp"

#include "IterativeAssembler.cpp"

BOOST_AUTO_TEST_SUITE(test_IterativeAssembler)

BOOST_AUTO_TEST_CASE(test_CircleDetector)
{
  IterativeAssemblerOptions assembleOpt;
  str_uint_map_t            wordCount;
  std::set<std::string>     repeatWords;

  wordCount["TACCA"] = 3;
  wordCount["CCACC"] = 3;
  wordCount["CACCA"] = 3;
  wordCount["ACCAC"] = 3;
  wordCount["CCACA"] = 3;
  wordCount["CACAC"] = 3;
  wordCount["ACACA"] = 3;
  wordCount["AAAAA"] = 2;

  getRepeatKmers(assembleOpt, wordCount, repeatWords);

  // the first circle
  BOOST_REQUIRE_EQUAL(repeatWords.count("ACCAC"), 1u);
  BOOST_REQUIRE_EQUAL(repeatWords.count("CACCA"), 1u);
  BOOST_REQUIRE_EQUAL(repeatWords.count("CCACC"), 1u);

  BOOST_REQUIRE_EQUAL(repeatWords.count("TACCA"), 0u);
  BOOST_REQUIRE_EQUAL(repeatWords.count("CCACA"), 0u);

  // the second circle
  BOOST_REQUIRE_EQUAL(repeatWords.count("CACAC"), 1u);
  BOOST_REQUIRE_EQUAL(repeatWords.count("ACACA"), 1u);

  // homopolymer: self-circle
  BOOST_REQUIRE_EQUAL(repeatWords.count("AAAAA"), 1u);
}

BOOST_AUTO_TEST_CASE(test_BasicAssembler)
{
  // test simple assembly functions at a single word size:
  IterativeAssemblerOptions assembleOpt;

  assembleOpt.minWordLength = 6;
  assembleOpt.maxWordLength = 6;
  assembleOpt.minCoverage   = 2;

  AssemblyReadInput reads;

  // clang-format off
  reads.emplace_back("ACGTGTATTACC");
  reads.emplace_back(  "GTGTATTACCTA");
  reads.emplace_back(      "ATTACCTAGTAC");
  reads.emplace_back(        "TACCTAGTACTC");
  reads.emplace_back("123456789123");
  // clang-format on

  AssemblyReadOutput readInfo;
  Assembly           contigs;

  runIterativeAssembler(assembleOpt, reads, readInfo, contigs);

  BOOST_REQUIRE_EQUAL(contigs.size(), 1u);
  BOOST_REQUIRE_EQUAL(contigs[0].seq, "GTGTATTACCTAGTAC");
  for (unsigned i(0); i < 4; ++i) {
    BOOST_REQUIRE(readInfo[i].isUsed);
    BOOST_REQUIRE_EQUAL(readInfo[i].contigIds[0], 0u);
  }
  BOOST_REQUIRE(!readInfo[4].isUsed);
}

BOOST_AUTO_TEST_CASE(test_IterativeKmer)
{
  // test simple assembly functions at a single word size:
  IterativeAssemblerOptions assembleOpt;

  assembleOpt.minWordLength = 3;
  assembleOpt.maxWordLength = 9;
  assembleOpt.wordStepSize  = 3;
  assembleOpt.minCoverage   = 1;

  AssemblyReadInput reads;

  // clang-format off
  reads.emplace_back("ACACACACGATG");
  reads.emplace_back(        "GATGTCTCTCTC");
  reads.emplace_back("123456789123");
  // clang-format on

  AssemblyReadOutput readInfo;
  Assembly           contigs;

  runIterativeAssembler(assembleOpt, reads, readInfo, contigs);

  BOOST_REQUIRE_EQUAL(contigs.size(), 1u);
  BOOST_REQUIRE_EQUAL(contigs[0].seq, "ACACACACGATGTCTCTCTC");
  for (unsigned i(0); i < 2; ++i) {
    BOOST_REQUIRE(readInfo[i].isUsed);
    BOOST_REQUIRE_EQUAL(readInfo[i].contigIds[0], 0u);
  }
  BOOST_REQUIRE(!readInfo[2].isUsed);
}

BOOST_AUTO_TEST_CASE(test_branching_basic)
{
  // test simple assembly functions at a single word size:
  IterativeAssemblerOptions assembleOpt;

  assembleOpt.minWordLength   = 6;
  assembleOpt.maxWordLength   = 6;
  assembleOpt.minCoverage     = 1;
  assembleOpt.minSupportReads = 1;
  assembleOpt.minUnusedReads  = 1;

  AssemblyReadInput reads;

  // clang-format off
  reads.emplace_back("ATATAGACGATG");
  reads.emplace_back(      "ACGATGTCTATCTT");
  reads.emplace_back(      "ACGATGTTGGCCTT");
  // clang-format on

  AssemblyReadOutput readInfo;
  Assembly           contigs;

  runIterativeAssembler(assembleOpt, reads, readInfo, contigs);

  BOOST_REQUIRE_EQUAL(contigs.size(), 2u);
  BOOST_REQUIRE_EQUAL(contigs[0].seq, "ATATAGACGATGTCTATCTT");
  BOOST_REQUIRE_EQUAL(contigs[1].seq, "ATATAGACGATGTTGGCCTT");

  BOOST_REQUIRE(readInfo[0].isUsed);
  BOOST_REQUIRE_EQUAL(readInfo[0].contigIds[0], 0u);
  BOOST_REQUIRE_EQUAL(readInfo[0].contigIds[1], 1u);

  BOOST_REQUIRE(readInfo[1].isUsed);
  BOOST_REQUIRE_EQUAL(readInfo[1].contigIds[0], 0u);

  BOOST_REQUIRE(readInfo[2].isUsed);
  BOOST_REQUIRE_EQUAL(readInfo[2].contigIds[0], 1u);
}

BOOST_AUTO_TEST_CASE(test_branching_iterative)
{
  // test simple assembly functions at a single word size:
  IterativeAssemblerOptions assembleOpt;

  assembleOpt.minWordLength   = 3;
  assembleOpt.maxWordLength   = 9;
  assembleOpt.wordStepSize    = 3;
  assembleOpt.minCoverage     = 1;
  assembleOpt.minSupportReads = 1;
  assembleOpt.minUnusedReads  = 1;

  AssemblyReadInput reads;

  // clang-format off
  reads.emplace_back("ACACACACGATG");
  reads.emplace_back(        "GATGGCCCCCCC");
  reads.emplace_back(        "GATGTCTCTCTC");
  // clang-format on

  AssemblyReadOutput readInfo;
  Assembly           contigs;

  runIterativeAssembler(assembleOpt, reads, readInfo, contigs);

  BOOST_REQUIRE_EQUAL(contigs.size(), 2u);
  BOOST_REQUIRE_EQUAL(contigs[0].seq, "ACACACACGATGGCCCCCCC");
  BOOST_REQUIRE_EQUAL(contigs[1].seq, "ACACACACGATGTCTCTCTC");

  BOOST_REQUIRE(readInfo[0].isUsed);
  BOOST_REQUIRE_EQUAL(readInfo[0].contigIds[0], 0u);
  BOOST_REQUIRE_EQUAL(readInfo[0].contigIds[1], 1u);

  BOOST_REQUIRE(readInfo[1].isUsed);
  BOOST_REQUIRE_EQUAL(readInfo[1].contigIds[0], 0u);

  BOOST_REQUIRE(readInfo[2].isUsed);
  BOOST_REQUIRE_EQUAL(readInfo[2].contigIds[0], 1u);
}

BOOST_AUTO_TEST_CASE(test_RefineConservativeRange)
{
  IterativeAssemblerOptions assembleOpt;
  assembleOpt.refineSeedLength        = 8;
  assembleOpt.refineMinOverlap        = 20;
  assembleOpt.refineMaxMismatchFrac   = 0.04;
  assembleOpt.minConservativeCoverage = 2;

  const std::string contigSeq("TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACG");
  BOOST_REQUIRE_EQUAL(contigSeq.size(), 60u);

  AssemblyReadInput reads;
  reads.push_back(contigSeq.substr(0, 45));
  reads.push_back(contigSeq.substr(10, 50));

  // a third read covering [5,50) with a single basecall error
  {
    std::string read(contigSeq.substr(5, 45));
    read[20] = ((read[20] == 'A') ? 'C' : 'A');
    reads.push_back(read);
  }

  // a read which shares a seed with the contig but otherwise disagrees with it must not count
  {
    std::string read(contigSeq.substr(0, 45));
    for (unsigned i(12); i < 45; i += 3) read[i] = ((read[i] == 'G') ? 'T' : 'G');
    reads.push_back(read);
  }

  // an unrelated read
  reads.emplace_back("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

  {
    // without the read carrying one error, only [10,45) is covered by two reads
    AssembledContig contig;
    contig.seq = contigSeq;
    refineContigConservativeRange(assembleOpt, reads, 2, contig);
    BOOST_REQUIRE_EQUAL(contig.conservativeRange.begin_pos(), 10);
    BOOST_REQUIRE_EQUAL(contig.conservativeRange.end_pos(), 45);
  }

  {
    // the read with a single error extends two-read coverage to [5,50), the remaining reads change nothing
    AssembledContig contig;
    contig.seq = contigSeq;
    refineContigConservativeRange(assembleOpt, reads, reads.size(), contig);
    BOOST_REQUIRE_EQUAL(contig.conservativeRange.begin_pos(), 5);
    BOOST_REQUIRE_EQUAL(contig.conservativeRange.end_pos(), 50);
  }

  {
    // a range found during contig extension is extended, but never reduced
    AssembledContig contig;
    contig.seq = contigSeq;
    contig.conservativeRange.set_range(0, 20);
    refineContigConservativeRange(assembleOpt, reads, reads.size(), contig);
    BOOST_REQUIRE_EQUAL(contig.conservativeRange.begin_pos(), 0);
    BOOST_REQUIRE_EQUAL(contig.conservativeRange.end_pos(), 50);
  }

  {
    // no change when fewer than two reads can be placed
    AssembledContig contig;
    contig.seq = contigSeq;
    contig.conservativeRange.set_range(7, 9);
    refineContigConservativeRange(assembleOpt, reads, 1, contig);
    BOOST_REQUIRE_EQUAL(contig.conservativeRange.begin_pos(), 7);
    BOOST_REQUIRE_EQUAL(contig.conservativeRange.end_pos(), 9);
  }
}

BOOST_AUTO_TEST_SUITE_END()
