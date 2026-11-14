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

#include "boost/test/unit_test.hpp"

#include "blt_exception.hpp"
#include "qscore_cache.hpp"
#include "qscore_snp.hpp"

BOOST_AUTO_TEST_SUITE(qscore_cache_test_suite)

// Quality scores above the cached range are valid input and must be treated as the highest cached score
BOOST_AUTO_TEST_CASE(test_high_qscore_is_capped)
{
  static const int maxQ(qphred_cache::MAX_QSCORE);

  for (const int highQ : {maxQ + 1, 93, 255}) {
    BOOST_REQUIRE_EQUAL(qphred_cache::get_error_prob(highQ), qphred_cache::get_error_prob(maxQ));
    BOOST_REQUIRE_EQUAL(qphred_cache::get_ln_error_prob(highQ), qphred_cache::get_ln_error_prob(maxQ));
    BOOST_REQUIRE_EQUAL(
        qphred_cache::get_ln_comp_error_prob(highQ), qphred_cache::get_ln_comp_error_prob(maxQ));
    BOOST_REQUIRE_EQUAL(qphred_cache::get_mapped_qscore(highQ, 60), qphred_cache::get_mapped_qscore(maxQ, 60));
  }

  const qscore_snp qscoreSnp(0.001);
  BOOST_REQUIRE_EQUAL(qscoreSnp.qphred_to_error_prob(93), qscoreSnp.qphred_to_error_prob(maxQ));
  BOOST_REQUIRE_EQUAL(qscoreSnp.qphred_to_ln_error_prob(93), qscoreSnp.qphred_to_ln_error_prob(maxQ));
  BOOST_REQUIRE_EQUAL(
      qscoreSnp.qphred_to_ln_comp_error_prob(93), qscoreSnp.qphred_to_ln_comp_error_prob(maxQ));

  // scores inside the cached range are unchanged
  BOOST_REQUIRE_CLOSE(qphred_cache::get_error_prob(30), 0.001, 0.0000001);
  BOOST_REQUIRE(qphred_cache::get_error_prob(maxQ - 1) > qphred_cache::get_error_prob(maxQ));
}

// Negative quality scores remain an error
BOOST_AUTO_TEST_CASE(test_negative_qscore_throws)
{
  BOOST_CHECK_THROW(qphred_cache::get_error_prob(-1), blt_exception);

  const qscore_snp qscoreSnp(0.001);
  BOOST_CHECK_THROW(qscoreSnp.qphred_to_error_prob(-1), blt_exception);
}

BOOST_AUTO_TEST_SUITE_END()
