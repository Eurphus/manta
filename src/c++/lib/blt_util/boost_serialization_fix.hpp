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
/// \brief Include this before any boost serialization header for a standard container
///
/// In boost 1.74 and 1.75 the serialization headers for standard containers (vector.hpp, map.hpp, set.hpp...)
/// use boost::serialization::library_version_type without including its definition, so that they fail to
/// compile unless the definition happens to be included first. The header providing the definition does not
/// exist before boost 1.74.

#pragma once

#include "boost/version.hpp"

#if (BOOST_VERSION >= 107400)
#include "boost/serialization/library_version_type.hpp"
#endif
