/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "TidyTemplatedExportIdentifier.hpp"
#include <algorithm>

std::string TidyTemplatedExportIdentifier(std::string identifier)
{
    // First remove spaces, so identifier now takes the form "pack<void(NameOfDerivedType<DIM>)>::type"
    std::string::iterator end_pos = std::remove(identifier.begin(), identifier.end(), ' ');
    identifier.erase(end_pos, identifier.end());

    // Then remove "pack<void(", so identifier now takes the form "NameOfDerivedType<DIM>)>::type"
    const std::string s1 = "pack<void(";
    std::string::size_type i = identifier.find(s1);
    if (i != identifier.npos)
    {
        identifier.erase(i, s1.length());
    }

    // Then replace "<" with "-", so identifier now takes the form "NameOfDerivedType-DIM>)>::type"
    const std::string s2 = "<";
    const std::string s3 = "-";
    i = identifier.find(s2);
    if (i != identifier.npos)
    {
        identifier.replace(i, s2.length(), s3);
    }

    // Finally remove ">)>::type", so that identifier now takes the form "NameOfDerivedType-DIM"
    const std::string s4 = ">)>::type";
    i = identifier.find(s4);
    if (i != identifier.npos)
    {
        identifier.erase(i, s4.length());
    }

    return identifier;
}
