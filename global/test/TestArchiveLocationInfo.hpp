/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef TESTARCHIVELOCATIONINFO_HPP_
#define TESTARCHIVELOCATIONINFO_HPP_

#include "ArchiveLocationInfo.hpp"

class TestArchiveLocationInfo : public CxxTest::TestSuite
{
public:
    void TestMethods()
    {
        //These throw because we are getting things before they are set.
        TS_ASSERT_THROWS_ANYTHING(ArchiveLocationInfo::GetArchiveDirectory());
        TS_ASSERT_THROWS_ANYTHING(ArchiveLocationInfo::GetMeshPathname());
        
        ArchiveLocationInfo::SetMeshPathname("archive_dir", "mesh_name");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshPathname(), "archive_dir/mesh_name");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveDirectory(), "archive_dir/");
        
        ArchiveLocationInfo::SetArchiveDirectory("new_archive_dir");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshPathname(), "new_archive_dir/mesh_name");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveDirectory(), "new_archive_dir/");
    }

    void TestProcessUniqueNaming()
    {
        ArchiveLocationInfo::SetArchiveDirectory("new_archive_dir");
        
        std::stringstream expected_filepath;
        expected_filepath << ArchiveLocationInfo::GetArchiveDirectory() << "fred";
        expected_filepath << "." << PetscTools::GetMyRank();
        
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetProcessUniqueFilePath("fred"), expected_filepath.str());
    }


};


#endif /*TESTARCHIVELOCATIONINFO_HPP_*/
