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


#ifndef TESTARCHIVINGHELPERCLASSES_HPP_
#define TESTARCHIVINGHELPERCLASSES_HPP_

#include <cstdlib> // For 'system'
#include <climits>
#include <sys/stat.h> // For chmod()
#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveLocationInfo.hpp"
#include "ArchiveOpener.hpp"
#include "ProcessSpecificArchive.hpp"
#include "PetscSetupAndFinalize.hpp"

// Save typing, and allow the use of these in cxxtest macros
typedef ArchiveOpener<boost::archive::text_iarchive, std::ifstream> InputArchiveOpener;
typedef ArchiveOpener<boost::archive::text_oarchive, std::ofstream> OutputArchiveOpener;

class TestArchivingHelperClasses : public CxxTest::TestSuite
{
public:

    void TestArchiveLocationInfoMethods() throw(Exception)
    {
        // These throw because we are getting things before they are set.
        TS_ASSERT_THROWS_THIS(ArchiveLocationInfo::GetArchiveDirectory(),"ArchiveLocationInfo::mDirPath has not been set");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshFilename(),"mesh"); //default value

        // To test exceptions (default value is now "mesh".)
        ArchiveLocationInfo::SetMeshFilename("");
        TS_ASSERT_THROWS_THIS(ArchiveLocationInfo::GetMeshFilename(),"ArchiveLocationInfo::mMeshFilename has not been set");

        ArchiveLocationInfo::SetMeshPathname("archive_dir", "mesh_name");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveDirectory(), "archive_dir/");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshFilename(), "mesh_name");

        ArchiveLocationInfo::SetArchiveDirectory("new_archive_dir");
        ArchiveLocationInfo::SetMeshFilename("new_mesh_name");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveDirectory(), "new_archive_dir/");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshFilename(), "new_mesh_name");

        // This throws because the full path isn't relative to test output
        TS_ASSERT_THROWS_THIS(ArchiveLocationInfo::GetArchiveRelativePath(),"Full path doesn\'t give a directory relative to CHASTE_TEST_OUTPUT");
        ArchiveLocationInfo::SetArchiveDirectory(OutputFileHandler::GetChasteTestOutputDirectory() + "relative_archive_dir");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveRelativePath(), "relative_archive_dir/");
        TS_ASSERT(ArchiveLocationInfo::GetIsDirRelativeToChasteTestOutput());
    }

    void TestArchiveLocationInfoProcessUniqueNaming() throw(Exception)
    {
        ArchiveLocationInfo::SetArchiveDirectory("new_archive_dir");

        std::stringstream expected_filepath;
        expected_filepath << ArchiveLocationInfo::GetArchiveDirectory() << "fred";
        expected_filepath << "." << PetscTools::GetMyRank();

        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetProcessUniqueFilePath("fred"), expected_filepath.str());
        
        std::string expected2 = ArchiveLocationInfo::GetArchiveDirectory() + "fred.12";
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetProcessUniqueFilePath("fred", 12), expected2);
    }

    void TestProcessSpecificArchive() throw(Exception)
    {
        TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_oarchive>::Get(),
                              "A ProcessSpecificArchive has not been set up.");
        TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_iarchive>::Get(),
                              "A ProcessSpecificArchive has not been set up.");

        // Set up an output archive pointer
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string arch_path = ArchiveLocationInfo::GetProcessUniqueFilePath("test.arch");
        std::ofstream ofs(arch_path.c_str());
        boost::archive::text_oarchive* p_arch = new boost::archive::text_oarchive(ofs);
        // Test the ProcessSpecificArchive Get and Set methods with this
        ProcessSpecificArchive<boost::archive::text_oarchive>::Set(p_arch);
        TS_ASSERT(ProcessSpecificArchive<boost::archive::text_oarchive>::Get()==p_arch);
        delete p_arch;

        // Set up an input archive pointer
        std::ifstream ifs(arch_path.c_str());
        boost::archive::text_iarchive* p_arch2 = new boost::archive::text_iarchive(ifs);
        // Test the ProcessSpecificArchive Get and Set methods with this
        ProcessSpecificArchive<boost::archive::text_iarchive>::Set(p_arch2);
        TS_ASSERT(ProcessSpecificArchive<boost::archive::text_iarchive>::Get()==p_arch2);
        delete p_arch2;

        // Clean up
        ProcessSpecificArchive<boost::archive::text_iarchive>::Set(NULL);
        ProcessSpecificArchive<boost::archive::text_oarchive>::Set(NULL);
    }

public:
    void TestArchiveOpenerReadAndWrite() throw(Exception)
    {
        std::string archive_dir = "archive";
        std::string archive_file = "archive_opener.arch";
        const unsigned test_int = 123;

        // Write
        {
            OutputArchiveOpener archive_opener_out(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = archive_opener_out.GetCommonArchive();
            boost::archive::text_oarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_oarchive>::Get();

            (*p_arch) & test_int; // All can write to the common archive - non-masters will write to /dev/null.
            (*p_process_arch) & test_int;

            // archive_opener_out will do a PetscTools::Barrier when it is destructed 
        }

        // Read
        {
            TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_oarchive>::Get(),
                                  "A ProcessSpecificArchive has not been set up.");
            TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_iarchive>::Get(),
                                  "A ProcessSpecificArchive has not been set up.");

            InputArchiveOpener archive_opener_in(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = archive_opener_in.GetCommonArchive();
            boost::archive::text_iarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_iarchive>::Get();

            unsigned test_int1, test_int2;
            (*p_arch) & test_int1;
            (*p_process_arch) & test_int2;
    
            TS_ASSERT_EQUALS(test_int1, test_int);
            TS_ASSERT_EQUALS(test_int2, test_int);
        }

        // Cover the case of an archive in the chaste folder (i.e. a path relative to the working directory).
        if (PetscTools::IsSequential())
        {
            // Make sure the file is not there before we do this
            int return_val = system("test -e testoutput/archive_opener.arch*");
            if (return_val == 0)
            {
                EXCEPTION("Unexpected file found; bailing out!");
            }
            { // Write
                OutputArchiveOpener archive_opener_relative("testoutput", "archive_opener.arch", false);
            }
            // Remove the file
            EXPECT0(system, "rm -f testoutput/archive_opener.arch*");
            
            { // Read
                InputArchiveOpener archive_opener_relative("apps/texttest/chaste/resume_bidomain/save_bidomain", "save_bidomain.arch", false);
            }
        }

        PetscTools::Barrier(); // Make sure all processes have finished this test before proceeding
    }

    // This test relies on TestArchiveOpenerReadAndWrite succeeding.
    void TestArchiveOpenerExceptions() throw(Exception)
    {
        std::string archive_dir = "archive";
        OutputFileHandler handler(archive_dir, false);
        handler.SetArchiveDirectory();
        std::string archive_base_name = "archive_opener.arch";

        // Remove the process-specific archive for this process
        std::string archive_path = ArchiveLocationInfo::GetProcessUniqueFilePath(archive_base_name);
        EXPECT0(system, "rm -f " + archive_path);
        TS_ASSERT_THROWS_CONTAINS(InputArchiveOpener archive_opener_in(archive_dir, archive_base_name),
                                  "Cannot load secondary archive file: ");

        // Remove the main archive
        if (PetscTools::AmMaster())
        {
            archive_path = handler.GetOutputDirectoryFullPath() + archive_base_name;
            EXPECT0(system, "rm -f " + archive_path);
        }
        PetscTools::Barrier();
        TS_ASSERT_THROWS_CONTAINS(InputArchiveOpener archive_opener_in(archive_dir, archive_base_name),
                                  "Cannot load main archive file: ");

        // Remove write permissions on the archive dir.
        if (PetscTools::AmMaster())
        {
            chmod(handler.GetOutputDirectoryFullPath().c_str(), 0444);
        }
        PetscTools::Barrier();
        // Now neither the master nor the slaves can write to their output files
        // This avoids hitting a PetscBarrier() in the ~ArchiveOpener() because they
        // all throw an error first.
        //
        // If this test starts hanging it is because these TS_ASSERT_THROWS_CONTAINS
        //  are not being thrown (rather than a real parallel calling problem).
        if (PetscTools::AmMaster())
        {
            TS_ASSERT_THROWS_CONTAINS(OutputArchiveOpener archive_opener_out(archive_dir, archive_base_name),
                                      "Failed to open main archive file for writing: ");
        }
        else
        {
            TS_ASSERT_THROWS_CONTAINS(OutputArchiveOpener archive_opener_out(archive_dir, archive_base_name),
                                      "Failed to open secondary archive file for writing: ");
        }
        PetscTools::Barrier();
        if (PetscTools::AmMaster())
        {   // Restore permissions on the folder before allowing processes to continue.
            chmod(handler.GetOutputDirectoryFullPath().c_str(), 0755);
        }
        PetscTools::Barrier();
    }
    
    void TestSpecifyingSecondaryArchive() throw (Exception)
    {
        std::string archive_dir = "archive";
        std::string archive_file = "specific_secondary.arch";
        const unsigned test_int = 321;
        const unsigned proc_id = PetscTools::GetMyRank();
        
        // Writing when specifying the secondary archive doesn't make sense
        {
            TS_ASSERT_THROWS_THIS(OutputArchiveOpener archive_opener_out(archive_dir, archive_file, true, UINT_MAX),
                                  "Specifying the secondary archive file ID doesn't make sense when writing.");
        }
        
        // Write normally so we can test reading
        {
            OutputArchiveOpener archive_opener_out(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = archive_opener_out.GetCommonArchive();
            boost::archive::text_oarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_oarchive>::Get();

            (*p_arch) & test_int; // All can write to the common archive - non-masters will write to /dev/null.
            (*p_process_arch) & proc_id;

            // archive_opener_out will do a PetscTools::Barrier when it is destructed 
        }

        // Read
        {
            TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_oarchive>::Get(),
                                  "A ProcessSpecificArchive has not been set up.");
            TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_iarchive>::Get(),
                                  "A ProcessSpecificArchive has not been set up.");

            InputArchiveOpener archive_opener_in(archive_dir, archive_file, true, 0);
            boost::archive::text_iarchive* p_arch = archive_opener_in.GetCommonArchive();
            boost::archive::text_iarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_iarchive>::Get();

            unsigned test_int1, test_int2;
            (*p_arch) & test_int1;
            (*p_process_arch) & test_int2;
    
            TS_ASSERT_EQUALS(test_int1, test_int);
            TS_ASSERT_EQUALS(test_int2, 0u);
        }
    }
};


#endif /*TESTARCHIVINGHELPERCLASSES_HPP_*/
