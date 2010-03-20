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

#ifndef CARDIACSIMULATIONARCHIVER_HPP_
#define CARDIACSIMULATIONARCHIVER_HPP_

#include <string>
#include <fstream>

// Must be included before any other serialisation headers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Exception.hpp"
#include "ArchiveOpener.hpp"
#include "OutputFileHandler.hpp"
#include "ArchiveLocationInfo.hpp"
#include "DistributedVectorFactory.hpp"
#include "PetscTools.hpp"

/**
 *  CardiacSimulationArchiver is a helper class for checkpointing of cardiac simulations.
 *
 *  The class is templated over the class defining the simulation (i.e. MonodomainProblem)
 */
template<class PROBLEM_CLASS>
class CardiacSimulationArchiver
{
public:
    /**
     * Archives a simulation in the directory specified.
     *
     * @note Must be called collectively, i.e. by all processes.
     *
     * @param simulationToArchive object defining the simulation to archive
     * @param rDirectory directory where the multiple files defining the checkpoint will be stored
     *     (relative to CHASTE_TEST_OUTPUT)
     * @param clearDirectory whether the directory needs to be cleared or not.
     */
    static void Save(PROBLEM_CLASS& simulationToArchive, const std::string& rDirectory, bool clearDirectory=true);


    /**
     * Unarchives a simulation from the directory specified.
     *
     * Does a migrate if necessary (this is actually just a wrapper around the
     * Migrate method now).
     *
     * @note Must be called collectively, i.e. by all processes.
     *
     * @param rDirectory directory where the multiple files defining the checkpoint are located
     *     (relative to CHASTE_TEST_OUTPUT)
     */
    static PROBLEM_CLASS* Load(const std::string& rDirectory);


    /**
     * Load a simulation, saved by any number of processes, into the correct
     * number of processes for those currently launched.
     *
     * @note Must be called collectively, i.e. by all processes.
     *
     * Uses the DistributedVectorFactory saved in the process 0 archive to work out
     * how many secondary archive files to read, and loads the cells and boundary
     * conditions from these too.
     *
     * Uses a dumb partition to work out how to distribute the mesh and cells over
     * the processes.
     * \todo Allow the use of METIS to give a better partitioning
     *
     * @param rDirectory directory where the multiple files defining the checkpoint are located
     *     (relative to CHASTE_TEST_OUTPUT)
     */
    static PROBLEM_CLASS* Migrate(const std::string &rDirectory);
};


template<class PROBLEM_CLASS>
void CardiacSimulationArchiver<PROBLEM_CLASS>::Save(PROBLEM_CLASS& simulationToArchive,
                                                    const std::string &rDirectory,
                                                    bool clearDirectory)
{
    // Clear directory if requested (and make sure it exists)
    OutputFileHandler handler(rDirectory, clearDirectory);

    // Nest the archive writing, so the ArchiveOpener goes out of scope before
    // the method ends.
    {
        // Open the archive files
        ArchiveOpener<boost::archive::text_oarchive, std::ofstream> archive_opener(rDirectory, "archive.arch", true);
        boost::archive::text_oarchive* p_main_archive = archive_opener.GetCommonArchive();

        // And save
        PROBLEM_CLASS* const p_simulation_to_archive = &simulationToArchive;
        (*p_main_archive) & p_simulation_to_archive;
    }

    // Write the info file
    if (PetscTools::AmMaster())
    {
        std::string info_path = handler.GetOutputDirectoryFullPath() + "archive.info";
        std::ofstream info_file(info_path.c_str());
        if (!info_file.is_open())
        {
            // Avoid deadlock...
            PetscTools::ReplicateBool(true);
            EXCEPTION("Unable to open archive information file: " + info_path);
        }
        PetscTools::ReplicateBool(false);
        unsigned archive_version = 0; /// \todo #1026 get a real version number!
        info_file << PetscTools::GetNumProcs() << " " << archive_version;
    }
    else
    {
        bool master_threw = PetscTools::ReplicateBool(false);
        if (master_threw)
        {
            EXCEPTION("Unable to open archive information file");
        }
    }
    // Make sure everything is written before any process continues.
    PetscTools::Barrier("CardiacSimulationArchiver::Save");
}

template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::Load(const std::string &rDirectory)
{
    return CardiacSimulationArchiver<PROBLEM_CLASS>::Migrate(rDirectory);
}


template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::Migrate(const std::string &rDirectory)
{
    // Load the info file
    OutputFileHandler handler(rDirectory, false);
    std::string info_path = handler.GetOutputDirectoryFullPath() + "archive.info";
    std::ifstream info_file(info_path.c_str());
    if (!info_file.is_open())
    {
        EXCEPTION("Unable to open archive information file: " + info_path);
    }
    unsigned num_procs, archive_version;
    info_file >> num_procs >> archive_version;

    PROBLEM_CLASS *p_unarchived_simulation;

    // Avoid the DistributedVectorFactory throwing a 'wrong number of processes' exception when loading,
    // and make it get the original DistributedVectorFactory from the archive so we can compare against
    // num_procs.
    DistributedVectorFactory::SetCheckNumberOfProcessesOnLoad(false);
    // Put what follows in a try-catch to make sure we reset this
    try
    {
        if (num_procs == PetscTools::GetNumProcs())
        {
            // We're not actually doing a migrate, and will re-use exactly the same mesh
            // partitioning etc. from before, so don't need any of the LoadExtraArchive
            // magic.  Indeed, we mustn't use it, or the mesh will get confused about
            // which nodes it owns.
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> archive_opener(rDirectory, "archive.arch", true);
            boost::archive::text_iarchive* p_main_archive = archive_opener.GetCommonArchive();
            (*p_main_archive) >> p_unarchived_simulation;

            // Paranoia checks
            DistributedVectorFactory* p_factory = p_unarchived_simulation->rGetMesh().GetDistributedVectorFactory();
            assert(p_factory != NULL);
            unsigned original_num_procs = p_factory->GetOriginalFactory()->GetNumProcs();
            if (original_num_procs != num_procs)
            {
                NEVER_REACHED;
                //assert(original_num_procs == num_procs);
            }
        }
        else
        {
            // We are migrating, and must re-distribute nodes, cells, etc.

            // Load the master and process-0 archive files.
            // This will also set up ArchiveLocationInfo for us.
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> archive_opener(rDirectory, "archive.arch", true, 0u);
            boost::archive::text_iarchive* p_main_archive = archive_opener.GetCommonArchive();
            (*p_main_archive) >> p_unarchived_simulation;

            // Work out how many more files to load
            DistributedVectorFactory* p_factory = p_unarchived_simulation->rGetMesh().GetDistributedVectorFactory();
            assert(p_factory != NULL);
            unsigned original_num_procs = p_factory->GetOriginalFactory()->GetNumProcs();
            assert(original_num_procs == num_procs); // Paranoia

            // Merge in the extra data
            for (unsigned archive_num=1; archive_num<original_num_procs; archive_num++)
            {
                std::string archive_path = ArchiveLocationInfo::GetProcessUniqueFilePath("archive.arch", archive_num);
                std::ifstream ifs(archive_path.c_str());
                boost::archive::text_iarchive archive(ifs);
                p_unarchived_simulation->LoadExtraArchive(archive, archive_version);
            }
        }
    }
    catch (Exception &e)
    {
        DistributedVectorFactory::SetCheckNumberOfProcessesOnLoad(true);
        throw e;
    }

    // Done.
    DistributedVectorFactory::SetCheckNumberOfProcessesOnLoad(true);
    return p_unarchived_simulation;
}

#endif /*CARDIACSIMULATIONARCHIVER_HPP_*/
