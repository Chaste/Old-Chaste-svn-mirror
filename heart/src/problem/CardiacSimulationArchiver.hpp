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
     * @param simulationToArchive object defining the simulation to archive
     * @param directory directory where the multiple files defining the checkpoint will be stored
     *     (relative to CHASTE_TEST_OUTPUT)
     * @param clearDirectory whether the directory needs to be cleared or not.
     */
    static void Save(PROBLEM_CLASS& simulationToArchive, const std::string& rDirectory, bool clearDirectory=true);


    /**
     * Unarchives a simulation from the directory specified
     * 
     * @param directory directory where the multiple files defining the checkpoint are located
     *     (relative to CHASTE_TEST_OUTPUT)
     */
    static PROBLEM_CLASS* Load(const std::string& rDirectory);
    
    /**
     * Load a simulation, saved by any number of processes, into a single process.
     * 
     * Uses the DistributedVectorFactory saved in the process 0 archive to work out
     * how many secondary archive files to read, and loads the cells and boundary
     * conditions from these too.
     * 
     * Throws an exception if not called from a sequential simulation, as otherwise
     * we'd get deadlock from PetscTools::Barrier() calls in the archive loader. 
     * 
     * @param directory directory where the multiple files defining the checkpoint are located
     *     (relative to CHASTE_TEST_OUTPUT)
     */
    static PROBLEM_CLASS* LoadAsSequential(const std::string &rDirectory);
    
    /**
     * Load a simulation, saved by a single process, into any number of processes.
     * 
     * Uses a dumb partition to work out how to distribute the mesh and cells over
     * the processes.
     * \todo Allow the use of METIS to give a better partitioning
     * 
     * @param directory directory where the files defining the checkpoint are located
     *     (relative to CHASTE_TEST_OUTPUT)
     */
    static PROBLEM_CLASS* LoadFromSequential(const std::string& rDirectory);
};


template<class PROBLEM_CLASS>
void CardiacSimulationArchiver<PROBLEM_CLASS>::Save(PROBLEM_CLASS& simulationToArchive,
                                                    const std::string &rDirectory,
                                                    bool clearDirectory)
{
    // Clear directory if requested (and make sure it exists)
    OutputFileHandler handler(rDirectory, clearDirectory);
    
    // Open the archive files
    ArchiveOpener<boost::archive::text_oarchive, std::ofstream> archive_opener(rDirectory, rDirectory + ".arch", true);
    boost::archive::text_oarchive* p_main_archive = archive_opener.GetCommonArchive();
    
    // And save
    PROBLEM_CLASS* const p_simulation_to_archive = &simulationToArchive;
    (*p_main_archive) & p_simulation_to_archive;
}

template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::Load(const std::string &rDirectory)
{
    // Open the archive files
    ArchiveOpener<boost::archive::text_iarchive, std::ifstream> archive_opener(rDirectory, rDirectory + ".arch", true);
    boost::archive::text_iarchive* p_main_archive = archive_opener.GetCommonArchive();

    // Load
    PROBLEM_CLASS *p_unarchived_simulation;
    (*p_main_archive) >> p_unarchived_simulation;

    return p_unarchived_simulation;
}

template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::LoadAsSequential(const std::string &rDirectory)
{
    // Check that we're running sequentially
    if (!PetscTools::IsSequential())
    {
        EXCEPTION("Cannot load sequentially when running in parallel.");
    }
    PROBLEM_CLASS *p_unarchived_simulation;
    
    // Avoid the DistributedVectorFactory throwing a 'wrong number of processes' exception when loading
    DistributedVectorFactory::SetCheckNumberOfProcessesOnLoad(false);
    // Put what follows in a try-catch to make sure we reset this
    try
    {
        // Load the master and process-0 archive files.
        // This will also set up ArchiveLocationInfo for us.
        ArchiveOpener<boost::archive::text_iarchive, std::ifstream> archive_opener(rDirectory, rDirectory + ".arch", true);
        boost::archive::text_iarchive* p_main_archive = archive_opener.GetCommonArchive();
        (*p_main_archive) >> p_unarchived_simulation;
        
        // Work out how many more files to load
        DistributedVectorFactory* p_factory = p_unarchived_simulation->rGetMesh().GetDistributedVectorFactory();
        unsigned original_num_procs = p_factory->GetOriginalFactory()->GetNumProcs();
        
        // Merge in the extra data
        for (unsigned archive_num=1; archive_num<original_num_procs; archive_num++)
        {
            std::string archive_path = ArchiveLocationInfo::GetProcessUniqueFilePath(rDirectory + ".arch", archive_num);
            std::ifstream ifs(archive_path.c_str());
            boost::archive::text_iarchive archive(ifs);
            p_unarchived_simulation->LoadExtraArchive(archive);
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


template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::LoadFromSequential(const std::string &rDirectory)
{
    PROBLEM_CLASS *p_unarchived_simulation;
    
    // Avoid the DistributedVectorFactory throwing a 'wrong number of processes' exception when loading
    DistributedVectorFactory::SetCheckNumberOfProcessesOnLoad(false);
    // Put what follows in a try-catch to make sure we reset this
    try
    {
        // Load both the master and process-0 archive files.
        // The deserialize methods will be clever enough to throw away data that they don't need.
        ArchiveOpener<boost::archive::text_iarchive, std::ifstream> archive_opener(rDirectory, rDirectory + ".arch", true, 0u);
        boost::archive::text_iarchive* p_main_archive = archive_opener.GetCommonArchive();
        (*p_main_archive) >> p_unarchived_simulation;
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

