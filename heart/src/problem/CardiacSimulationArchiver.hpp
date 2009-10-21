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

// Must be included before any other serialisation headers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "OutputFileHandler.hpp"
#include "ArchiveLocationInfo.hpp"

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
     *  Archives a simulation in the directory specified.
     * 
     *  @param simulationToArchive object defining the simulation to archive
     *  @param directory directory where the multiple files defining the checkpoint will be stored
     *  @param clearDirectory whether the directory needs to be cleared or not.
     */
    static void Save(PROBLEM_CLASS& simulationToArchive, std::string directory, bool clearDirectory=true);


    /**
     *  Unarchives a simulation from the directory specified
     * 
     *  @param directory directory where the multiple files defining the checkpoint are located
     */
    static PROBLEM_CLASS* Load(std::string directory);
    
     /**
     *  Archives a simulation in the directory specified.
     * 
     *  @param simulationToArchive object defining the simulation to archive
     *  @param directory directory where the multiple files defining the checkpoint will be stored
     *  @param clearDirectory whether the directory needs to be cleared or not.
     */
    static void SaveAsSequential(PROBLEM_CLASS& simulationToArchive, std::string directory, bool clearDirectory=true);
 
    /**
     * Take a parallel archive and convert it to a sequential one
     * @param inputDirectory
     * @param outputDirectory
     * @param clearDirectory whether the directory needs to be cleared or not.
     */
    static void MigrateToSequential(std::string inputDirectory, std::string outputDirectory, bool clearDirectory=true);
};


template<class PROBLEM_CLASS>
void CardiacSimulationArchiver<PROBLEM_CLASS>::Save(PROBLEM_CLASS& simulationToArchive, std::string directory, bool clearDirectory)
{
    OutputFileHandler handler(directory, clearDirectory);
    handler.SetArchiveDirectory();
    std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath(directory + ".arch");

    std::ofstream ofs(archive_filename.c_str());
    boost::archive::text_oarchive output_arch(ofs);
    PROBLEM_CLASS* const p_simulation_to_archive = &simulationToArchive;
    output_arch & p_simulation_to_archive;
}

template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::Load(std::string directory)
{
    OutputFileHandler handler(directory, false);
    handler.SetArchiveDirectory();
    std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath(directory + ".arch");

    std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
    if(!ifs.is_open())
    {
        EXCEPTION("Cannot load file: " + archive_filename);
    }
    boost::archive::text_iarchive input_arch(ifs);

    PROBLEM_CLASS *p_unarchived_simulation;
    input_arch >> p_unarchived_simulation;

    return p_unarchived_simulation;
}

template<class PROBLEM_CLASS>
void CardiacSimulationArchiver<PROBLEM_CLASS>::SaveAsSequential(PROBLEM_CLASS& simulationToArchive, std::string directory, bool clearDirectory)
{
    OutputFileHandler handler(directory, clearDirectory);
    handler.SetArchiveDirectory();
    std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath(directory + ".arch");

    std::ofstream ofs(archive_filename.c_str());
    boost::archive::text_oarchive output_arch(ofs);
    PROBLEM_CLASS* const p_simulation_to_archive = &simulationToArchive;
   ///\todo #1159 
    output_arch & p_simulation_to_archive;
}

template<class PROBLEM_CLASS>
void CardiacSimulationArchiver<PROBLEM_CLASS>::MigrateToSequential(std::string inputDirectory, std::string outputDirectory, bool clearDirectory)
{
    OutputFileHandler handler(inputDirectory, false);
    handler.SetArchiveDirectory();
    std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath(inputDirectory + ".arch");
    
    //Assume that the archive has been written for the right number of processors
    if (PetscTools::IsSequential())
    {
        EXCEPTION("Archive doesn't need to be migrated since it is already sequential");
    }
    std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
    assert(ifs.is_open());
    boost::archive::text_iarchive input_arch(ifs);

    PROBLEM_CLASS *p_unarchived_simulation;
    input_arch >> p_unarchived_simulation;
    
    SaveAsSequential(*p_unarchived_simulation, outputDirectory, clearDirectory);
    delete p_unarchived_simulation;
}
#endif /*CARDIACSIMULATIONARCHIVER_HPP_*/
