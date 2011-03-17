/*

Copyright (C) University of Oxford, 2005-2011

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

#include "CellMLToSharedLibraryConverter.hpp"

#include <sstream>
#include <unistd.h> // For getpid()
#include <sys/stat.h> // For mkdir()
#include <ctime>

#include "Exception.hpp"
#include "ChasteBuildRoot.hpp"
#include "PetscTools.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "GetCurrentWorkingDirectory.hpp"

CellMLToSharedLibraryConverter::CellMLToSharedLibraryConverter(bool preserveGeneratedSources,
                                                               std::string component)
    : mPreserveGeneratedSources(preserveGeneratedSources),
      mComponentName(component)
{
}

DynamicCellModelLoader* CellMLToSharedLibraryConverter::Convert(const FileFinder& rFilePath,
                                                                bool isCollective)
{
    DynamicCellModelLoader* p_loader;
    std::string absolute_path = rFilePath.GetAbsolutePath();
    // Check the file exists
    if (!rFilePath.Exists())
    {
        EXCEPTION("Dynamically loadable cell model '" + absolute_path + "' does not exist.");
    }
    // Find out whether rFilePath is a .cellml or .so
    size_t dot_position = absolute_path.find_last_of(".");
    if (dot_position == std::string::npos)
    {
        EXCEPTION("File does not have an extension: " + absolute_path);
    }
    std::string extension = absolute_path.substr(dot_position+1);
    if (extension == "cellml")
    {
        // Split the path into folder and leaf
        size_t slash_position = absolute_path.find_last_of("/\\");
        assert(slash_position != std::string::npos);
        std::string folder = absolute_path.substr(0, slash_position+1); // Include trailing slash
        std::string leaf = absolute_path.substr(slash_position+1, dot_position-slash_position); // Include dot
        std::string so_path = folder + "lib" + leaf + "so";
        // Does the .so file already exist (and was it modified after the .cellml?)
        FileFinder so_file(so_path, RelativeTo::Absolute);
        if (!so_file.Exists() || rFilePath.IsNewerThan(so_file))
        {
            if (!isCollective)
            {
                EXCEPTION("Unable to convert .cellml to .so unless called collectively, due to possible race conditions.");
            }
            ConvertCellmlToSo(absolute_path, folder, leaf);
        }
        // Load the .so
        p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(so_file);
    }
    else if (extension == "so")
    {
        // Just load the .so
        p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(rFilePath);
    }
    else
    {
        EXCEPTION("Unsupported extension '." + extension + "' of file '" + absolute_path + "'; must be .so or .cellml");
    }

    return p_loader;
}

void CellMLToSharedLibraryConverter::ConvertCellmlToSo(const std::string& rCellmlFullPath,
                                                       const std::string& rCellmlFolder,
                                                       const std::string& rModelLeafName)
{
    std::string tmp_folder, build_folder;
    std::string old_cwd = GetCurrentWorkingDirectory();
    try
    {
        // Need to create a .so file from the CellML...
        if (PetscTools::AmMaster())
        {
            // Create a temporary folder within heart/dynamic
            std::stringstream folder_name;
            folder_name << "dynamic/tmp_" << getpid() << "_" << time(NULL);
            tmp_folder = std::string(ChasteBuildRootDir()) + mComponentName + "/" + folder_name.str();
            build_folder = std::string(ChasteBuildRootDir()) + mComponentName + "/build/" + ChasteBuildDirName() + "/" + folder_name.str();
            int ret = mkdir(tmp_folder.c_str(), 0700);
            EXCEPT_IF_NOT(ret == 0);
            // Copy the .cellml file (and any relevant others) into the temporary folder
            size_t dot_pos = rCellmlFullPath.rfind('.');
            std::string cellml_base = rCellmlFullPath.substr(0, dot_pos);
            EXPECT0(system, "cp " + cellml_base + "* " + tmp_folder);
            // If there's a config file, copy that too
            std::string config_path = rCellmlFullPath.substr(0, rCellmlFullPath.length() - 7) + "-conf.xml";
            if (FileFinder(config_path, RelativeTo::Absolute).Exists())
            {
                EXPECT0(system, "cp " + config_path + " " + tmp_folder);
            }
            // Change to Chaste source folder
            EXPECT0(chdir, ChasteBuildRootDir());
            // Run scons to generate C++ code and compile it to a .so
            EXPECT0(system, "scons --warn=no-all dyn_libs_only=1 build=" + ChasteBuildType() + " " + tmp_folder);
            EXCEPT_IF_NOT(FileFinder(tmp_folder + "/lib" + rModelLeafName + "so", RelativeTo::Absolute).Exists());
            // CD back
            EXPECT0(chdir, old_cwd);
            // Copy the .so to the same folder as the original .cellml file
            EXPECT0(system, "cp " + tmp_folder + "/lib" + rModelLeafName + "so " + rCellmlFolder);
            if (mPreserveGeneratedSources)
            {
                // Copy generated source code as well
                EXPECT0(system, "cp " + build_folder + "/*.?pp " + rCellmlFolder);
            }
            // Delete the temporary folders
            EXPECT0(system, "rm -r " + build_folder);
            EXPECT0(system, "rm -r " + tmp_folder);
        }
    }
    catch (Exception& e)
    {
        PetscTools::ReplicateException(true);
        // Delete the temporary folders
        EXPECT0(system, "rm -rf " + build_folder); // -f because folder might not exist
        EXPECT0(system, "rm -r " + tmp_folder);
        EXPECT0(chdir, old_cwd);
        EXCEPTION("Conversion of CellML to Chaste shared object failed. Error was: " + e.GetMessage());
    }
    // This also has the effect of a barrier, ensuring all processes wait for the
    // shared library to be created.
    PetscTools::ReplicateException(false);
}
