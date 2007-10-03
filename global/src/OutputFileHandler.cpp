#include "OutputFileHandler.hpp"

#include <cstdlib>
#include <sys/stat.h>

#ifndef SPECIAL_SERIAL
#include <petsc.h>
#endif //SPECIAL_SERIAL
#include "Exception.hpp"


OutputFileHandler::OutputFileHandler(const std::string &rDirectory,
                                     bool rCleanOutputDirectory)
{
    // Are we the master process?  Only the master should do any writing to disk
#ifndef SPECIAL_SERIAL
    PetscTruth is_there;
    PetscInitialized(&is_there);
    if (is_there == PETSC_TRUE)
    {
        PetscInt my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0)
        {
            mAmMaster=true;
        }
        else
        {
            mAmMaster=false;
        }
    }
    else
#endif //SPECIAL_SERIAL
    {
        // Not using PETSc, so we're definitely the only process
        mAmMaster = true;
    }
    mDirectory = GetTestOutputDirectory(rDirectory);
    
    // Clean the output dir?
    if (rCleanOutputDirectory && mAmMaster &&
        rDirectory != "" && rDirectory.find("..") == std::string::npos)
    {
        std::string directory_to_move_to = GetTestOutputDirectory("last_cleaned_directory");
        system(("rm -rf " + directory_to_move_to).c_str());
        system(("mv " + mDirectory + " " + directory_to_move_to).c_str());
        //system(("rm -rf " + mDirectory).c_str());
        // Re-create the directory
        mkdir(mDirectory.c_str(), 0775);
    }
}


std::string OutputFileHandler::GetTestOutputDirectory(std::string directory)
{
    char *chaste_test_output = getenv("CHASTE_TEST_OUTPUT");
    std::string directory_root;
    if (chaste_test_output == NULL || *chaste_test_output == 0)
    {
        // Default to within the current directory
        directory_root = "./";
    }
    else
    {
        directory_root = std::string(chaste_test_output);
        // Add a trailing slash if not already there
        if (! ( *(directory_root.end()-1) == '/'))
        {
            directory_root = directory_root + "/";
        }
    }
    directory = directory_root + directory;
    // Make sure it exists (ish)
    if (mAmMaster)
    {
        system(("mkdir -p " + directory).c_str());
    }
    
    // Add a trailing slash if not already there
    if (! ( *(directory.end()-1) == '/'))
    {
        directory = directory + "/";
    }
    return directory;
}


std::string OutputFileHandler::GetTestOutputDirectory()
{
    return mDirectory;
}


out_stream OutputFileHandler::OpenOutputFile(std::string filename,
                                             std::ios_base::openmode mode)
{
    out_stream p_output_file(new std::ofstream((mDirectory+filename).c_str(), mode));
    if (!p_output_file->is_open())
    {
        EXCEPTION("Could not open file " + filename + " in " + mDirectory);
    }
    return p_output_file;
}


out_stream OutputFileHandler::OpenOutputFile(std::string fileName,
                                             unsigned number,
                                             std::string fileFormat,
                                             std::ios_base::openmode mode)
{
    std::stringstream string_stream;
    string_stream << fileName << number << fileFormat;
    return OpenOutputFile(string_stream.str(), mode);
}                                              
