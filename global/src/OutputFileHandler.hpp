#ifndef OUTPUTFILEHANDLER_HPP_
#define OUTPUTFILEHANDLER_HPP_

#include <string>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <sys/stat.h>

#include "Exception.hpp"
#include <petsc.h>

typedef std::auto_ptr<std::ofstream> out_stream;

/**
 * This file abstracts stuff that needs to be done when creating output files for tests.
 * It defines some helpful functions, so that things that are otherwise repeated in lots
 * of places are just done here.
 */
class OutputFileHandler
{
private:
    std::string mDirectory; ///< The directory to store output files in
    bool mAmMaster; ///< Are we the master process?
    
public:
    /**
     * Create an OutputFileHandler that will create output files in the given directory.
     * The directory name should be relative to the place where Chaste test output is
     * stored.  If the user needs to know where this is, use the GetTestOutputDirectory
     * method.
     * 
     * Will check that the directory exists and create it if needed.
     * 
     * @param rDirectory  the directory to put output files in.
     * @param rCleanOutputDirectory  whether to remove any existing files in the output directory
     */
    OutputFileHandler(const std::string &rDirectory,
                      bool rCleanOutputDirectory = true)
    {
        // Are we the master process?  Only the master should do any writing to disk
        PetscTruth is_there;
        PetscInitialized(&is_there);
        if (is_there == PETSC_TRUE)
        {
            int my_rank;
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
        {
            // Not using PETSc, so we're definitely the only process
            mAmMaster = true;
        }
        mDirectory = GetTestOutputDirectory(rDirectory);
        
        // Clean the output dir?
        if (rCleanOutputDirectory && mAmMaster &&
            rDirectory != "" && rDirectory.find("..") == std::string::npos)
        {
            // Remove the directory itself rather than contents, to avoid
            // problems with too long command lines
            system(("rm -rf " + mDirectory).c_str());
            // Re-create the directory
            mkdir(mDirectory.c_str(), 0775);
        }
    }
    
    /**
     * Check that the desired output directory exists and is writable by us.
     * Create it if needed.
     * Return the full pathname of the output directory.
     * 
     * The environment variable CHASTE_TEST_OUTPUT will be examined.  If it is set
     * and non-empty it is taken to be a directory where test output should be stored.
     * Otherwise a default location of testoutput/ within the current directory is used.
     * 
     * @param directory  pathname of the output directory, relative to where Chaste
     *         output will be stored (user shouldn't care about this).
     * @return  full pathname to the output directory
     */
    std::string GetTestOutputDirectory(std::string directory)
    {
        char *chaste_test_output = getenv("CHASTE_TEST_OUTPUT");
        std::string directory_root;
        if (chaste_test_output == NULL || *chaste_test_output == 0)
        {
            // Default to within the Chaste directory
            directory_root = "testoutput/";
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
    /**
     * Return the full pathname to the directory this object will create files
     * in.
     */
    std::string GetTestOutputDirectory()
    {
        return mDirectory;
    }
    
    /**
     * Open an output file in our directory, and check it was opened successfully.
     * Throws an Exception if not.
     * 
     * @param filename  the name of the file to open, relative to the output directory.
     * @param mode  optionally, flags to use when opening the file (defaults are as for
     *         std::ofstream).
     * @return  a managed pointer to the opened file stream.
     */
    out_stream OpenOutputFile(std::string filename,
                              std::ios_base::openmode mode=std::ios::out | std::ios::trunc)
    {
        out_stream p_output_file(new std::ofstream((mDirectory+filename).c_str(), mode));
        if (!p_output_file->is_open())
        {
            EXCEPTION("Could not open file " + filename + " in " + mDirectory);
        }
        return p_output_file;
    }
};

#endif /*OUTPUTFILEHANDLER_HPP_*/
