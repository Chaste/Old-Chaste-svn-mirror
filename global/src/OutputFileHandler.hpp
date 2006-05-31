#ifndef OUTPUTFILEHANDLER_HPP_
#define OUTPUTFILEHANDLER_HPP_

#include <string>
#include <fstream>
#include <cstdlib>

#include "Exception.hpp"

/**
 * This file abstracts stuff that needs to be done when creating output files for tests.
 * It defines some helpful functions, so that things that are otherwise repeated in lots
 * of places are just done here.
 */
class OutputFileHandler
{
private:
	std::string mDirectory; ///< The directory to store output files in
	
public:
	/**
     * Create an OutputFileHandler that will create output files in the given directory.
     * The directory name should be relative to the place where Chaste test output is
     * stored.  If the user needs to know where this is, use the GetTestOutputDirectory
     * method.
     * 
     * Will check that the directory exists and create it if needed.
     * 
     * @param directory  the directory to put output files in.
     */
    OutputFileHandler(std::string directory)
    {
        mDirectory = GetTestOutputDirectory(directory);
        
        // Add a trailing slash if not already there
        if (! ( *(mDirectory.end()-1) == '/'))
        {
            mDirectory = mDirectory + "/";
        }
    }
    
    /**
     * Check that the desired output directory exists and is writable by us.
     * Create it if needed.
     * Return the full pathname of the output directory.
     * 
     * @param directory  pathname of the output directory, relative to where Chaste
     *         output will be stored (user shouldn't care about this)
     * @return  full pathname to the output directory
     */
    std::string GetTestOutputDirectory(std::string directory)
    {
        // Find the current user's name
        std::string username = std::string(getenv("USER"));
        // Use it to have a separate output dir for each user
        directory = "/tmp/" + username + "/testoutput/" + directory;
        // Make sure it exists (ish)
        system(("mkdir -p " + directory).c_str());
        return directory;
    }
    
    /**
     * Open an output file in our directory, and check it was opened successfully.
     * Throws an Exception if not.
     * 
     * @param filename  the name of the file to open, relative to the output directory.
     * @param mode  optionally, flags to use when opening the file (defaults are as for
     *         std::ofstream).
     * @return  a pointer to the opened file stream.
     */
    std::ofstream* OpenOutputFile(std::string filename,
                                  std::ios_base::openmode mode=std::ios::out | std::ios::trunc)
    {
        std::ofstream *p_output_file = new std::ofstream((mDirectory+filename).c_str(), mode);
        if (!p_output_file->is_open())
        {
            delete p_output_file;
            throw Exception("Could not open file " + filename + " in " +
                            mDirectory);
        }
        return p_output_file;
    }
};

#endif /*OUTPUTFILEHANDLER_HPP_*/
