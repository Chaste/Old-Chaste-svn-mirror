#ifndef LOGFILE_HPP_
#define LOGFILE_HPP_

#include <string>
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include <cassert>
//#include <iostream>

/**
 *  A singleton log file class. Allows the user to define log file in the test, which
 *  can be written to (without being passed around in the code).
 * 
 *  Usage (in test):
 *  LogFile* p_log = LogFile::Instance();
 *  p_log->SetDirectoryAndFile("dir","file");
 *  // run simulatiom
 *  LogFile::Close();
 * 
 *  Usage (in source):
 *  (*LogFile::Instance()) << "Info to be written to the log file\n";
 * 
 *  If no log file is set up,  
 *  (*LogFile::Instance()) << "Info to be written to the log file\n";
 *  does nothing.
 */
class LogFile
{
private:
    /** the static single instance */
    static LogFile* mpInstance;
    /** whether a directory and filename has been set */
    bool mFileSet;
    /** the file to be written to */
    out_stream mpOutStream;
    
    unsigned mLevel;

    static const unsigned mMaxLoggingLevel = 2;

public:
    /**
     *  Get the single instance of the LogFile object. 
     */
    static LogFile* Instance()
    {
        if (mpInstance == NULL)
        {
            mpInstance = new LogFile; // default construtor which doesn't write
        }
        return mpInstance;
    }
    
    static unsigned Level()
    {
        if (mpInstance == NULL)
        {
            return 0;
        }
        else
        {
            return mpInstance->mLevel;
        }
    }

    /**
     *  Constructor. Should never be called directly, call LogFile::Instance() instead.
     */
    LogFile()
    {
        mFileSet = false;
        mLevel = 0;
    }

    /**
     *  Set the logging level, the directory (relative to TEST_OUTPUT) and the file 
     *  the log should be written to (file defaults to "log.txt").
     *  
     *  The level should be a number between 0 and LogFile
     * 
     *  Note: we intentionally do NOT check or throw an exception if a file has already
     *  been set (ie Close() wasn't called the last time a log was used).
     *  
     *  The directory is never cleaned.
     */
    void Set(unsigned level, std::string directory, std::string fileName="log.txt")
    {
        if(level > mMaxLoggingLevel)
        {
            std::stringstream string_stream;
            string_stream << "Requested level " << Level 
                          << " should have been less than or equal to " << mMaxLoggingLevel;
            EXCEPTION(string_stream.str());
        }
        mLevel = level;

        OutputFileHandler handler(directory, false);
        mpOutStream = handler.OpenOutputFile(fileName);
        mFileSet = true;
        
        // write header in the log file..?
    }
    
    static unsigned MaxLoggingLevel()
    {
        return mMaxLoggingLevel;
    } 

    /**
     *  Close the currently open file, and delete the single LogFile instance
     */
    static void Close()
    {
        if (mpInstance)
        {
            mpInstance->mpOutStream->close();
            delete mpInstance;
            mpInstance = NULL;
        }
    }
    
    /**
     *  Whether SetDirectoryAndFile() or SetDirectory() has been called
     */
    bool IsFileSet()
    {
        return mFileSet;
    }
    
   
    /**
     *  Overloaded << operator, to write to the log file, if one has been set, and
     *  does nothing if not
     */
    template <class T>
    LogFile& operator<<(T message)
    {
        if(mFileSet)
        {
            (*mpOutStream) << message << std::flush;
        }

        return *this;
    }
};

LogFile* LogFile::mpInstance = NULL;

#define LOG(level, message) assert(level>0); if(level <= LogFile::Level()) { (*LogFile::Instance()) << message; }

#endif /*LOGFILE_HPP_*/
