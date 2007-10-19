#ifndef LOGFILE_HPP_
#define LOGFILE_HPP_

#include <string>
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include <cassert>

/**
 *  A singleton log file class. Allows the user to define log file in the test, which
 *  can be written to (without being passed around in the code).
 * 
 *  Usage (in test):
 *  
 *  // begining of test
 *  LogFile* p_log = LogFile::Instance();
 *  p_log->Set(level, "dir","file");
 *  p_log->WriteHeader("type_of_sim"); // optional
 * 
 *  // at end of simulation
 *  LogFile::Close();
 * 
 *  Here 'level' is a number between 0 and LogFile::MaxLoggingLevel, with zero
 *  meaning no logging and MaxLoggingLevel meaning full logging.
 * 
 *  Usage (in source) - use the macro 'LOG'
 *  LOG(1, "Info to be written to the log file\n" << "More info\n");
 *  LogFile::Instance()->WriteElapsedTime(); // optional
 * 
 *  This checks whether the given level (here '1') is greater or equal to the given
 *  logging level, in which case it writes to the current log file. If there is
 *  no log file is set up it does nothing.
 * 
 *  Note the log file can be written to directly, without any level-checking, using  
 *  (*LogFile::Instance()) << "Info to be written to the log file\n";
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

    time_t mInitTime;

    /** the level of logging required for this particular log file */
    unsigned mLevel;

    /** the max level of logging */
    static const unsigned mMaxLoggingLevel = 2;

    /**
     *  Constructor. Should never be called directly, call LogFile::Instance() instead.
     */
    LogFile();

public:
    /**
     *  Get the single instance of the LogFile object. 
     */
    static LogFile* Instance();
    
    static unsigned Level();

    /**
     *  Set the logging level, the directory (relative to TEST_OUTPUT) and the file 
     *  the log should be written to (file defaults to "log.txt").
     *  
     *  The level should be a number between 0 and LogFile::MaxLoggingLevel() (which is the 
     *  same as LogFile::mMaxLoggingLevel)
     * 
     *  Note: we intentionally do NOT check or throw an exception if a file has already
     *  been set (ie Close() wasn't called the last time a log was used).
     *  
     *  The directory is never cleaned.
     */
    void Set(unsigned level, std::string directory, std::string fileName="log.txt");

    /** Get the maximum allowed logging level */
    static unsigned MaxLoggingLevel();

    /** 
     *  Write a header in the log file, stating the (given) type of simulation and the 
     *  date and time
     *  
     *  @simulationType The type of simulation, eg "Bidomain" or "Crypt" or 
     *  "Cardiac Electromechanics". Defaults to empty
     */
    void WriteHeader(std::string simulationType="");
    
    /** 
     *  Write the elapsed time since the simulation began (since the log file was created) 
     *  @param pre a string (eg spacings) to write before the elapsed time line
     */
    void WriteElapsedTime(std::string pre="");

    /**
     *  Close the currently open file, and delete the single LogFile instance
     */
    static void Close();

    /**
     *  Whether Set() has been called
     */
    bool IsFileSet();    
   
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

#ifndef NDEBUG
    // define the log macro
    #define LOG(level, message) assert(level>0); if(level <= LogFile::Level()) { (*LogFile::Instance()) << message << "\n"; }
#else
    // do nothing
    #define LOG(level, message) 
#endif


#endif /*LOGFILE_HPP_*/
