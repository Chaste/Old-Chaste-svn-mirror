#include "LogFile.hpp"

LogFile* LogFile::mpInstance = NULL;

LogFile::LogFile()
{
    mFileSet = false;
    mLevel = 0;
}


LogFile* LogFile::Instance()
{
    if (mpInstance == NULL)
    {
            mpInstance = new LogFile; // default construtor which doesn't write
    }
    return mpInstance;
}

unsigned LogFile::Level()
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


void LogFile::Set(unsigned level, std::string directory, std::string fileName)
{
    if(level > mMaxLoggingLevel)
    {
        std::stringstream string_stream;
        string_stream << "Requested level " << level 
                      << " should have been less than or equal to " << mMaxLoggingLevel;
        EXCEPTION(string_stream.str());
    }
    mLevel = level;

    OutputFileHandler handler(directory, false);
    mpOutStream = handler.OpenOutputFile(fileName);
    mFileSet = true;
    
    // write header in the log file..?
}

unsigned LogFile::MaxLoggingLevel()
{
    return mMaxLoggingLevel;
} 


void LogFile::Close()
{
    if (mpInstance)
    {
        mpInstance->mpOutStream->close();
        delete mpInstance;
        mpInstance = NULL;
    }
}


bool LogFile::IsFileSet()
{
    return mFileSet;
}
