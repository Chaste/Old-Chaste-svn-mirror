#ifndef INPUTFILEREADER_HPP_
#define INPUTFILEREADER_HPP_

#include <string>
#include <vector>
#include <fstream>

#include "Exception.hpp"


class InputFileReader
{
private:
    std::string mFileName;  

    std::string RemoveLeadingWhitespace(std::string st)
    {
        std::string::size_type index = st.find_first_not_of(" \t\n"); // find first char not equal to ' ' or '\t'
        if( index != std::string::npos )
        {
            return st.substr(index);
        }
        else
        {
            return st;
        }
    }
         
    std::string RemoveTrailingWhitespace(std::string st)
    {
        std::string::size_type index = st.find_last_not_of(" \t\n"); // find first char not equal to ' ' or '\t'
        if( index != std::string::npos )
        {
            return st.substr(0,index+1);
        }
        else
        {
            return st;
        }
    }

    
    /** 
     *  Find a string in the file.
     * 
     *  This method searches for an occurance of the given field name in the file. The field name
     *  has to be at the beginning of a line, excluding whitespace, and must be followed by
     *  a colon  (or whitespace and then a colon). Exceptions are thrown if either of these 
     *  are not true, or if the field is found twice. If the field name is not found, 
     *  an exception is only thrown if the final argument is passed as true.
     * 
     *  @param the name of the field being searched for. 
     *  @param a boolean which will be set to true if the field was found, and false otherwise
     *  @param a boolean saying whether to quit (throw Exception) if the string is not found
     * 
     *  @return everything else in the line in the file containing that string, excluding the colon 
     *  and anything following a hash (taken to represent a comment) 
     */
    std::string Find(std::string fieldName, bool& foundFieldName, bool quitIfNotFound)
    {
        if(fieldName.size()==0)
        {
            throw Exception("InputFileReader.hpp: error, asked to find empty string");
        }
 
        std::ifstream file(mFileName.c_str(), std::ios::in);
        // throw exception if file can't be opened
        if (!file.is_open())
        {
            throw Exception("InputFileReader.hpp: Couldn't open file " + mFileName);
        }

        foundFieldName = false;
        std::string line;
        std::string return_string;

        bool end_of_file_reached=false;
        while(!end_of_file_reached)
        {
            // get the next line and check if end of file has been reached
            end_of_file_reached = std::getline(file, line).eof();

            // remove whitespace at beginning of the line
            line = RemoveLeadingWhitespace(line);
            
            // remove anything after a hash (taken as a comment)
            std::string::size_type index = line.find_first_of('#');
            if(index != std::string::npos)
            {
                line = line.substr(0,index);
            }
            
            // is the beginning of the line the same as the name being search for?
            if(line.substr(0, fieldName.length())==fieldName)
            {
                if(foundFieldName)
                {
                    // appears we have found the name twice, throw Exception
                    file.close();
                    throw Exception("InputFileReader.hpp: error, found '" + fieldName +
                                    "'twice in file " + mFileName);
                }
                    
                // prepare return string
                return_string = line.substr(fieldName.length(), line.length());

                // remove leading whitespace from returned string 
                return_string = RemoveLeadingWhitespace(return_string);
                    
                // check next character is a colon
                if(return_string.substr(0,1)!=":")
                {
                    file.close();
                    throw Exception("InputFileReader.hpp: error, found '"+fieldName+"' in file "
                                    +mFileName+", but it was not followed by a colon");
                }
                
                // remove the colon, then ready to return
                return_string = return_string.substr(1,return_string.length());
                foundFieldName = true;                  
            }
            // name is not at beginning, check the name isn't elsewhere in the line
            // before we move to the next line
            else if(line.find(fieldName)!=std::string::npos)
            {
                file.close();
                if(foundFieldName)
                {
                    // appears we have found the name twice, throw Exception
                    throw Exception("InputFileReader.hpp: error, found '" + fieldName +
                                    "'twice in file " + mFileName);
                }

                throw Exception("InputFileReader.hpp: error, found '"+fieldName+"' in file "
                                    +mFileName+", but it was not at the beginning of a line, "
                                    +"or commented out");
            }
        }
        file.close();
        
        if( (!foundFieldName) && (quitIfNotFound) )
        {
            // unable to find name, and have been asked to quit if so, so throw
            throw Exception("InputFileReader.hpp: error, unable to find '"
                            + fieldName + "' in file " + mFileName);
        }
        return return_string;  // will be the empty string if found = false          
    }

    
public :
    /** 
     *  Constructor
     * 
     *  @param The name of the file to be read. 
     */
    InputFileReader(std::string fileName)
     : mFileName(fileName) // save the file name
    {
    }


    /** 
     *  Read a vector of variables from the file.
     * 
     *  This method searches for an occurance of the given field name in the file. 
     *  The field name can be any string, but must be at the beginning of a line, excluding whitespace, 
     *  and must be followed by a colon (or whitespace and then a colon). Exceptions are 
     *  thrown if either of these are not true, or if the field is found twice. The method reads 
     *  and returns the data after the field name and the colon.
     *  
     *  @param the name of the field being searched for.
     *  @param the number of data entries after the field name to be read  
     *  @param a boolean which will be set to true if the field was found, and false otherwise
     *  @param optional boolean saying whether to quit (throw Exception) if the string is not found. Defaults to true.
     * 
     *  @return std::vector of the type specified containing the data read. If the field was not found but
     *   quitIfNotFound was passed in as false, the empty vector is returned.
     * 
     *  Method should be called using one of
     *    data = ReadVector<int>   ("MyField", 4, found)
     *    data = ReadVector<double>("MyField", 4, found)
     *    data = ReadVector<string>("MyField", 4, found)
     */
    template <class T>
    std::vector<T> ReadVector(std::string fieldName, unsigned numberOfEntries, bool& found, bool quitIfNotFound = true)
    {
        found = false;
        // the returned vector
        std::vector<T> ret;

        // search for <fieldName> in the file. found will be set to be true if it is
        // found and data will be the rest of the appropriate line in the file
        std::string data = Find(fieldName, found, quitIfNotFound);
        if(found)
        {
            // get data from the string
            std::stringstream line_stream(data);
            for(unsigned i=0; i<numberOfEntries; i++)
            {
                T value;
                line_stream >> value;
                if(line_stream.fail())
                {
                    std::stringstream ss;
                    ss << "InputFileReader.hpp, error while attempting to read " << numberOfEntries
                       << " data entries corresponding to '" << fieldName << "' in file "
                       << mFileName;
                    throw Exception(ss.str());
                }
                ret.push_back(value);
            }
        }
        return ret;  // empty vector if the name wasn't found in the file
    }
    
        
    /** 
     *  Read a string from the file.
     * 
     *  This method searches for an occurance of the given field name in the file, and reads
     *  the data corresponding to that field.
     *  The field name can be any string, but must be at the beginning of a line, excluding whitespace, 
     *  and must be followed by a colon (or whitespace and then a colon). Exceptions are 
     *  thrown if either of these are not true, or if the field is found twice. The method reads 
     *  and returns the data after the field name and the colon.
     *  
     *  @param the name of the field being searched for.
     *  @param a boolean which will be set to true if the field was found, and false otherwise
     *  @param optional boolean saying whether to quit (throw Exception) if the string is not found. Defaults to true.
     * 
     *  @return string containing the data read. If the field was not found but
     *   quitIfNotFound was passed in as false, the empty string is returned.
     */        
    std::string ReadString(std::string fieldName, bool& found, bool quitIfNotFound = true)
    {
        std::vector<std::string> ret = ReadVector<std::string>(fieldName, 1, found, quitIfNotFound);
        if(found)
        {
            return ret[0];
        }
        else
        {
            return "";
        }
    }

    /** 
     *  Exactly the same as ReadString but reads an integer. If the field was not found but
     *  quitIfNotFound was passed in as false, 0 is returned.
     */        
    int ReadInt(std::string fieldName, bool& found, bool quitIfNotFound = true)
    {
        std::vector<int> ret = ReadVector<int>(fieldName, 1, found, quitIfNotFound);
        if(found)
        {
            return ret[0];
        }
        else
        {
            return 0;
        }
    }

    /** 
     *  Exactly the same as ReadString but reads an double. If the field was not found but
     *  quitIfNotFound was passed in as false, 0 is returned.
     */        
     double ReadDouble(std::string fieldName, bool& found, bool quitIfNotFound = true)
     {
        std::vector<double> ret = ReadVector<double>(fieldName, 1, found, quitIfNotFound);
        if(found)
        {
            return ret[0];
        }
        else
        {
            return 0;
        }
    }
};

#endif /*INPUTFILEREADER_HPP_*/
