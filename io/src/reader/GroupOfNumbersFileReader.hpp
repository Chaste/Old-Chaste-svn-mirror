#ifndef GROUPOFNUMBERSFILEREADER_HPP_
#define GROUPOFNUMBERSFILEREADER_HPP_

#include <string>
#include <vector>
#include <fstream>
#include "Exception.hpp"

template <class T>
class GroupOfNumbersFileReader
{
private :
    std::vector<T> mData;
    
public :
    GroupOfNumbersFileReader(std::string fileName)
    {
        // open file
        std::ifstream file(fileName.c_str(), std::ios::in);
        // throw exception if file can't be opened
        if (!file.is_open())
        {
            throw Exception("GroupOfNumbersFileReader.hpp: Couldn't open file " + fileName);
        }
        
        // read data
        while(!file.eof())
        {
            T value;
            file >> value;
            //std::cout << file.fail() << " " << file.eof() << " " << value << "\n";
            if(!file.fail())
            {
                mData.push_back(value);
            }
            else if (!file.eof())
            {
                // failed but no because reached end of file => error
                throw Exception("GroupOfNumbersFileReader.hpp, error reading data from " 
                                 + fileName);
            }
        }
    }
    
    std::vector<T> GetData()
    {
        return mData;
    }
};
#endif /*GROUPOFNUMBERSFILEREADER_HPP_*/
