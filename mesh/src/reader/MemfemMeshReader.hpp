/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _MEMFEMMESHREADER_HPP_
#define _MEMFEMMESHREADER_HPP_

/**
 * Concrete version of the AbstractMeshReader class.
 * A MemfemMeshReader takes the base name of a set of Memfem
 * mesh files (ie. the path and name of the files without the suffices).
 * Once constructed the public methods of the AbstractMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */

#include "AbstractMeshReader.cpp"
#include "Exception.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MemfemMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:
    std::vector<std::vector<double> > TokenizeStringsToDoubles(
        std::vector<std::string> rawData);
        
    std::vector<std::vector<unsigned> > TokenizeStringsToInts(
        std::vector<std::string> rawData,
        unsigned dimensionOfObject,
        bool readHeader);
        
public:
    MemfemMeshReader(std::string pathBaseName);
    virtual ~MemfemMeshReader();
};

#endif //_MEMFEMMESHREADER_HPP_
