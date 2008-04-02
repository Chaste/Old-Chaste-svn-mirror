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

/**
 * Concrete version of the AbstractMeshReader class.
 * A TrianglesMeshReader takes the base name of a set of Triangles or
 * Tetgen mesh files (ie. the path and name of the files without the suffices).
 * Once constructed the public methods of the AbstractMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */
#ifndef _TRIANGLESMESHREADER_H_
#define _TRIANGLESMESHREADER_H_

#include "AbstractMeshReader.cpp"
#include "Exception.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:
    std::vector<std::vector<double> > TokenizeStringsToDoubles(
        std::vector<std::string> rawData);
        
    /**
     *  Convert strings for vectors of unsigned
     *  @param rawData the raw data. Each string should look like: index value0 value1 .. valueN marker
     *  Here marker doesn't have to be present, it is ignored unless onlyMarked=true
     *  @dimensionOfObject The number of values
     *  @onlyMarked Set this to true to look at the marker and ignore any strings 
     *  for which the marker is set to zero.
     */
    std::vector<std::vector<unsigned> > TokenizeStringsToInts(
        std::vector<std::string> rawData,
        unsigned dimensionOfObject,
        bool onlyMarked=false);

    void ReadFacesAsElements(std::string pathBaseName);
    void ReadEdgesAsFaces(std::string pathBaseName);
    
public:
    TrianglesMeshReader(std::string pathBaseName);
    
    virtual ~TrianglesMeshReader();
};

#endif //_TRIANGLESMESHREADER_H_
