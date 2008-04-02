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
 * A FemlabMeshReader takes the file names of a set of Femlab mesh files.
 * Once constructed the public methods of the AbstractMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */
#ifndef _FEMLABMESHREADER_H_
#define _FEMLABMESHREADER_H_

#include "AbstractMeshReader.cpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class FemlabMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>

{
private:
    std::vector<std::vector<double> > TokenizeStringsToDoubles(
        std::vector<std::string> rawData);
        
    std::vector<std::vector<unsigned> > TokenizeStringsToInts(
        std::vector<std::string> rawData,
        unsigned dimensionOfObject);
public:
    FemlabMeshReader(std::string pathBaseName, std::string nodeFileName, std::string elementFileName, std::string edgeFileName);
    virtual ~FemlabMeshReader();
};

#endif //_FEMLABMESHREADER_H_
