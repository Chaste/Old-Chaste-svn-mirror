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

#ifndef _MESHALYZERMESHWRITER_HPP_
#define _MESHALYZERMESHWRITER_HPP_

#include "AbstractMeshWriter.cpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MeshalyzerMeshWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
public:
    MeshalyzerMeshWriter(const std::string &rDirectory,
                         const std::string &rBaseName,
                         const bool &rCleanDirectory=true,
                         const bool &rSetCoolGraphics=false);
    void WriteFiles();
    virtual ~MeshalyzerMeshWriter();
    
};

#endif //_MESHALYZERMESHWRITER_HPP_
