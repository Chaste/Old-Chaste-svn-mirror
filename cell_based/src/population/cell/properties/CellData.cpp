/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "CellData.hpp"
#include "Debug.hpp"

void CellData::SetItem(const std::string& variableName, double data)
{
    mCellData[variableName] = data;
}

void CellData::SetItem(unsigned variableNumber, double data)
{
    switch (variableNumber)
    {
        case 0:
            SetItem("Var0", data);
            break;
        case 1:
            SetItem("Var1", data);
            break;
        case 2:
            SetItem("Var2", data);
            break;
        default:
            NEVER_REACHED;
    }
}

double CellData::GetItem(const std::string& variableName) const
{
    //Note that mCellData[variableName] is not const.
    //If variableName is not a key, then mCellData[variableName] will create a new item in the map
    //and increase the size by one.  Using a const_iterator ensures that the map remains const.
    std::map<std::string, double>::const_iterator it = mCellData.find(variableName);
    if (it == mCellData.end())
    {
        PRINT_VECTOR(GetKeys());///\todo This is useful for #2115
        EXCEPTION("The item " << variableName << " is not stored");
    }
    if (it->second == DOUBLE_UNSET)
    {
        EXCEPTION("The item " << variableName << " has not yet been set");
    }
    return(it->second);
}

double CellData::GetItem(unsigned variableNumber) const
{
    switch (variableNumber)
    {
        case 0:
            return GetItem("Var0");
            break;
        case 1:
            return GetItem("Var1");
            break;
        case 2:
            return GetItem("Var2");
            break;
        default:
            NEVER_REACHED;
    }

}

unsigned CellData::GetNumItems() const
{
    return mCellData.size();
}

std::vector<std::string> CellData::GetKeys() const
{
    std::vector<std::string> keys;
    //Note: Does the order of the keys matter?  If so, then sort the vector + Doxygen
    for (std::map<std::string, double>::const_iterator it = mCellData.begin(); it != mCellData.end(); ++it)
    {
        keys.push_back(it->first);
    }

    // We assume that the iterator is returning sorted keys sort(keys.begin(), keys.end());
    return keys;
}


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellData)
