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

#ifndef SPACECONVERGENCETESTER_HPP_
#define SPACECONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class SpaceConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM, PROBLEM_DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->MeshNum=0;
    }
    void UpdateConvergenceParameters()
    {
        this->MeshNum++;
    
    }
    bool GiveUpConvergence()
    {
        switch(DIM)
        {
            case 1:
            {
                return this->MeshNum>6;
                break;
            }
            case 2:
            {
                return this->MeshNum>6;
                break;
            }
            case 3:
            {
                return this->MeshNum>4;
                break;
            }
            default:
                assert(0);
                return true;//To keep Inntel compiler happy
        }
        return true;//To keep Intel compiler happy
    }
    double Abscissa()
    {
        unsigned mesh_size = (unsigned) pow(2, this->MeshNum+2); // number of elements in each dimension
        return this->mMeshWidth/(double) mesh_size;
    }
    
    int GetMeshNum()
    {
        return (int) this->MeshNum; //unsigned -> int is just cosmetic here.  (The test looks prettier).
    }
    double GetSpaceStep()
    {
        unsigned mesh_size = (unsigned) pow(2, this->MeshNum+2);// number of elements in each dimension
        double scaling = this->mMeshWidth/(double) mesh_size;
        return scaling;
    }
    
};

#endif /*SPACECONVERGENCETESTER_HPP_*/
