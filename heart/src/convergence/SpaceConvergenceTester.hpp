/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
