/*

Copyright (C) University of Oxford, 2008

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


#ifndef ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_
#define ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_

#include <vector>
#include "AbstractDealiiAssembler.hpp" // just to include some dealii classes

/**
 *  Defines an interface and provides quad point information, for use
 *  by AbstractCardiacElectroMechanicsProblem
 */
template<unsigned DIM>
class AbstractCardiacMechanicsAssembler
{
protected :
    static const unsigned mNumQuadPointsInEachDimension = 3;

    unsigned mTotalQuadPoints;
    unsigned mCurrentQuadPointGlobalIndex;

    /**
     *  The x stretch, quad point-wise. NOTE: the i-th entry of this vector is
     *  assumed to be the i-th quad point obtained by looping over cells in the obvious
     *  way and then looping over quad points
     */
    std::vector<double> mLambda;

public :
    AbstractCardiacMechanicsAssembler(Triangulation<DIM>* pMesh)
    {
        assert(pMesh!=NULL);

        // set up quad point info
        QGauss<DIM>   quadrature_formula(mNumQuadPointsInEachDimension);
        mTotalQuadPoints = quadrature_formula.n_quadrature_points * pMesh->n_active_cells();
        mCurrentQuadPointGlobalIndex = 0;

        mLambda.resize(mTotalQuadPoints, 1.0);
    }

    virtual ~AbstractCardiacMechanicsAssembler()
    {
    }


    unsigned GetNumQuadPointsInEachDimension()
    {
        return mNumQuadPointsInEachDimension;
    }

    /**
     *  Get the total number of quadrature points (equal to the n.q^d, where n=number of cells
     *  q is the number of quad points in each dimension, and d=DIM).
     */
    unsigned GetTotalNumQuadPoints()
    {
        return mTotalQuadPoints;
    }

    /**
     *  Get the total number of quadrature points in each element (ie num_quad_in_each_dir^DIM)
     */
    unsigned GetNumQuadPointsPerElement()
    {
        if(DIM==1)
        {
            return mNumQuadPointsInEachDimension;
        }
        if(DIM==2)
        {
            return mNumQuadPointsInEachDimension*mNumQuadPointsInEachDimension;
        }
        else //DIM==3
        {
            return mNumQuadPointsInEachDimension*mNumQuadPointsInEachDimension*mNumQuadPointsInEachDimension;
        }
    }


    /**
     *  Get lambda (the stretch ratio).
     *
     *  NOTE: the i-th entry of this vector is
     *  assumed to be the i-th quad point obtained by looping over cells in the obvious
     *  way and then looping over quad points
     */
    std::vector<double>& rGetLambda()
    {
        return mLambda;
    }

    /**
     *  Solve between given times. A static explicit assembler will ignore
     *  the times
     */
    virtual void Solve(double currentTime, double nextTime, double timestep)=0;

    /**
     *  Set the appropriate forcing quantity, by gauss point (active tension
     *  if explicit, [Ca]_i if implicit
     */
    virtual void SetForcingQuantity(std::vector<double>& forcingQuantity)=0;
};


#endif /*ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_*/
