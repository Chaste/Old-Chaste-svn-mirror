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

#ifndef BIDOMAINFASTSLOWPDE_HPP_
#define BIDOMAINFASTSLOWPDE_HPP_

#include "AbstractCardiacFastSlowPde.hpp"
#include "BidomainPde.hpp"


/**
 *  Class which specifies and solves a bidomain problem using the FAST/SLOW
 *  CELL ALGORITHM [J. Whiteley] on a coarse/fine mesh.
 */ 
template <unsigned DIM>
class BidomainFastSlowPde : public AbstractCardiacFastSlowPde<DIM>, public BidomainPde<DIM>
{
friend class TestCardiacFastSlowPde;


public:
    /**
     *  Note: we have BidomainFastSlowPde inheriting from both AbstractCardiacFastSlowPde and
     * BidomainPde, both of which inherit from AbstractCardiacPde (a 'dreaded diamond'). Therefore,
     * both of those have to use virtual inheritence of AbstractCardiacPde, and also the constructor
     * of AbstractCardiacPde (ie the top base class) has to be explicitly called before the other
     * constructors.
     * @param pCellFactory Cell factory associated with a fine mesh
     * @param rMixedMesh Coarse/fine mesh
     * @param slowCurrentsTimeStep Time step on which to interpolate slow currents from the coarse mesh onto the fine mesh (multiple of PDE timestep)
     */
     BidomainFastSlowPde(AbstractCardiacCellFactory<DIM>* pCellFactory,
                        MixedTetrahedralMesh<DIM,DIM>& rMixedMesh,
                        double slowCurrentsTimeStep)
            :  AbstractCardiacPde<DIM>(pCellFactory),
               AbstractCardiacFastSlowPde<DIM>(pCellFactory, rMixedMesh, slowCurrentsTimeStep),
               BidomainPde<DIM>(pCellFactory)
    {
    }
};


#endif /*BIDOMAINFASTSLOWPDE_HPP_*/
