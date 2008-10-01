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

#ifndef NONLINEARELASTICITYTOOLS_HPP_
#define NONLINEARELASTICITYTOOLS_HPP_

#include "ConformingTetrahedralMesh.hpp"
    
/**
 *  A class of helper methods for problems which use NonlinearElasticityAssembler
 */
template<unsigned DIM>
class NonlinearElasticityTools
{
public:
    /** 
     *  Collect all the nodes which satisfy x[k] = c, for given k and c, in order
     *  to be set as fixed (or displacement boundary condition) nodes. Note that
     *  this method does not check if the nodes on the required surface are actually
     *  boundary nodes. It does however throw an exception if no nodes on the given
     *  surface are found.
     */
    static std::vector<unsigned> GetNodesByComponentValue(ConformingTetrahedralMesh<DIM,DIM>& rMesh,
                                                          unsigned component,
                                                          double value)
    {
        std::vector<unsigned> fixed_nodes;
        double tol = 1e-8;
        for(unsigned i=0; i<rMesh.GetNumNodes(); i++)
        {
            if( fabs(rMesh.GetNode(i)->rGetLocation()[component] - value)<1e-8)
            {
                fixed_nodes.push_back(i);
            }
        }
        
        if(fixed_nodes.size()==0)
        {
            std::stringstream error;
            error << "Could not find any nodes on requested surface (note: tolerance = "<<tol<<")";
            EXCEPTION(error.str());
        }
        
        return fixed_nodes;
    }
};
        


#endif /*NONLINEARELASTICITYTOOLS_HPP_*/
