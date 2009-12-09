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
#ifndef SIMPLEBATHPROBLEMSETUP_HPP_
#define SIMPLEBATHPROBLEMSETUP_HPP_

/**
 * @file
 * 
 * Some helper classes and functions for setting up a simple bath problem for testing.
 */

#include "AbstractCardiacCellFactory.hpp"
#include "SimpleStimulus.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "HeartRegionCodes.hpp"

/**
 * A simple cell factory for bath problems, applying a SimpleStimulus for
 * 0.5ms at a single point.
 */
template<unsigned DIM>
class BathCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    /** The stimulus to apply */
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    /** and where to apply it */
    c_vector<double,DIM> mStimulatedPoint;

public:
    /**
     * Create the cell factory.
     * @param stimulusMagnitude
     * @param stimulatedPoint spatial co-ordinates of where to stimulate.
     *    Must correspond to a node location.
     */
    BathCellFactory(double stimulusMagnitude, c_vector<double,DIM> stimulatedPoint)
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(stimulusMagnitude, 0.5)),
          mStimulatedPoint(stimulatedPoint)
    {
    }

    /**
     * Create a new cell.
     * @param node  where to put the cell
     */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        // paranoia - check this is really a tissue node
        assert(this->GetMesh()->GetNode(node)->GetRegion() == HeartRegionCode::TISSUE);

        // stimulate centre node normally..
        bool is_centre;

        if (DIM==1)
        {
            is_centre = (fabs(this->GetMesh()->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6);
        }
        else if (DIM==2)
        {
            is_centre = (    (fabs(this->GetMesh()->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6)
                          && (fabs(this->GetMesh()->GetNode(node)->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6) );
        }
        else
        {
            is_centre = (    (fabs(this->GetMesh()->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6)
                          && (fabs(this->GetMesh()->GetNode(node)->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6)
                          && (fabs(this->GetMesh()->GetNode(node)->GetPoint()[2]-mStimulatedPoint(2)) < 1e-6) );
        }

        if (is_centre)
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

/**
 * Load a 2d mesh, and set everything outside a central circle to be bath.
 *
 * @param rMeshPath relative path to the mesh
 * @param centreX X co-ord of tissue centre
 * @param centreY Y co-ord of tissue centre
 * @param radius radius of tissue
 */
template<class MeshType>
MeshType* Load2dMeshAndSetCircularTissue(const std::string& rMeshPath,
                                         double centreX, double centreY, double radius)
{
    TrianglesMeshReader<2,2> reader(rMeshPath);
    MeshType* p_mesh = new MeshType;
    p_mesh->ConstructFromMeshReader(reader);

    for (typename MeshType::ElementIterator it = p_mesh->GetElementIteratorBegin();
         it != p_mesh->GetElementIteratorEnd();
         ++it)
    {
        double x = it->CalculateCentroid()[0];
        double y = it->CalculateCentroid()[1];
        if ( (x-centreX)*(x-centreX) + (y-centreY)*(y-centreY) > radius*radius )
        {
            it->SetRegion(HeartRegionCode::BATH);
        }
    }
    
    return p_mesh;
}

#endif /*SIMPLEBATHPROBLEMSETUP_HPP_*/
