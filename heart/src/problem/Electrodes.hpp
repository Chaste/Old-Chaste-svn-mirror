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

#ifndef ELECTRODES_HPP_
#define ELECTRODES_HPP_

#include "AbstractTetrahedralMesh.hpp"
#include "DistributedVector.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"


/**
 *  A class for setting up the boundary conditions associated with electrodes.
 *  There are two modes: grounding the second electrode, in which case the first
 *  electrode has an input flux (Neumann boundary condition) of the specified
 *  magnitude, the extracellular potential is fixed on the second electrode; or
 *  not grounded, in which case the opposite electrode has an equal and opposite
 *  output flux.
 *
 *  This class assumes the given mesh is cuboid, and the electrodes are taken to be
 *  the specified opposite surfaces.
 */
template<unsigned DIM>
class Electrodes
{
friend class TestBidomainWithBathAssembler;

private:
    /** Whether the second electrode is grounded */
    bool mGroundSecondElectrode;
    /** The created bcc, which BidomainProblem will use */
    BoundaryConditionsContainer<DIM,DIM,2>* mpBoundaryConditionsContainer;
    /** The time the electrodes are switched off */
    double mEndTime;
    /** Whether the electrodes are currently switched on */
    bool mAreActive;

public:
    /** Constructor.
     *  @param rMesh The mesh, assumed to be a cuboid.
     *  @param groundSecondElectrode Whether to ground the second electrode (see class documentation)
     *  @param index The value i when applying the electrodes to x_i=a and x_i=b (a<b)
     *  @param lowerValue The value a when applying the electrodes to x_i=a and x_i=b (a<b) (should
     *    be the minimum value of x_i for the given mesh)
     *  @param upperValue The value b when applying the electrodes to x_i=a and x_i=b (a<b) (should
     *    be the maximum value of x_i for the given mesh)
     *  @param magnitude Magnitude of the stimulus
     *  @param duration Duration of the stimulus. Note, start time currently assumed to be zero.
     */
    Electrodes(AbstractTetrahedralMesh<DIM,DIM>& rMesh,
               bool groundSecondElectrode,
               unsigned index, double lowerValue, double upperValue,
               double magnitude, double duration); // implemented in cpp

    /**
     *  Delete the set up bcc
     */
    ~Electrodes()
    {
        delete mpBoundaryConditionsContainer;
    }

    /**
     *  Get the boundary conditions container in which is set up the Neumann
     *  fluxes for the first electrode, and the opposite fluxes for the second
     *  electrode if the the second electrode isn't grounded
     */
    BoundaryConditionsContainer<DIM,DIM,2>* GetBoundaryConditionsContainer()
    {
        assert(mAreActive);
        return mpBoundaryConditionsContainer;
    }

    /**
     *  Whether it is time to switch off the electrodes yet. THIS ONLY RETURNS
     *  TRUE ONCE - the first appropriate time. After that the electrodes assume
     *  they have been switched off and therefore this returns false.
     * 
     * @param time  the current time
     */
    bool SwitchOff(double time)
    {
        if(mAreActive && time>mEndTime)
        {
            mAreActive = false;
            return true;
        }

        return false;
    }
};

#endif /*ELECTRODES_HPP_*/
