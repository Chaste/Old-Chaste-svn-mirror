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


#ifndef GENERALPLANESTIMULUSCELLFACTORY_HPP_
#define GENERALPLANESTIMULUSCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"

template <class CELL, unsigned DIM>
class GeneralPlaneStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    /** The stimulus that gets applied. */
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    GeneralPlaneStimulusCellFactory(unsigned numEleAcross, double meshWidth, bool useMeshWidthAsMag=false)
        : AbstractCardiacCellFactory<DIM>()
    {
        ///\todo The useMeshWidth is temporary, while we are sorting out
        ///3D stimulus.  It is to be removed later (along with StimulusConvergenceTester)
        /// scale stimulus depending on space_step of elements
        ///\todo It looks like the value of the stimulus is specific to 3D
        if (useMeshWidthAsMag)
        {
            #define COVERAGE_IGNORE
            mpStimulus.reset(new SimpleStimulus(meshWidth, 0.5));
            #undef COVERAGE_IGNORE
        }
        else
        {
            double stimulus_magnitude=-1e7;//wrt mesh 4 which has 64 elements in 1D
            switch(DIM)
            {
                case 1:
                {
                    stimulus_magnitude*=numEleAcross/(64.0);
                    //Justification: elements go half size with each refinement
                    stimulus_magnitude*=meshWidth/(0.2);
                    break;
                }
                case 2:
                {
                    stimulus_magnitude*=numEleAcross/(64.0);
                    //Justification: Triangles go quarter size with each refinement, but there are twice as many nodes on boundary
                    stimulus_magnitude*=meshWidth/(0.2);
                    break;
                }
                default: //3D
                {
                    stimulus_magnitude*=numEleAcross/(64.0);
                    //Hypothesis: Triangles go eighth size with each refinement, but there are four-times as many nodes on boundary
                    stimulus_magnitude*=meshWidth/(0.2);
                    break;
                }
            }
            //std::cout<<"Mag is "<<stimulus_magnitude<<"\n";
            mpStimulus.reset(new SimpleStimulus(stimulus_magnitude, 0.5));
        }
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        double x = this->GetMesh()->GetNode(node)->GetPoint()[0];
        if (x*x<=1e-10)
        {
            return new CELL(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CELL(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

#endif /*GENERALPLANESTIMULUSCELLFACTORY_HPP_*/
