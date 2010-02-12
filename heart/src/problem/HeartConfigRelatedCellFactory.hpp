/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef HEARTCONFIGRELATEDCELLFACTORY_HPP_
#define HEARTCONFIGRELATEDCELLFACTORY_HPP_

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCellFactory.hpp"
#include "HeartConfig.hpp"

#include "AbstractStimulusFunction.hpp"
#include "MultiStimulus.hpp" // Included here for archiving - see below.
#include "SimpleStimulus.hpp"

#include "ChasteCuboid.hpp"

/*
 * Even though these classes are only used in the .cpp file, they need to be
 * included here for serialization to work - the archiving code needs to see
 * the CHASTE_CLASS_EXPORT incantations.
 */
#include "BackwardEulerFoxModel2002Modified.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "FoxModel2002Modified.hpp"
#include "FaberRudy2000Version3.hpp"
#include "FaberRudy2000Version3Optimised.hpp"
#include "DiFrancescoNoble1985OdeSystem.hpp"
#include "Mahajan2008OdeSystem.hpp"
#include "TenTusscher2006OdeSystem.hpp"
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "Maleckar2009OdeSystem.hpp"
#include "AbstractChasteRegion.hpp"

/**
 *  \todo: dox, coverage, maybe own tests
 */
template<unsigned SPACE_DIM>
class HeartConfigRelatedCellFactory : public AbstractCardiacCellFactory<SPACE_DIM>
{
private:
    /** Default cardiac cell model to be used in all tissue (except heterogeneous regions)*/
    cp::ionic_model_selection_type mDefaultIonicModel;
    /** List of axis-aligned box regions which contain heterogeneous cardiac ionic model types*/
    std::vector<ChasteCuboid<SPACE_DIM> > mIonicModelRegions;
    /** List of ionic model (size matches that of mIonicModelRegions)*/
    std::vector<cp::ionic_model_selection_type> mIonicModelsDefined;
    
    /** List of axis-aligned box regions which represent areas to stimulate*/
    std::vector<ChasteCuboid<SPACE_DIM> > mStimulatedAreas;
    /** List of intracellular current stimuli to apply (size matches that of mStimulatedAreas)*/
    std::vector<boost::shared_ptr<SimpleStimulus> > mStimuliApplied;
    
    /** 
     *  List of regions which represent areas in which to give parametric heterogeneity (scaling gating parameters)
     *  This vector will be filled in by the HeartConfig::GetCellHeterogeneity method if the user requested
     *  to specify the heterogeneity areas by cuboids, or, alternatively, by the FillInCellularTransmuralAreas method
     *  in the cell factory called by the problem class AFTER setting the mesh 
     *  (which is needed for the calculations of the distance maps for the calculations of heterogeneities).
     *  
     *  When creating a cardiac cell for each node (CreateCardiacCellForTissueNode) the code will check whether 
     *  that node is contained in the heterogeneity area or not. 
     * 
     */
    std::vector<AbstractChasteRegion<SPACE_DIM>* > mCellHeterogeneityAreas;
    /** List of scale factors for Gks scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorGks;
    /** List of scale factors for Ito scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorIto;
    /** List of scale factors for Gkr scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorGkr;
    
public:
    /** Default constructor */
    HeartConfigRelatedCellFactory();
    
    /** Destructor*/
    ~HeartConfigRelatedCellFactory();

    /**
     * Create the correct tissue cell for a given region in the mesh
     * @param intracellularStimulus is computed in CreateCardiacCellForTissueNode determined by the list of stimulation regions
     * @param nodeIndex is the global index within the mesh
     */
    AbstractCardiacCell* CreateCellWithIntracellularStimulus(
            boost::shared_ptr<AbstractStimulusFunction> intracellularStimulus,
            unsigned nodeIndex);

    /**
     * Create the correct stimulated tissue cell for a given region in the mesh
     * The stimulus is determined in this method (using the list of stimulation regions).  
     * The cardiac cell type (and parameters) are 
     * determined in the CreateCellWithIntracellularStimulus method
     * @param nodeIndex is the global index within the mesh
     */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex);
    
    /**
     * Helper method to calculate and fill in the heterogeneities areas (mCellHeterogeneityAreas)
     */
    void FillInCellularTransmuralAreas();

};


#endif /*HEARTCONFIGRELATEDCELLFACTORY_HPP_*/
