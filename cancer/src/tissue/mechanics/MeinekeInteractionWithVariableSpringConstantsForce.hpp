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
#ifndef MEINEKEINTERACTIONWITHVARIABLESPRINGCONSTANTSFORCE_HPP_
#define MEINEKEINTERACTIONWITHVARIABLESPRINGCONSTANTSFORCE_HPP_

#include "MeinekeInteractionForce.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 * A subclass of MeinekeInteractionForce with variable spring constants.
 */
template<unsigned DIM>
class MeinekeInteractionWithVariableSpringConstantsForce : public MeinekeInteractionForce<DIM>
{
    friend class TestForces;
    
private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<MeinekeInteractionForce<DIM> >(*this);
        archive & mUseEdgeBasedSpringConstant;
        archive & mUseMutantSprings;
        archive & mMutantMutantMultiplier;
        archive & mNormalMutantMultiplier;
        archive & mUseBCatSprings;
        archive & mUseApoptoticSprings;
    }

protected :
    
    /** Whether to use spring constant proportional to cell-cell contact length/area (defaults to false) */
    bool mUseEdgeBasedSpringConstant;

    /** Whether to use different stiffnesses depending on whether either cell is a mutant */
    bool mUseMutantSprings;

    /** Multiplier for spring stiffness if mutant */
    double mMutantMutantMultiplier;

    /** Multiplier for spring stiffness if mutant */
    double mNormalMutantMultiplier;

    /** Use springs which are dependent on beta-catenin levels */
    bool mUseBCatSprings;

    /** Use springs which are dependent on whether cells are apoptotic */
    bool mUseApoptoticSprings;

public :

    /**
     * Constructor.
     */
    MeinekeInteractionWithVariableSpringConstantsForce();
    
    /**
     * Destructor.
     */
    ~MeinekeInteractionWithVariableSpringConstantsForce();

    /**
     * Use an edge-based spring constant
     */
    void SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant);

    /**
     * Use Different spring strengths depending on two cells:
     * Normal-normal, Normal-mutant, mutant-mutant
     */
    void SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier=2, double normalMutantMultiplier=1.5);

    /**
     * Use the amount of B-Catenin on an edge to find spring constant.
     */
    void SetBCatSprings(bool useBCatSprings);

    /**
     * Set spring stiffness to be dependent on whether cells are necrotic
     */
    void SetApoptoticSprings(bool useApoptoticSprings);    

    /**
     * Return a multiplication factor for the spring constant, which 
     * may depend on whether the given pair of neighbouring cells are 
     * e.g. undergoing apoptosis, have mutations, or experience variable 
     * levels of beta catenin.
     */     
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex, 
                                                      unsigned nodeBGlobalIndex, 
                                                      AbstractTissue<DIM>& rTissue, 
                                                      bool isCloserThanRestLength);
    
    /**
     * Overridden AddForceContribution method.
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractTissue<DIM>& rTissue);
 
};

#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeinekeInteractionWithVariableSpringConstantsForce)


#endif /*MEINEKEINTERACTIONWITHVARIABLESPRINGCONSTANTSFORCE_HPP_*/
