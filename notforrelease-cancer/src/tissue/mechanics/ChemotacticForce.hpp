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
#ifndef CHEMOTACTICFORCE_HPP_
#define CHEMOTACTICFORCE_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

/// \todo This class needs documenting (see #736)
template<unsigned DIM>
class ChemotacticForce  : public AbstractForce<DIM>
{
friend class TestForcesNotForRelease;

private:

    /// \todo This method needs documenting (see #736)
    double GetChemotacticForceMagnitude(const double concentration, const double concentrationGradientMagnitude);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

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

    /** Use springs which are dependent on whether cells are necrotic */
    bool mUseApoptoticSprings;

public:

    /**
     * Constructor.
     */
    ChemotacticForce();

    /**
     * Destructor.
     */
    ~ChemotacticForce();

    /**
     * Calculates the forces on each node
     *
     * @return the force components on each node. Of size NUM_NODES x DIM.
     *
     * Fc = chi(C,|gradC|) gradC/|gradC|  (if |gradC|>0, else Fc = 0)
     *
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                 AbstractTissue<DIM>& rTissue);

};

#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemotacticForce)

#endif /*CHEMOTACTICFORCE_HPP_*/
