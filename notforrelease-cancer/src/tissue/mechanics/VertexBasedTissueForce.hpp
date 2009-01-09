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
#ifndef VERTEXBASEDTISSUEFORCE_HPP_
#define VERTEXBASEDTISSUEFORCE_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

/**
 * A force class for use in vertex-based tissue simulations
 */
template<unsigned DIM>
class VertexBasedTissueForce  : public AbstractForce<DIM>
{
friend class TestForcesNotForRelease;

private:

    /// \todo Should probably put default values for these variables into CancerParameters,
    ///       and generalize mCellCellAdhesionEnergyParameter to be cell-type dependent
    ///       (see #861)    
     
    /** The target area of a 'fully-grown' cell in the tissue. */
    double mTissueCellTargetArea;
    
    /** The deformation energy parameter (denoted by lambda in my notes - see #861). */
    double mDeformationEnergyParameter;
    
    /** The membrane surface energy parameter (denoted by beta in my notes - see #861). */
    double mMembraneSurfaceEnergyParameter;
    
    /** The cell-cell adhesion energy parameter (denoted by gamma in my notes - see #861). */
    double mCellCellAdhesionEnergyParameter;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mTissueCellTargetArea;
        archive & mDeformationEnergyParameter;
        archive & mMembraneSurfaceEnergyParameter;
        archive & mCellCellAdhesionEnergyParameter;               
    }

public:

    /**
     * Constructor.
     * 
     * @param tissueCellTargetArea the cell target area (defaults to 1)
     *                             \todo should probably make this a CancerParameter (see #861) 
     */
    VertexBasedTissueForce(double tissueCellTargetArea=1.0);

    /**
     * Destructor.
     */
    ~VertexBasedTissueForce();

    /**
     * Overridden AddForceContribution() method.
     * 
     * Calculates the force on each node in the vertex-based tissue.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rTissue reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractTissue<DIM>& rTissue);

    /**
     * Set method for mDeformationEnergyParameter.
     * 
     * @param deformationEnergyParameter the value for mDeformationEnergyParameter
     */
    void SetDeformationEnergyParameter(double deformationEnergyParameter);
    
    /**
     * Set method for mMembraneSurfaceEnergyParameter.
     * 
     * @param membraneSurfaceEnergyParameter the value for mMembraneSurfaceEnergyParameter
     */
    void SetMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter);

    /**
     * Set method for mCellCellAdhesionEnergyParameter.
     * 
     * @param cellCellAdhesionEnergyParameter the value for mCellCellAdhesionEnergyParameter
     */
    void SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter);
};

                              
#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedTissueForce)

#endif /*VERTEXBASEDTISSUEFORCE_HPP_*/
