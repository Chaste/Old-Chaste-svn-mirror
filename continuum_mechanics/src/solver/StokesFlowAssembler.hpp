/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef STOKESFLOWASSEMBLER_HPP_
#define STOKESFLOWASSEMBLER_HPP_

#include "AbstractContinuumMechanicsAssembler.hpp"
#include "StokesFlowProblemDefinition.hpp"

///\todo #1956 dox, cpp, test with non-zero body force

/**
 *  TODO DOX
 */
template<unsigned DIM>
class StokesFlowAssembler : public AbstractContinuumMechanicsAssembler<DIM,true,true>
{
private:
    /** TODO DOX  */
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;

    /** Number of nodes per element. */
    static const unsigned NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic

    /** TODO DOX  */
    static const unsigned SPATIAL_BLOCK_SIZE_ELEMENTAL = DIM*NUM_NODES_PER_ELEMENT;
    /** TODO DOX  */
    static const unsigned PRESSURE_BLOCK_SIZE_ELEMENTAL = NUM_VERTICES_PER_ELEMENT;

    /** TODO DOX  */
    StokesFlowProblemDefinition<DIM>* mpProblemDefinition;

    /** TODO DOX  */
    c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialSpatialMatrixTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        ChastePoint<DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> ret = zero_matrix<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL);

        for (unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
        {
            unsigned spatial_dim1 = index1%DIM;
            unsigned node_index1 = (index1-spatial_dim1)/DIM;

            for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
            {
                unsigned spatial_dim2 = index2%DIM;
                unsigned node_index2 = (index2-spatial_dim2)/DIM;

                if (spatial_dim1 == spatial_dim2)
                {
                    double grad_quad_phi_grad_quad_phi = 0.0;
                    for (unsigned k=0; k<DIM; k++)
                    {
                        grad_quad_phi_grad_quad_phi += rGradQuadPhi(k, node_index1) * rGradQuadPhi(k, node_index2);
                    }

                    ret(index1,index2) += mpProblemDefinition->GetViscosity() * grad_quad_phi_grad_quad_phi;
                }
            }
        }
        return ret;

    }

    /** TODO DOX  */
    c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputeSpatialPressureMatrixTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
        ChastePoint<DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ret = zero_matrix<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL);

        for (unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
        {
            unsigned spatial_dim1 = index1%DIM;
            unsigned node_index1 = (index1-spatial_dim1)/DIM;

            for (unsigned index2=0; index2<NUM_VERTICES_PER_ELEMENT; index2++)
            {
                ret(index1,index2) += -rGradQuadPhi(spatial_dim1, node_index1) * rLinearPhi(index2);
            }
        }

        return ret;
    }

    // We don't implement this method - so it is a zero block
    //c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressurePressureMatrixTerm(
    //    c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
    //    c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
    //    ChastePoint<DIM>& rX,
    //    Element<DIM,DIM>* pElement)


    /** TODO DOX  */
    c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialVectorTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        ChastePoint<DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ret = zero_vector<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL);

        c_vector<double,DIM> body_force;
        switch (mpProblemDefinition->GetBodyForceType())
        {
            case CONSTANT_BODY_FORCE:
            {
                body_force = mpProblemDefinition->GetConstantBodyForce();
                break;
            }
            default:
                NEVER_REACHED;
        }

        for (unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
        {
            unsigned spatial_dim = index%DIM;
            unsigned node_index = (index-spatial_dim)/DIM;

            ret(index) += body_force(spatial_dim) * rQuadPhi(node_index);
        }

        return ret;
    }

    // We don't implement this method - so it is a zero block of the vector:
    //c_vector<double,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressureVectorTerm(
    //        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
    //        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
    //        ChastePoint<DIM>& rX,
    //        Element<DIM,DIM>* pElement)

public:
    /** TODO DOX  */
    StokesFlowAssembler(QuadraticMesh<DIM>* pMesh,
                        StokesFlowProblemDefinition<DIM>* pProblemDefinition)
        : AbstractContinuumMechanicsAssembler<DIM,true,true>(pMesh),
          mpProblemDefinition(pProblemDefinition)
    {
    }
};

#endif // STOKESFLOWASSEMBLER_HPP_
