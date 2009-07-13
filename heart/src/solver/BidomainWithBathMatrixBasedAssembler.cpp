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

#include "BidomainWithBathMatrixBasedAssembler.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>

#include "ConstBoundaryCondition.hpp"
#include "HeartRegionCodes.hpp"


/////////////////////////////////////////////////////////////////////
// BidomainWithBathRhsMatrixAssembler
/////////////////////////////////////////////////////////////////////

template<unsigned DIM>
c_matrix<double,2*(DIM+1),2*(DIM+1)> BidomainWithBathRhsMatrixAssembler<DIM>::ComputeMatrixTerm(
            c_vector<double, DIM+1> &rPhi,
            c_matrix<double, DIM, DIM+1> &rGradPhi,
            ChastePoint<DIM> &rX,
            c_vector<double,2> &u,
            c_matrix<double,2,DIM> &rGradU /* not used */,
            Element<DIM,DIM>* pElement)
{
    c_matrix<double,2*(DIM+1),2*(DIM+1)> ret = zero_matrix<double>(2*(DIM+1), 2*(DIM+1));

    if (pElement->GetRegion() != HeartRegionCode::BATH)
    {
        c_matrix<double, DIM+1, DIM+1> basis_outer_prod = outer_prod(rPhi, rPhi);

        // even rows, even columns
        matrix_slice<c_matrix<double, 2*DIM+2, 2*DIM+2> >
        slice00(ret, slice (0, 2, DIM+1), slice (0, 2, DIM+1));
        slice00 =  basis_outer_prod;

        // odd rows, even columns: are zero
        // even rows, odd columns: are zero
        // odd rows, odd columns: are zero
    }

    return ret;
}

template<unsigned DIM>
c_vector<double,2*(DIM+1)> BidomainWithBathRhsMatrixAssembler<DIM>::ComputeVectorTerm(
    c_vector<double, DIM+1> &rPhi,
    c_matrix<double, DIM, DIM+1> &rGradPhi,
    ChastePoint<DIM> &rX,
    c_vector<double,2> &u,
    c_matrix<double, 2, DIM> &rGradU /* not used */,
    Element<DIM,DIM>* pElement)

{
    #define COVERAGE_IGNORE
    NEVER_REACHED;
    return zero_vector<double>(2*(DIM+1));
    #undef COVERAGE_IGNORE
}


template<unsigned DIM>
c_vector<double, 2*DIM> BidomainWithBathRhsMatrixAssembler<DIM>::ComputeVectorSurfaceTerm(
    const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
    c_vector<double, DIM> &rPhi,
    ChastePoint<DIM> &rX )
{
    #define COVERAGE_IGNORE
    NEVER_REACHED;
    return zero_vector<double>(2*DIM);
    #undef COVERAGE_IGNORE
}


template<unsigned DIM>
BidomainWithBathRhsMatrixAssembler<DIM>::BidomainWithBathRhsMatrixAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh)
    : AbstractLinearAssembler<DIM,DIM,2,false,BidomainWithBathRhsMatrixAssembler<DIM> >()
{
    this->mpMesh = pMesh;

    // this needs to be set up, though no boundary condition values are used in the matrix
    this->mpBoundaryConditions = new BoundaryConditionsContainer<DIM,DIM,2>;
    this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(pMesh);

    //DistributedVector::SetProblemSize(this->mpMesh->GetNumNodes()); WOULD BE WRONG -- we need the maintain an uneven distribution, if given
    Vec template_vec = this->mpMesh->GetDistributedVectorFactory()->CreateVec(2);
    this->mpLinearSystem = new LinearSystem(template_vec);
    VecDestroy(template_vec);


    this->AssembleSystem(false,true);
}

template<unsigned DIM>
BidomainWithBathRhsMatrixAssembler<DIM>::~BidomainWithBathRhsMatrixAssembler()
{
    delete this->mpBoundaryConditions;
}

template<unsigned DIM>
Mat* BidomainWithBathRhsMatrixAssembler<DIM>::GetMatrix()
{
    return &(this->mpLinearSystem->rGetLhsMatrix());
}

/////////////////////////////////////////////////////////////////////
// BidomainMatrixBasedAssembler
/////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BidomainWithBathMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>::BidomainWithBathMatrixBasedAssembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            BidomainPde<SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
            unsigned numQuadPoints)
    : BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints),
      BidomainMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints),
      BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints)


{
    // construct matrix using the helper class
    mpBidomainWithBathRhsMatrixAssembler = new BidomainWithBathRhsMatrixAssembler<SPACE_DIM>(pMesh);
    this->mpMatrixForMatrixBasedRhsAssembly = mpBidomainWithBathRhsMatrixAssembler->GetMatrix();

    /// \todo: at this point we'll have a BidomainWithBathRhsMatrixAssembler object and a BidomainRhsMatrixAssembler
    /// object in memory. This is a waste of memory since both construct and store a matrix for RHS assembly.
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BidomainWithBathMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>::~BidomainWithBathMatrixBasedAssembler()
{
    delete mpBidomainWithBathRhsMatrixAssembler;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainWithBathMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>::ConstructVectorForMatrixBasedRhsAssembly(
        Vec currentSolution)
{
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
    DistributedVector::Stripe distributed_current_solution_vm(distributed_current_solution, 0);

    // dist stripe for z
    DistributedVector dist_vec_matrix_based = p_factory->CreateDistributedVector(this->mVectorForMatrixBasedRhsAssembly);
    DistributedVector::Stripe dist_vec_matrix_based_vm(dist_vec_matrix_based, 0);
    DistributedVector::Stripe dist_vec_matrix_based_phie(dist_vec_matrix_based, 1);

    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm  = HeartConfig::Instance()->GetCapacitance();

    for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
         index != dist_vec_matrix_based.End();
         ++index)
    {
        if (this->mpMesh->GetNode(index.Global)->GetRegion() != HeartRegionCode::BATH)
        {
            double V = distributed_current_solution_vm[index];
            double F = - Am*this->mpBidomainPde->rGetIionicCacheReplicated()[index.Global]
                       - this->mpBidomainPde->rGetIntracellularStimulusCacheReplicated()[index.Global];

            dist_vec_matrix_based_vm[index] = Am*Cm*V*this->mDtInverse + F;
        }
        else
        {
            dist_vec_matrix_based_vm[index] = 0.0;
        }

        dist_vec_matrix_based_phie[index] = 0.0;
    }

    dist_vec_matrix_based.Restore();

    VecAssemblyBegin(this->mVectorForMatrixBasedRhsAssembly);
    VecAssemblyEnd(this->mVectorForMatrixBasedRhsAssembly);
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class BidomainWithBathMatrixBasedAssembler<1,1>;
template class BidomainWithBathMatrixBasedAssembler<2,2>;
template class BidomainWithBathMatrixBasedAssembler<3,3>;

template class BidomainWithBathRhsMatrixAssembler<1>;
template class BidomainWithBathRhsMatrixAssembler<2>;
template class BidomainWithBathRhsMatrixAssembler<3>;
