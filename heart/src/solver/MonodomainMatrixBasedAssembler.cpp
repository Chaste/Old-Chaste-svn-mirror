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

#include "MonodomainMatrixBasedAssembler.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//   NOTE: The main class in this file, MonodomainMatrixBasedAssembler, is defined at the
//   bottom, after MonodomainRhsMatrixAssembler
//
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////
// MonodomainRhsMatrixAssembler
/////////////////////////////////////////////////////////////////////
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> MonodomainRhsMatrixAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
    c_vector<double, ELEMENT_DIM+1> &rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
    ChastePoint<SPACE_DIM> &rX,
    c_vector<double,1> &u,
    c_matrix<double,1,SPACE_DIM> &rGradU /* not used */,
    Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    return outer_prod(rPhi, rPhi);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,1*(ELEMENT_DIM+1)> MonodomainRhsMatrixAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(
    c_vector<double, ELEMENT_DIM+1> &rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
    ChastePoint<SPACE_DIM> &rX,
    c_vector<double,1> &u,
    c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
    Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    #define COVERAGE_IGNORE
    NEVER_REACHED;
    return zero_vector<double>(SPACE_DIM+1);
    #undef COVERAGE_IGNORE
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, ELEMENT_DIM> MonodomainRhsMatrixAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeVectorSurfaceTerm(
    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
    c_vector<double, ELEMENT_DIM> &rPhi,
    ChastePoint<SPACE_DIM> &rX )
{
    #define COVERAGE_IGNORE
    NEVER_REACHED;
    return zero_vector<double>(SPACE_DIM);
    #undef COVERAGE_IGNORE
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainRhsMatrixAssembler<ELEMENT_DIM,SPACE_DIM>::MonodomainRhsMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,1,false,MonodomainRhsMatrixAssembler<ELEMENT_DIM,SPACE_DIM> >()
{
    this->mpMesh = pMesh;

    // this needs to be set up, though no boundary condition values are used in the matrix
    this->mpBoundaryConditions = new BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>;
    this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(pMesh);

    //This linear system needs a distribution from the DistributedVector class
    Vec temp_vec = pMesh->GetDistributedVectorFactory()->CreateVec();
    unsigned preallocation = this->mpMesh->CalculateMaximumNodeConnectivityPerProcess();
    this->mpLinearSystem = new LinearSystem(temp_vec, preallocation);
    VecDestroy(temp_vec);
    this->AssembleSystem(false,true);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainRhsMatrixAssembler<ELEMENT_DIM, SPACE_DIM>::~MonodomainRhsMatrixAssembler()
{
    delete this->mpBoundaryConditions;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Mat* MonodomainRhsMatrixAssembler<ELEMENT_DIM,SPACE_DIM>::GetMatrix()
{
    return &(this->mpLinearSystem->rGetLhsMatrix());
}





/////////////////////////////////////////////////////////////////////
// MonodomainMatrixBasedAssembler
/////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>::MonodomainMatrixBasedAssembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            MonodomainPde<ELEMENT_DIM,SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBcc,
            unsigned numQuadPoints)
    : MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints)
{
    // construct matrix using the helper class
    mpMonodomainRhsMatrixAssembler = new MonodomainRhsMatrixAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh);
    this->mpMatrixForMatrixBasedRhsAssembly = mpMonodomainRhsMatrixAssembler->GetMatrix();

    // set variables on parent class so that we do matrix-based assembly, and allocate
    // memory for the vector 'z'
    this->mUseMatrixBasedRhsAssembly = true;
    this->mVectorForMatrixBasedRhsAssembly = pMesh->GetDistributedVectorFactory()->CreateVec();

    // Tell pde there's no need to replicate ionic caches
    pPde->SetCacheReplication(false);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>::~MonodomainMatrixBasedAssembler()
{
    delete mpMonodomainRhsMatrixAssembler;
    VecDestroy(this->mVectorForMatrixBasedRhsAssembly);
}

//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//void MonodomainMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>::ConstructVectorForMatrixBasedRhsAssembly(Vec existingSolution)
//{
//    // copy V to z
//    VecCopy(existingSolution, this->mVectorForMatrixBasedRhsAssembly);
//
//    // set up a vector which has the nodewise force term (ie A*I_ionic+I_stim)
//    Vec force_term_at_nodes = this->mpMesh->GetDistributedVectorFactory()->CreateVec();
//    PetscInt lo, hi;
//    VecGetOwnershipRange(force_term_at_nodes, &lo, &hi);
//    double* p_force_term;
//    VecGetArray(force_term_at_nodes, &p_force_term);
//    for (int global_index=lo; global_index<hi; global_index++)
//    {
//        unsigned local_index = global_index - lo;
//        const Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(global_index);
//        p_force_term[local_index] = this->mpMonodomainPde->ComputeNonlinearSourceTermAtNode(*p_node, 0.0);
//    }
//    VecRestoreArray(force_term_at_nodes, &p_force_term);
//    VecAssemblyBegin(force_term_at_nodes);
//    VecAssemblyEnd(force_term_at_nodes);
//
//    double one=1.0;
//    double scaling=  this->mpMonodomainPde->ComputeDuDtCoefficientFunction(ChastePoint<SPACE_DIM>())
//                    *this->mDtInverse;
//
//#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
//    // VecAXPBY(a,b,x,y) does y = ax + by
//    VecAXPBY(&one,
//             &scaling,
//             force_term_at_nodes,
//             this->mVectorForMatrixBasedRhsAssembly);
//#else
//
//    // VecAXPBY(y,a,b,x) does y = ax + by
//    VecAXPBY(this->mVectorForMatrixBasedRhsAssembly,
//             one,
//             scaling,
//             force_term_at_nodes);
//#endif
//
//    VecAssemblyBegin(this->mVectorForMatrixBasedRhsAssembly);
//    VecAssemblyEnd(this->mVectorForMatrixBasedRhsAssembly);
//    VecDestroy(force_term_at_nodes);
//}







template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>::ConstructVectorForMatrixBasedRhsAssembly(
        Vec existingSolution)
{
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();

    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(existingSolution);

    // dist stripe for z (return value)
    DistributedVector dist_vec_matrix_based = p_factory->CreateDistributedVector(this->mVectorForMatrixBasedRhsAssembly);

    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm  = HeartConfig::Instance()->GetCapacitance();

    for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
         index!= dist_vec_matrix_based.End();
         ++index)
    {
        double V = distributed_current_solution[index];
        double F = - Am*this->mpMonodomainPde->rGetIionicCacheReplicated()[index.Global]
                   - this->mpMonodomainPde->rGetIntracellularStimulusCacheReplicated()[index.Global];

        dist_vec_matrix_based[index] = Am*Cm*V*this->mDtInverse + F;
    }

    dist_vec_matrix_based.Restore();
}



/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class MonodomainMatrixBasedAssembler<1,1>;
template class MonodomainMatrixBasedAssembler<1,2>;
template class MonodomainMatrixBasedAssembler<1,3>;
template class MonodomainMatrixBasedAssembler<2,2>;
template class MonodomainMatrixBasedAssembler<3,3>;

template class MonodomainRhsMatrixAssembler<1,1>;
template class MonodomainRhsMatrixAssembler<1,2>;
template class MonodomainRhsMatrixAssembler<1,3>;
template class MonodomainRhsMatrixAssembler<2,2>;
template class MonodomainRhsMatrixAssembler<3,3>;
