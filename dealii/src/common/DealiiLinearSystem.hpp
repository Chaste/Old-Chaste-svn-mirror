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

#ifndef DEALIILINEARSYSTEM_HPP_
#define DEALIILINEARSYSTEM_HPP_

#include "AbstractDealiiAssembler.hpp"
#include "UblasCustomFunctions.hpp"
#include "QuadraticMesh.hpp"

class DealiiLinearSystem
{
friend class TestDealiiLinearSystem;
    
private:
    SparsityPattern mSparsityPattern;
    SparseMatrix<double>  mLhsMatrix;
    Vector<double> mRhsVector;
    Vector<double> mLhsVector;

public:
    DealiiLinearSystem(unsigned size)
        : mSparsityPattern(size,size), // second param is num_nonzeros!!
          mRhsVector(size)
    {
        assert(size>0);

        for(unsigned i=0; i<size; i++)
        {
            for(unsigned j=0; j<size; j++)
            {
                mSparsityPattern.add(i,j);
            }
        }

        mSparsityPattern.compress();

        mLhsMatrix.reinit(mSparsityPattern);
        ZeroRhsVector();
        ZeroLhsMatrix();
    }
    
    template<unsigned DIM>
    DealiiLinearSystem(QuadraticMesh<DIM>& rQuadMesh)
        : mSparsityPattern(DIM*rQuadMesh.GetNumNodes() + rQuadMesh.GetNumVertices(), 200)
    {
        unsigned num_nodes = rQuadMesh.GetNumNodes(); 
        unsigned num_dofs = DIM*num_nodes + rQuadMesh.GetNumVertices(); 
        mRhsVector.reinit(num_dofs);
        
        for(unsigned i=0; i<rQuadMesh.GetNumElements(); i++)
        {
            Element<DIM,DIM>& r_elem = *(rQuadMesh.GetElement(i));
            for(unsigned j=0; j<r_elem.GetNumNodes(); j++)
            {
                bool pressure_node1 = (r_elem.GetNodeGlobalIndex(j) < rQuadMesh.GetNumVertices());

                unsigned last1 = pressure_node1 ? DIM+1 : DIM;
                for(unsigned dim1=0; dim1<last1; dim1++)
                {
                    unsigned index1 = DIM*r_elem.GetNodeGlobalIndex(j)  + dim1;
                    if(dim1==DIM) //pressure
                    {
                        index1 = DIM*num_nodes + r_elem.GetNodeGlobalIndex(j);
                    }
                    assert(index1<num_dofs);
                    for(unsigned k=0; k<r_elem.GetNumNodes(); k++)
                    {
                        bool pressure_node2 = (r_elem.GetNodeGlobalIndex(k) < rQuadMesh.GetNumVertices());
        
                        unsigned last2 = pressure_node2 ? DIM+1 : DIM;
                        for(unsigned dim2=0; dim2<last2; dim2++)
                        {
                            unsigned index2 = DIM*r_elem.GetNodeGlobalIndex(k)  + dim2;
                            if(dim2==DIM) //pressure
                            {
                                index2 = DIM*num_nodes + r_elem.GetNodeGlobalIndex(k);
                            }
                            assert(index2<num_dofs);
                          //  std::cout << i << ": " << r_elem.GetNodeGlobalIndex(j) << ", " << r_elem.GetNodeGlobalIndex(k)
                            //          << ": adding " << index1 << " and " << index2 << "\n";
                            mSparsityPattern.add(index1,index2);
                        }
                    }
                }
            }
        }
                        
        mSparsityPattern.compress();

        mLhsMatrix.reinit(mSparsityPattern);
        
        
        ZeroRhsVector();
        ZeroLhsMatrix();
    }
    
    void ZeroRhsVector()
    {
        mRhsVector = 0;
    }

    void ZeroLhsMatrix()
    {
        mLhsMatrix = 0;
    }
    
    void Solve()
    {
        mLhsVector.reinit(mLhsMatrix.m());
        SparseDirectUMFPACK direct_solver;
        direct_solver.initialize(mLhsMatrix);
        direct_solver.vmult(mLhsVector, mRhsVector);
    }
    
    Vector<double>& rGetLhsVector()
    {
        return mLhsVector;
    }
    
    void ZeroMatrixRow(int row)
    {
//// not sure how to use the iterator......
//        for(SparseMatrix<double>::iterator iter = mLhsMatrix.begin(row);
//            iter != mLhsMatrix.end(row);
//            ++iter)
//        {
//            mLhsMatrix.set(row,*iter,0.0);
//        }

        for(unsigned j=0; j<mLhsMatrix.n(); j++)
        {
            mLhsMatrix.set(row,j,0.0);
        }
    }

    void SetMatrixElement(unsigned index1, unsigned index2, double value)
    {
        mLhsMatrix.set(index1, index2, value);
    }
    
    void SetRhsVectorElement(unsigned index, double value)
    {
        mRhsVector(index) = value;
    }
    
    double GetRhsVectorNorm()
    {
        return mRhsVector.l2_norm();
    }
    
    template<size_t SIZE>
    void AddLhsMultipleValues(unsigned pIndices[SIZE], c_matrix<double,SIZE,SIZE> aElem)
    {
        for(unsigned i=0; i<SIZE; i++)
        {
            unsigned index1 = pIndices[i];
            assert(index1<mLhsMatrix.m());
            for(unsigned j=0; j<SIZE; j++)
            {
                unsigned index2 = pIndices[j];
                assert(index2<mLhsMatrix.m());
                mLhsMatrix.add(index1,index2,aElem(i,j));
            }
        } 
    }

    template<size_t SIZE>
    void AddRhsMultipleValues(unsigned pIndices[SIZE], c_vector<double,SIZE> bElem)
    {
        for(unsigned i=0; i<SIZE; i++)
        {
            unsigned index = pIndices[i];
            assert(index<mLhsMatrix.m());
            
            mRhsVector(index) += bElem(i);
        }
    }
};
#endif /*DEALIILINEARSYSTEM_HPP_*/
