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

#ifndef ABSTRACTMONODOMAINSOLVER_HPP_
#define ABSTRACTMONODOMAINSOLVER_HPP_


#include "MonodomainTissue.hpp"
#include "MonodomainAssembler.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"

/**
 *  Abstract Monodomain class containing some common functionality
 *  Inherits from AbstractDynamicLinearPdeSolver so child classes must
 *  implement SetupLinearSystem()
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractMonodomainSolver 
   : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>
{
protected:
    /** Boundary conditions */    
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* mpBoundaryConditions;
   
    /** Monodomain tissue class (collection of cells, and conductivities) */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* mpMonodomainTissue;
    
    /**
     *  The monodomain assembler, used to set up the LHS matrix;
     *  and RHS vector for BasicMonodomainSolver (but not MatrixBasedMonodomainSolver)
     */
    MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainAssembler;
    
    /**
     *  Number of quadrature points per dimension (only saved so it can be
     *  passed to the assembler
     */
    unsigned mNumQuadPoints;        

public:
    /**
     *  Overloaded PrepareForSetupLinearSystem() methods which
     *  gets the cell models to solve themselves
     * 
     *  @param currentSolution solution at current time
     */
    void PrepareForSetupLinearSystem(Vec currentSolution);

    /**
     *  Overloaded InitialiseForSolve
     * 
     *  @param initialSolution initial solution
     */
    virtual void InitialiseForSolve(Vec initialSolution);


    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    AbstractMonodomainSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                             MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                             BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                             unsigned numQuadPoints = 2);
    
    /** 
     *  Destructor - just deletes assembler 
     */
    virtual ~AbstractMonodomainSolver();
};    

#endif /*ABSTRACTMONODOMAINSOLVER_HPP_*/
