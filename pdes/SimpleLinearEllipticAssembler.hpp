#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
//#include "AbstractBoundaryConditions.hpp"
#include <vector>
#include "petscvec.h"
#include "AbstractLinearSolver.hpp"

#include <iostream>


template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleLinearEllipticAssembler : public AbstractLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{
    
private:
    LinearSystem *mpAssembledLinearSystem;

public:
    Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractLinearEllipticPde<SPACE_DIM> *pPde, 
//                       BoundaryConditionsContainer &rBoundaryConditions,
                       AbstractLinearSolver *solver)
	{
		// Linear system in n unknowns, where n=#nodes
        mpAssembledLinearSystem	= new LinearSystem(rMesh.GetNumNodes());
        
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter = rMesh.GetFirstElement();
        while (iter != rMesh.GetLastElement())
        {
            Element<ELEMENT_DIM, SPACE_DIM> element = *iter;
            // This assumes linear basis functions in 1d
            int node1 = element.GetNodeGlobalIndex(0);
            int node2 = element.GetNodeGlobalIndex(1);
            
            double x1 = element.GetNodeLocation(0,0);
            double x2 = element.GetNodeLocation(1,0);
            double h = x2-x1;
            double one_over_h_squared = 1.0/h/h;
            double integral_of_d = h; // Will depend on pPde->ComputeDiffusionTerm
            
            mpAssembledLinearSystem->AddToMatrixElement(node1,node1, one_over_h_squared*integral_of_d);
            mpAssembledLinearSystem->AddToMatrixElement(node1,node2,-one_over_h_squared*integral_of_d);
            mpAssembledLinearSystem->AddToMatrixElement(node2,node1,-one_over_h_squared*integral_of_d);
            mpAssembledLinearSystem->AddToMatrixElement(node2,node2, one_over_h_squared*integral_of_d);
            
            mpAssembledLinearSystem->AssembleIntermediateMatrix();  
            
            // Will depend on pPde->Compute(Linear|Nonlinear)SourceTerm
            mpAssembledLinearSystem->AddToRhsVectorElement(node1,0.5*h);
            mpAssembledLinearSystem->AddToRhsVectorElement(node2,0.5*h);      
         
            iter++;
        }
        
//        for(int i=0; i<rBoundaryConditions.size(); i++)
//        {
//            rBoundaryConditions[i]->ApplyLinearBoundaryConditions(*mpAssembledLinearSystem);   
//        }
		mpAssembledLinearSystem->SetMatrixElement(0, 0, 1.0);
    	mpAssembledLinearSystem->SetMatrixElement(0, 1, 0.0);
    	mpAssembledLinearSystem->SetRhsVectorElement(0, 0.0);
    
        mpAssembledLinearSystem->AssembleFinalMatrix();
        
        Vec sol = mpAssembledLinearSystem->Solve(solver);       
        return sol;
	}
};



#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
