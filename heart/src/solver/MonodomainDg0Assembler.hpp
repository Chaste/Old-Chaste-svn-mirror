#ifndef _MONODOMAINDG0ASSEMBLER_HPP_
#define _MONODOMAINDG0ASSEMBLER_HPP_


//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "ConformingTetrahedralMesh.cpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "GaussianQuadratureRule.hpp"
#include "MonodomainPde.hpp"


/**
 *  MonodomainDg0Assembler
 *
 *  This is essentially the same as the SimpleDg0ParabolicAssembler (which it inherits from),
 *  except that the source term (ie ionic current + stimulus) is interpolated from
 *  their nodal values, instead of computed at the gauss point, since they are only
 *  known at the nodes.
 *
 *  Also, the MonodomainAssembler automatically creates zero neumann boundary conditions
 *  when constructed and therefore does not need to take in a BoundaryConditionsContainer.
 *
 *  The user should call Solve() from the superclass AbstractDynamicAssemblerMixin.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainDg0Assembler
    : public SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, false, MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> >
{
public:
    static const unsigned E_DIM = ELEMENT_DIM;
    static const unsigned S_DIM = SPACE_DIM;
    static const unsigned P_DIM = 1u;
private:
    double mSourceTerm;
    
    MonodomainPde<SPACE_DIM>* mpMonodomainPde;

    // Save typing
    typedef MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> SelfType;
    typedef SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, false, SelfType> BaseClassType;
    
    /// Allow the AbstractStaticAssembler to call our private/protected methods using static polymorphism.
    friend class AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, 1u, false, BaseClassType>;
    
protected:

    /**
     *  ComputeVectorTerm()
     * 
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness vector.
     * 
     *  Here, the SimpleDg0ParabolicAssembler version of this method is 
     *  overloaded using the interpolated source term
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        return  rPhi * (mSourceTerm + this->mDtInverse *
                        mpMonodomainPde->ComputeDuDtCoefficientFunction(rX) * u(0));
    }
    
    
    void ResetInterpolatedQuantities( void )
    {
        mSourceTerm=0;
    }
    
    
    void IncrementInterpolatedQuantities(double phi_i, const Node<SPACE_DIM> *pNode)
    {
        mSourceTerm += phi_i * mpMonodomainPde->ComputeNonlinearSourceTermAtNode(*pNode, this->mCurrentSolutionOrGuessReplicated[ pNode->GetIndex() ] );
    }
    
    
    virtual void PrepareForAssembleSystem(Vec currentSolution, double currentTime)
    {
        mpMonodomainPde->SolveCellSystems(currentSolution, currentTime, currentTime+this->mDt);
    }
    
    
    
public:
    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     */
    MonodomainDg0Assembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                           MonodomainPde<SPACE_DIM>* pPde,
                           BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBcc,
                           unsigned numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
            BaseClassType(pMesh, pPde, NULL /*bcs - set below*/, numQuadPoints)
    {
        mpMonodomainPde = pPde;
        
        this->mpBoundaryConditions = pBcc;
        
        this->SetMesh(pMesh);
        
        this->SetMatrixIsConstant();
    }
    
    ~MonodomainDg0Assembler()
    {
    }
};

/**
 * Specialization of AssemblerTraits for the MonodomainDg0Assembler.
 *
 * This is always a concrete class, but only defines some of the methods.
 * For others it thus has to know which base class defines them.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct AssemblerTraits<MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> >
{
    typedef MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CVT_CLS;
    typedef SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, false, MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> >
            CMT_CLS;
    typedef MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> INTERPOLATE_CLS;
};

#endif //_MONODOMAINDG0ASSEMBLER_HPP_
