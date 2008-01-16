#ifndef PARABOLICFLAGGEDMESHASSEMBLER_HPP_
#define PARABOLICFLAGGEDMESHASSEMBLER_HPP_

#include "SimpleDg0ParabolicAssembler.hpp"
#include "AbstractFlaggedMeshAssemblerMixin.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"

/**
 * A simple flagged mesh assembler for parabolic PDEs.
 */
template<unsigned DIM>
class ParabolicFlaggedMeshAssembler
    : public SimpleDg0ParabolicAssembler<DIM, DIM, true, ParabolicFlaggedMeshAssembler<DIM> >,
      public AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>
{
public:
    static const unsigned E_DIM = DIM;
    static const unsigned S_DIM = DIM;
    static const unsigned P_DIM = 1u;

    typedef SimpleDg0ParabolicAssembler<DIM, DIM, true, ParabolicFlaggedMeshAssembler<DIM> > BaseClassType;

private:
    friend class TestFlaggedMeshAssembler;
    
    /**
     * This needs to be implemented so the compiler knows which subclass version of the
     * method to call, in this case the flagged version.
     * 
     * This is a minor pain, but there's no way around it.
     */
    void AssembleSystem(bool assembleVector, bool assembleMatrix, Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>::AssembleSystem(assembleVector, assembleMatrix, currentSolutionOrGuess, currentTime);
    }
    
    double GetCurrentSolutionOrGuessValue(unsigned nodeIndex, unsigned indexOfUnknown)
    {      
        return AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>::GetCurrentSolutionOrGuessValue(nodeIndex, indexOfUnknown);
    }
  
public :
    ParabolicFlaggedMeshAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                                  AbstractLinearParabolicPde<DIM>* pPde,
                                  FlaggedMeshBoundaryConditionsContainer<DIM,1>* pBoundaryConditions,
                                  unsigned numQuadPoints = 2) :
            AbstractAssembler<DIM,DIM,1>(),
            BaseClassType(pMesh,pPde,NULL,numQuadPoints),
            AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>(pBoundaryConditions)
    {
        this->SetMatrixIsConstant(false);
    }
};

template<unsigned DIM>
struct AssemblerTraits<ParabolicFlaggedMeshAssembler<DIM> >
{
    typedef typename ParabolicFlaggedMeshAssembler<DIM>::BaseClassType CVT_CLS;
    typedef typename ParabolicFlaggedMeshAssembler<DIM>::BaseClassType CMT_CLS;
    typedef AbstractAssembler<DIM, DIM, 1u> INTERPOLATE_CLS;
};

#endif /*PARABOLICFLAGGEDMESHASSEMBLER_HPP_*/
