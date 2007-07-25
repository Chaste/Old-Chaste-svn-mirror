#ifndef PARABOLICFLAGGEDMESHASSEMBLER_HPP_
#define PARABOLICFLAGGEDMESHASSEMBLER_HPP_

#include "SimpleDg0ParabolicAssembler.hpp"
#include "AbstractFlaggedMeshAssemblerMixin.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"

/**
 * A simple flagged mesh assembler for parabolic PDEs.
 */
template<unsigned DIM>
class ParabolicFlaggedMeshAssembler : public SimpleDg0ParabolicAssembler<DIM,DIM>, public AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>
{
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
    
public :
    ParabolicFlaggedMeshAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                                  AbstractLinearParabolicPde<DIM>* pPde,
                                  FlaggedMeshBoundaryConditionsContainer<DIM,1>* pBoundaryConditions,
                                  unsigned numQuadPoints = 2) :
            AbstractAssembler<DIM,DIM,1>(),
            SimpleDg0ParabolicAssembler<DIM,DIM>(pMesh,pPde,NULL,numQuadPoints),
            AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>(pBoundaryConditions)
    {
    }
};
#endif /*PARABOLICFLAGGEDMESHASSEMBLER_HPP_*/
