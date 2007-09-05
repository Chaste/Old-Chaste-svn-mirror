#ifndef ELLIPTICFLAGGEDMESHASSEMBLER_HPP_
#define ELLIPTICFLAGGEDMESHASSEMBLER_HPP_

#include "SimpleLinearEllipticAssembler.hpp"
#include "AbstractFlaggedMeshAssemblerMixin.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"

/**
 * A simple flagged mesh assembler for elliptic PDEs.
 */
template<unsigned DIM>
class EllipticFlaggedMeshAssembler : public SimpleLinearEllipticAssembler<DIM,DIM>, public AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>
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
        
    double GetCurrentSolutionOrGuessValue(unsigned nodeIndex, unsigned indexOfUnknown)
    {      
        return AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>::GetCurrentSolutionOrGuessValue(nodeIndex, indexOfUnknown);
    }
public :
    EllipticFlaggedMeshAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                                  AbstractLinearEllipticPde<DIM>* pPde,
                                  FlaggedMeshBoundaryConditionsContainer<DIM,1>* pBoundaryConditions,
                                  unsigned numQuadPoints = 2) :
            AbstractAssembler<DIM,DIM,1>(),
            SimpleLinearEllipticAssembler<DIM,DIM>(pMesh,pPde,NULL,numQuadPoints),
            AbstractFlaggedMeshAssemblerMixin<DIM,DIM,1>(pBoundaryConditions)
    {
    }
};
#endif /*ELLIPTICFLAGGEDMESHASSEMBLER_HPP_*/
