#ifndef CARDIACMECHASSEMBLER_HPP_
#define CARDIACMECHASSEMBLER_HPP_

#include "FiniteElasticityAssembler.cpp"

template<unsigned DIM>
class CardiacMechAssembler : public FiniteElasticityAssembler<DIM>
{
private:
    bool mAllocatedMaterialLawMemory;
    
    /** The active tension, node-wise */
    std::vector<double> mActiveTension;

    /** Overloaded method for assembling system, which takes into account the active tensions */
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                           Vector<double>&       elementRhs,
                           FullMatrix<double>&   elementMatrix,
                           bool                  assembleResidual,
                           bool                  assembleJacobian);

public:
    /**
     *  Constructor
     *  
     *  @param pMesh. A pointer to the mesh. Should have a surface set as the fixed surface
     *  @param outputDirectory. The output directory, relative to TEST_OUTPUT
     *  @param pMaterialLaw. The material law for the tissue. Defaults to NULL, in which case
     *   a default material law is used.
     */
    CardiacMechAssembler(Triangulation<DIM>* pMesh, 
                              std::string outputDirectory,
                              AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw = NULL);
    ~CardiacMechAssembler();
    
    /** Set the current active tensions, node-wise */
    void SetActiveTension(std::vector<double> activeTension);
    
    /** 
     *  Get lambda (the stretch ratio). 
     *  @param lambda A std::vector which will be filled in. Does not have to be correctly-sized
     */
    void GetLambda(std::vector<double>& lambda);    
};

#endif /*CARDIACMECHASSEMBLER_HPP_*/
