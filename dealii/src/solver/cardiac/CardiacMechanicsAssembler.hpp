#ifndef CARDIACMECHANICSASSEMBLER_HPP_
#define CARDIACMECHANICSASSEMBLER_HPP_
    
#include "FiniteElasticityAssembler.cpp"

template<unsigned DIM>
class CardiacMechanicsAssembler : public FiniteElasticityAssembler<DIM>
{
private:
    /** The active tensions, node-wise */
    std::vector<double> mActiveTensions;

    /** Overloaded method for assembling system, which takes into account the active tensions */
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                           Vector<double>&       elementRhs,
                           FullMatrix<double>&   elementMatrix,
                           bool                  assembleResidual,
                           bool                  assembleJacobian);

public:
    CardiacMechanicsAssembler(Triangulation<DIM>* pMesh, std::string outputDirectory);
    ~CardiacMechanicsAssembler();
    
    /** Set the current active tensions, node-wise */
    void SetActiveTensions(std::vector<double> activeTensions);
};

#endif /*CARDIACMECHANICSASSEMBLER_HPP_*/
