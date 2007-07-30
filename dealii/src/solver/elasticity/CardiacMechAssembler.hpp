#ifndef CARDIACMECHASSEMBLER_HPP_
#define CARDIACMECHASSEMBLER_HPP_

#include "FiniteElasticityAssembler.cpp"

template<unsigned DIM>
class CardiacMechAssembler : public FiniteElasticityAssembler<DIM>
{
private:
    static const unsigned mNumQuadPointsInEachDimension = 3;
    unsigned mTotalQuadPoints;
    unsigned mCurrentQuadPointGlobalIndex;

    bool mAllocatedMaterialLawMemory;

    /** The active tension, node-wise */
    std::vector<double> mActiveTension;
    std::vector<double> mLambda;
    
    /** Overloaded method for assembling system, which takes into account the active tensions */
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                           Vector<double>&       elementRhs,
                           FullMatrix<double>&   elementMatrix,
                           bool                  assembleResidual,
                           bool                  assembleJacobian);
                           
    void SetUpQuadraturePointsInfo();                           

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
    
    /** 
     *  Get the total number of quadrature points (equal to the n.q^d, where n=number of cells
     *  q is the number of quad points in each dimension, and d=DIM). 
     */
    unsigned GetTotalNumQuadPoints();
    
    /** 
     *  Set the current active tensions, by quadrature point. Quad points don't have indices,
     *  so these values should be in the order given by looping over cells and then looping
     *  over quad points
     */
    void SetActiveTension(std::vector<double> activeTension);
    
    /** 
     *  Get lambda (the stretch ratio). 
     */
    std::vector<double>& GetLambda();    
};

#endif /*CARDIACMECHASSEMBLER_HPP_*/
