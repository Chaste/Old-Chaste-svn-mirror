#ifndef CARDIACMECHANICSASSEMBLER_HPP_
#define CARDIACMECHANICSASSEMBLER_HPP_

#include "FiniteElasticityAssembler.cpp"
#include "AbstractCardiacMechanicsAssembler.hpp"

template<unsigned DIM>
class CardiacMechanicsAssembler : public FiniteElasticityAssembler<DIM>, public AbstractCardiacMechanicsAssembler<DIM>
{
friend class TestImplicitCardiacMechanicsAssembler;

protected:
    bool mAllocatedMaterialLawMemory;

//    virtual double GetActiveTensionAtCurrentQuadPoint(double lam)
//    {
//        return mActiveTension[this->mCurrentQuadPointGlobalIndex];
//    }

    /**
     *  Storage space for dTdE when T and E are in the rotated fibre-sheet frame
     */
    FourthOrderTensor<DIM>  mDTdE_fibre;

    /**
     *  The matrix P using JonW's convention. Orthogonal
     */
    Tensor<2,DIM> mFibreSheetMat;
    
    /**
     *  The transpose of P, which is also the inverse of P
     */
    Tensor<2,DIM> mTransFibreSheetMat;

    /** 
     *  The active tension, quad point-wise. NOTE: the i-th entry of this vector is
     *  assumed to be the i-th quad point obtained by looping over cells in the obvious
     *  way and then looping over quad points 
     */
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
    CardiacMechanicsAssembler(Triangulation<DIM>* pMesh, 
                              std::string outputDirectory,
                              AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw = NULL);
    virtual ~CardiacMechanicsAssembler();
    
        
    /**
     *  Specify a constant fibre-sheet rotation matrix
     * 
     *  This is really a temporary method until the fibre-sheet direction can be read in
     */
    virtual void SetFibreSheetMatrix(Tensor<2,DIM> fibreSheetMat)
    {
        // check orthogonal
        Tensor<2,DIM> P_times_transP = fibreSheetMat * transpose(fibreSheetMat);
        for(unsigned i=0; i<DIM; i++)
        {
            for(unsigned j=0; j<DIM; j++)
            {
                double expected = i==j ? 1.0 : 0.0;
                if (fabs(P_times_transP[i][j] - expected) > 1e-9)
                {
                    EXCEPTION("Fibre-sheet matrix passed in does not seem to be orthogonal");
                }
            }
        }
        
        mFibreSheetMat = fibreSheetMat;
        mTransFibreSheetMat = transpose(mFibreSheetMat);
    }
    
    virtual void Solve(double currentTime, double nextTime, double timestep);

    
    /** 
     *  Set the current active tensions, by quadrature point. Quad points don't have indices,
     *  so these values should be in the order given by looping over cells and then looping
     *  over quad points
     */
    virtual void SetForcingQuantity(std::vector<double>& activeTension);
};

#endif /*CARDIACMECHANICSASSEMBLER_HPP_*/
