#ifndef ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_
#define ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_

#include <vector>
#include "AbstractDealiiAssembler.hpp" // just to include some dealii classes


template<unsigned DIM>
class AbstractCardiacMechanicsAssembler
{
protected :
    static const unsigned mNumQuadPointsInEachDimension = 3;

    unsigned mTotalQuadPoints;
    unsigned mCurrentQuadPointGlobalIndex;    

    /** 
     *  The x stretch, quad point-wise. NOTE: the i-th entry of this vector is
     *  assumed to be the i-th quad point obtained by looping over cells in the obvious
     *  way and then looping over quad points 
     */
    std::vector<double> mLambda;

public :
    AbstractCardiacMechanicsAssembler(Triangulation<DIM>* pMesh)
    {
        assert(pMesh!=NULL);

        // set up quad point info
        QGauss<DIM>   quadrature_formula(mNumQuadPointsInEachDimension);
        mTotalQuadPoints = quadrature_formula.n_quadrature_points * pMesh->n_active_cells();
        mCurrentQuadPointGlobalIndex = 0;
        
        mLambda.resize(mTotalQuadPoints, 1.0);
    }
    
    virtual ~AbstractCardiacMechanicsAssembler()
    {
    }


    unsigned GetNumQuadPointsInEachDimension()
    {
        return mNumQuadPointsInEachDimension;
    }
    
    /** 
     *  Get the total number of quadrature points (equal to the n.q^d, where n=number of cells
     *  q is the number of quad points in each dimension, and d=DIM). 
     */
    unsigned GetTotalNumQuadPoints()
    {
        return mTotalQuadPoints;
    }

    /**
     *  Get the total number of quadrature points in each element (ie num_quad_in_each_dir^DIM)
     */
    unsigned GetNumQuadPointsPerElement()
    {
        if(DIM==1)
        {
            return mNumQuadPointsInEachDimension;
        }
        if(DIM==2)
        { 
            return mNumQuadPointsInEachDimension*mNumQuadPointsInEachDimension;
        }
        else //DIM==3
        {
            return mNumQuadPointsInEachDimension*mNumQuadPointsInEachDimension*mNumQuadPointsInEachDimension;
        }
    }
    
        
    /** 
     *  Get lambda (the stretch ratio). 
     * 
     *  NOTE: the i-th entry of this vector is
     *  assumed to be the i-th quad point obtained by looping over cells in the obvious
     *  way and then looping over quad points 
     */
    std::vector<double>& rGetLambda()
    {
        return mLambda;
    }       

    virtual void Solve(double currentTime, double nextTime, double timestep)=0;

    virtual void SetForcingQuantity(std::vector<double>& forcingQuantity)=0;
};


#endif /*ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_*/
