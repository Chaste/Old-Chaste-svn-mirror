#ifndef ABSTRACTONEDIMCARDIACMECHANICSASSEMBLER_HPP_
#define ABSTRACTONEDIMCARDIACMECHANICSASSEMBLER_HPP_

#include "AbstractElasticityAssembler.hpp"
#include "PoleZero3dIn1dLaw.hpp"

class AbstractOneDimCardiacMechanicsAssembler : public AbstractElasticityAssembler<1>
{
protected :
    FE_Q<1> mFe;
    double mDensity;
    unsigned mNumNewtonIterations;

    // cardiac mech
    static const unsigned mNumQuadPointsInEachDimension = 3;
    unsigned mTotalQuadPoints;
    unsigned mCurrentQuadPointGlobalIndex;

    unsigned mEndNodeDof;
    PoleZero3dIn1dLaw mLaw;

    void ApplyDirichletBoundaryConditions()
    {
        for(unsigned j=0; j<mSystemMatrix.n(); j++)
        {
            this->mSystemMatrix.set(mEndNodeDof, j, 0.0);
        }
        this->mSystemMatrix.set(mEndNodeDof, mEndNodeDof, 1.0);
        this->mRhsVector(mEndNodeDof) = this->mCurrentSolution(mEndNodeDof);
    }

    void DistributeDofs()
    {
        this->mDofHandler.distribute_dofs(mFe);
    }
    
    std::vector<double> mLambda;

public:
    AbstractOneDimCardiacMechanicsAssembler(Triangulation<1>* pMesh)
        : AbstractElasticityAssembler<1>(pMesh),
          mFe(1)
    {
        DistributeDofs();
        InitialiseMatricesVectorsAndConstraints();
        mDofsPerElement = mFe.dofs_per_cell;
        
        mDensity = 1.0;
        mNumNewtonIterations = 0;
        
        // set up quad point info
        QGauss<1>   quadrature_formula(mNumQuadPointsInEachDimension);
        mTotalQuadPoints = quadrature_formula.n_quadrature_points *this->mpMesh->n_active_cells();
        mCurrentQuadPointGlobalIndex = 0;
        
        bool found = false;
        DofVertexIterator<1> vertex_iter(this->mpMesh, &this->mDofHandler);
        while (!vertex_iter.ReachedEnd())
        {
            Point<1> posn = vertex_iter.GetVertex();
            if( posn[0]==0)
            {
                mEndNodeDof = vertex_iter.GetDof(0);
                found = true;
                break;
            }
            vertex_iter.Next();
        }
        assert(found); // check have found the end node..
        
        mLaw.SetUpStores();
        
        mLambda.resize(mTotalQuadPoints,1.0);
    }    


    virtual void Solve(double currentTime, double nextTime, double timestep) //params are a bit of a hack for refactoring at the moment..
    {
        // compute residual
        this->AssembleSystem(true, false);
        double norm_resid = this->CalculateResidualNorm();
        std::cout << "\nNorm of residual is " << norm_resid << "\n";
        
        mNumNewtonIterations = 0;
        unsigned counter = 1;
    
        // use the larger of the tolerances formed from the absolute or
        // relative possibilities
        double tol = NEWTON_ABS_TOL;
        if ( tol < NEWTON_REL_TOL*norm_resid )
        {
            tol = NEWTON_REL_TOL*norm_resid;
        }
        std::cout << "Solving with tolerance " << tol << "\n";
        
        while (norm_resid > tol)
        {
            std::cout <<  "\n-------------------\n"
                      <<   "Newton iteration " << counter
                      << ":\n-------------------\n";
            
            this->TakeNewtonStep();
            this->AssembleSystem(true, false);
            norm_resid = this->CalculateResidualNorm();
            
            std::cout << "Norm of residual is " << norm_resid << "\n";
            
            //WriteOutput(counter);
            mNumNewtonIterations = counter;
            
            counter++;
            if (counter==20)
            {
                EXCEPTION("Not converged after 20 newton iterations, quitting");
            }
        }
    
        if (norm_resid > tol)
        {
            EXCEPTION("Failed to converge");
        }
    }

    unsigned GetTotalNumQuadPoints()
    {
        return mTotalQuadPoints;
    }
    
    unsigned GetNumQuadPointsPerElement()
    {
        return mNumQuadPointsInEachDimension; 
    }
    
    unsigned GetNumQuadPointsInEachDimension()
    {
        return mNumQuadPointsInEachDimension;
    }
    
    virtual void SetForcingQuantity(std::vector<double>& forcingQuantity)=0;

    std::vector<double>& GetLambda()
    {
        return mLambda;
    }       
};

#endif /*ABSTRACTONEDIMCARDIACMECHANICSASSEMBLER_HPP_*/
