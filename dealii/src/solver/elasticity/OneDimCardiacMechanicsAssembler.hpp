#ifndef _ONEDIMCARDIACMECHANICSASSEMBLER_HPP_
#define _ONEDIMCARDIACMECHANICSASSEMBLER_HPP_

#include "AbstractElasticityAssembler.hpp"
#include "PoleZero3dIn1dLaw.hpp"

class OneDimCardiacMechanicsAssembler : public AbstractElasticityAssembler<1>
{
protected :
    FE_Q<1> mFe;
    double mDensity;
    unsigned mNumNewtonIterations;

    // cardiac mech
    static const unsigned mNumQuadPointsInEachDimension = 3;
    unsigned mTotalQuadPoints;
    unsigned mCurrentQuadPointGlobalIndex;
    std::vector<double> mActiveTension;
    std::vector<double> mLambda;

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
    

public:
    OneDimCardiacMechanicsAssembler(Triangulation<1>* pMesh)
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
        
        mLambda.resize(mTotalQuadPoints, 1.0);

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
    }    


    virtual void Solve(double startTime, double endTime, double timestep) //params are a bit of a hack for refactoring at the moment..
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
    
        // we have solved for a deformation so note this
        //mADeformedHasBeenSolved = true;
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

    virtual void SetForcingQuantity(std::vector<double>& activeTension)
    {
        assert(activeTension.size() == mTotalQuadPoints);
        mActiveTension = activeTension;
    }

    std::vector<double>& GetLambda()
    {
        return mLambda;
    }   

private:

    void AssembleOnElement(DoFHandler<1>::active_cell_iterator  elementIter,
                           Vector<double>&        elementRhs,
                           FullMatrix<double>&    elementMatrix,
                           bool                   assembleResidual,
                           bool                   assembleJacobian)
    {
        // if mCurrentQuadPointGlobalIndex is greater than the total num of quad points something
        // very bad has happened. 
        assert(mCurrentQuadPointGlobalIndex <= mTotalQuadPoints);
        
        if(mCurrentQuadPointGlobalIndex==mTotalQuadPoints)
        {
            // if we are not back to the first cell something bad has happened
            assert( elementIter == this->mDofHandler.begin_active() );
            
            mCurrentQuadPointGlobalIndex = 0;
        }
    
    
        static QGauss<1>   quadrature_formula(mNumQuadPointsInEachDimension);
        const unsigned n_q_points    = quadrature_formula.n_quadrature_points;
        
        
        // would want this to be static too (slight speed up), but causes errors
        // in debug mode (upon destruction of the class, in 2d, or something)
        FEValues<1> fe_values(mFe, quadrature_formula,
                              UpdateFlags(update_values    |
                                          update_gradients |
                                          update_q_points  |     // needed for interpolating u and u' on the quad point
                                          update_JxW_values));

        const unsigned dofs_per_element = mFe.dofs_per_cell;
        
        static std::vector< Vector<double> >                  local_solution_values(n_q_points);
        static std::vector< std::vector< Tensor<1,1> > >    local_solution_gradients(n_q_points);
        
        static Tensor<2,1> identity;
        
        static bool first = true;
        
        
        if (first)
        {
            for (unsigned q_point=0; q_point<n_q_points; q_point++)
            {
                local_solution_values[q_point].reinit(1+1);
                local_solution_gradients[q_point].resize(1+1);
            }
            for (unsigned i=0; i<1; i++)
            {
                for (unsigned j=0; j<1; j++)
                {
                    identity[i][j] = i==j ? 1.0 : 0.0;
                }
            }
        }
            
        elementMatrix = 0;
        elementRhs = 0;
    
        fe_values.reinit(elementIter); // compute fe values for this element
        fe_values.get_function_values(this->mCurrentSolution, local_solution_values);
        fe_values.get_function_grads(this->mCurrentSolution, local_solution_gradients);

        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            const std::vector< Tensor<1,1> >& grad_u = local_solution_gradients[q_point];
            
            double F = 1 + grad_u[0][0];;
            double C = F*F;
            double T;
            double dTdE;            

            // get the active tension at this quad point
            double active_tension = mActiveTension[mCurrentQuadPointGlobalIndex];
            mLambda[mCurrentQuadPointGlobalIndex] = F; 

            // Compute T(C), dTdE(C)...
            double E = 0.5*(C-1);
            T = mLaw.GetT(E) + active_tension/C; 
            dTdE = mLaw.GetDTdE(E) - 2*active_tension/(C*C);
    

            for (unsigned i=0; i<dofs_per_element; i++)
            {
                if (assembleJacobian)
                {
                    for (unsigned j=0; j<dofs_per_element; j++)
                    {
                        elementMatrix(i,j) +=   T                                 
                                              * fe_values.shape_grad(j,q_point)[0]
                                              * fe_values.shape_grad(i,q_point)[0]
                                              * fe_values.JxW(q_point);

                        elementMatrix(i,j) +=   0.5
                                              * dTdE 
                                              * (
                                                    fe_values.shape_grad(j,q_point)[0]
                                                  * F  
                                                  +
                                                    fe_values.shape_grad(j,q_point)[0]
                                                  * F  
                                                )
                                              * F   
                                              * fe_values.shape_grad(i,q_point)[0]
                                              * fe_values.JxW(q_point);
                    }
                }
                
                if (assembleResidual)
                {
                    elementRhs(i) +=   T
                                     * F
                                     * fe_values.shape_grad(i,q_point)[0]
                                     * fe_values.JxW(q_point);
                }
            }
    
            mCurrentQuadPointGlobalIndex++;
        }
        
        first = false;
    }
};

#endif /*_ONEDIMCARDIACMECHANICSASSEMBLER_HPP_*/
