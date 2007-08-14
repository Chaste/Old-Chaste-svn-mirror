#ifndef _1DCARDIACMECHANICSASSEMBLER_HPP_
#define _1DCARDIACMECHANICSASSEMBLER_HPP_

#include "FiniteElasticityAssembler.hpp" // just for tolerances...
#include "AbstractDealiiAssembler.hpp"

/* Todo:
 * 
 * add T(C) function
 * compute dTdE
 * 
 * change algorithm...
 */

class PoleZero3dIn1dLaw
{
    static const double TOL = 1e-8;
    
    std::vector<double> mPStore;
    std::vector<double> mE22Store;
    std::vector<double> mE33Store;
    bool mUseStore;
    
    double mLastE22;
    double mLastE33;
    
    std::vector<double> k;
    std::vector<double> a;
    std::vector<double> b;

    double dWdE(double Eii, unsigned index)
    {
        assert(index>=0 && index<=2);
        if(Eii >= a[index])
        {
            std::cout << index << " " << Eii << "\n";
            assert(0);
        }
        
        return   k[index] 
               * Eii 
               * (2+b[index]*Eii/(a[index]-Eii))
               / pow(a[index]-Eii,b[index]);
    }

    double d2WdE2(double Eii, unsigned index)
    {
        assert(index>=0 && index<=3);
        assert(Eii < a[index]);

        return    k[index] 
                * pow(a[index]-Eii, -b[index]-2) 
                * (   2*(a[index]-Eii)*(a[index]-Eii)
                    + 4*b[index]*Eii*(a[index]-Eii)
                    + b[index]*(b[index]+1)*Eii*Eii  );
    }

    
    void CalcResidual(Tensor<1,3>& resid, const double E, double p, double E22, double E33)
    {
        resid[0] = (2*E+1)*(2*E22+1)*(2*E33+1) - 1;
        resid[1] = dWdE(E22,1) - p/(2*E22+1);
        resid[2] = dWdE(E33,2) - p/(2*E33+1);
    }
    
    double NormResidual(Tensor<1,3>& resid)
    {
        return sqrt( resid[0]*resid[0] + resid[1]*resid[1] + resid[2]*resid[2] );
    }


    void CalcJacobian(Tensor<2,3>& jac, const double E, double p, double E22, double E33)
    {
        jac[0][0] = 0;
        jac[0][1] = (2*E+1) * 2 * (2*E33+1);
        jac[0][2] = (2*E+1) * (2*E22+1) * 2;

        jac[1][0] = - 1/(2*E22+1);
        jac[1][1] = d2WdE2(E22,1) + 2*p/( (2*E22+1)*(2*E22+1) );
        jac[1][2] = 0.0;

        jac[2][0] = - 1/(2*E33+1);
        jac[2][1] = 0.0;
        jac[2][2] = d2WdE2(E33,2) + 2*p/( (2*E33+1)*(2*E33+1) );
    }

        
    double SolveSystem(const double E)
    {
        Tensor<1,3> resid;
        Tensor<2,3> jac;
        
        double p = 0.0;
        double C22 = 1/sqrt(2*E+1);
        double E22 = 0.5*(C22-1);   
        double E33 = E22;
        
        if(mUseStore && E>=-0.3)
        {
            unsigned index = (unsigned)floor( (E+0.3)*mPStore.size() );
            assert(index < mPStore.size());
            p = mPStore[index];
            E22 = mE22Store[index];
            E33 = mE33Store[index];
        }
        
        CalcResidual(resid, E, p, E22, E33);
        double norm = NormResidual(resid);
        //std::cout << "E,||r||,p,E22,E33 = " << E << " " << norm << " " << p << " " << E22 << " " << E33 << "\n";
        unsigned counter = 0;

//        Tensor<2,3> num_jac;
//        double h=0.001;
//
//        Tensor<1,3> resid2;
//        CalcResidual(resid2, E, p+h, E22, E33);
//        for(unsigned i=0; i<3; i++)
//        {
//            num_jac[i][0] = (resid2[i]-resid[i])/h;
//        }
//                 
//
//        CalcResidual(resid2, E, p, E22+h, E33);
//        for(unsigned i=0; i<3; i++)
//        {
//            num_jac[i][1] = (resid2[i]-resid[i])/h;
//        }
//
//        CalcResidual(resid2, E, p, E22, E33+h);
//        for(unsigned i=0; i<3; i++)
//        {
//            num_jac[i][2] = (resid2[i]-resid[i])/h;
//        }
//
//
//            CalcJacobian(jac, E, p, E22, E33);
//
//            for(unsigned i=0; i<3; i++)
//            {
//                for(unsigned j=0; j<3; j++)
//                {
//                    std::cout << i << " " << j << " " << num_jac[i][j]-jac[i][j] << ", " << num_jac[i][j] << " " << jac[i][j] << "\n";
//                }
//            }        
//
//
//
//assert(0);


        while ((norm > TOL) && (counter<20))
        {
            CalcJacobian(jac, E, p, E22, E33);
            Tensor<2,3> invJac = invert(jac);
            
            std::vector<double> update(3,0.0);

            for(unsigned i=0; i<3; i++)
            {
                for(unsigned j=0; j<3; j++)
                {
                   update[i] += invJac[i][j]*resid[j];
                }
            }        

            double damping = ChooseBestUpdate(E, update, p, E22, E33);

            p  -= damping*update[0];
            E22 -= damping*update[1];
            E33 -= damping*update[2];
            
            CalcResidual(resid, E, p, E22, E33);
            norm = NormResidual(resid);
            //std::cout << "E,||r||,p,E22,E33 = " << E << " " <<  norm << " " << p << " " << E22 << " " << E33 << "\n\n\n";

            counter++;
        }

        assert(counter<20);

        mLastE22 = E22;
        mLastE33 = E33;

        return p;
    }

    double ChooseBestUpdate(double E, std::vector<double>& update, double p, double E22, double E33)
    {
        std::vector<double> try_vars = update;
        double best_damping = 0.0;

        Tensor<1,3> resid;
        CalcResidual(resid, E, p, E22, E33);
        double best_norm = NormResidual(resid);

        for(unsigned i=1; i<=10; i++)
        {
            double damping = i/10.0;

            try_vars[0] = p   - damping*update[0];
            try_vars[1] = E22 - damping*update[1];
            try_vars[2] = E33 - damping*update[2];
            
            CalcResidual(resid, E, try_vars[0], try_vars[1], try_vars[2]);
            double norm = NormResidual(resid);
            
            if(norm < best_norm)
            {
                best_norm = norm;
                best_damping = damping;
            }
        }
        
        if(best_damping == 0)
        {
            assert(0);
        }
        
        //std::cout << "Best s = " << best_damping << "\n";
        
        return best_damping;
    }

            
            
            
            


public : 
    PoleZero3dIn1dLaw()
    {
        k.resize(3);
        a.resize(3);
        b.resize(3);
        
        k[0] = 2;
        k[1] = 2;
        k[2] = 2;
        a[0] = 0.476;
        a[1] = 0.619;
        a[2] = 0.943;
        b[0] = 1.5;
        b[1] = 1.5;
        b[2] = 0.442;
        
        mUseStore = false;
    }
    
    void SetUpStores()
    {
        unsigned num = 1000; 
        for(unsigned i=0; i<=num; i++)
        {
            double E =  -0.3 + i*0.3/num;
            double p = SolveSystem(E);
            
            mPStore.push_back(p);
            mE22Store.push_back(mLastE22);
            mE33Store.push_back(mLastE33);
        }
            
        mUseStore = true;       
    }
        

    double GetT(const double E)
    {
        if(E>=0)
        {
            return dWdE(E,0); 
        }
        else
        {
            double p = SolveSystem(E);
            return -p/(2*E+1);
        }
    }
    
    double GetDTdE(const double E)
    {
        double h = 0.001;
        return (GetT(E+h)-GetT(E))/h;
    }
};




class OneDimCardiacMechanicsAssembler : public AbstractDealiiAssembler<1>
{
    FE_Q<1> mFe;
    double mDensity;
    unsigned mNumNewtonIterations;
    std::vector<Vector<double> > mDeformedPosition;
    std::vector<Vector<double> > mUndeformedPosition;


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
        : AbstractDealiiAssembler<1>(pMesh),
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


    void Solve()
    {
        //WriteOutput(0);
    
        // if nothing has been solved for yet, form an initial guess which is
        // the zero deformation solution (other the current solution is the best
        // initial guees)
        //if(!mADeformedHasBeenSolved)
        //{
        //    FormInitialGuess();
        //}
    
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

    void SetActiveTension(std::vector<double> activeTension)
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
        
        //AbstractIncompressibleMaterialLaw<1>* p_material_law = GetMaterialLawForElement(elementIter);
        

        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            //const std::vector< Tensor<1,1> >& grad_u_p = local_solution_gradients[q_point];
            const std::vector< Tensor<1,1> >& grad_u = local_solution_gradients[q_point];
            
//            double p = local_solution_values[q_point](this->PRESSURE_COMPONENT_INDEX);            
//            static Tensor<2,DIM> F;
//            static Tensor<2,DIM> C;
//            static Tensor<2,DIM> inv_C;
//            static Tensor<2,DIM> inv_F;
//            static Tensor<2,DIM> T;
//            
//            for (unsigned i=0; i<DIM; i++)
//            {
//                for (unsigned j=0; j<DIM; j++)
//                {
//                    F[i][j] = identity[i][j] + grad_u_p[i][j];
//                }
//            }
//            
//            C = transpose(F) * F;
//            
//            inv_C = invert(C);
//            inv_F = invert(F);
//    
//            double detF = determinant(F);

            double F = 1 + grad_u[0][0];;
            double C = F*F;
            //double inv_C = 1.0/C;
            //double inv_F = 1.0/F;
            double T;
            double dTdE;            
            //double detF = F;
            

            /*************************************
             *  The cardiac-specific code
             ************************************/
    
            // get the active tension at this quad point
            double active_tension = mActiveTension[mCurrentQuadPointGlobalIndex];
    
//            static Tensor<2,DIM> C_fibre;          // C when transformed to fibre-sheet axes
//            static Tensor<2,DIM> inv_C_fibre;      // C^{-1} transformed to fibre-sheet axes
//            static SymmetricTensor<2,DIM> T_fibre; // T when transformed to fibre-sheet axes
//            
//            // transform C and invC
//            C_fibre = mTransFibreSheetMat * C * mFibreSheetMat;
//            inv_C_fibre = mTransFibreSheetMat * inv_C * mFibreSheetMat;
//    
//            // store the stretch in the fibre direction
//            mLambda[mCurrentQuadPointGlobalIndex] = sqrt(C_fibre[0][0]);

            mLambda[mCurrentQuadPointGlobalIndex] = F; //=sqrt(C);

//            p_material_law->ComputeStressAndStressDerivative(C_fibre,inv_C_fibre,p,T_fibre,mDTdE_fibre,assembleJacobian);
//            // amend the stress and dTdE using the active tension
//            T_fibre[0][0] += active_tension/C_fibre[0][0];
//            mDTdE_fibre(0,0,0,0) -= 2*active_tension/(C_fibre[0][0]*C_fibre[0][0]);  
//
//            for(unsigned M=0; M<DIM; M++) 
//            {
//                for(unsigned N=0; N<DIM; N++)
//                {
//                    T[M][N] = 0;        
//                    for(unsigned al=0; al<DIM; al++) 
//                    {
//                        for(unsigned be=0; be<DIM; be++)
//                        {
//                            T[M][N] +=                T_fibre [al][be]
//                                        *      mFibreSheetMat [M][al]
//                                        * mTransFibreSheetMat [be][N];
//                        }
//                    }
//                }
//            }            
//            static FourthOrderTensor<DIM> temp1;
//            static FourthOrderTensor<DIM> temp2;
//            static FourthOrderTensor<DIM> temp3;
//    
//            temp1.SetAsProduct(mDTdE_fibre, mFibreSheetMat, 0);
//            temp2.SetAsProduct(temp1,       mFibreSheetMat, 1);
//            temp3.SetAsProduct(temp2,       mFibreSheetMat, 2);
//            
//            this->dTdE.SetAsProduct(temp3, mFibreSheetMat, 3);


// Compute T(C), dTdE(C)...
//assert(0); 
double E = 0.5*(C-1);
T = mLaw.GetT(E) + active_tension/C; //      = 0.5*(C-1) + active_tension;
dTdE = mLaw.GetDTdE(E) - active_tension/(C*C);
    
            
            /********************************
             * end of cardiac specific code
             ********************************/
    
            for (unsigned i=0; i<dofs_per_element; i++)
            {
                //const unsigned component_i = this->mFeSystem.system_to_component_index(i).first;
                
                if (assembleJacobian)
                {
                    for (unsigned j=0; j<dofs_per_element; j++)
                    {
//                        const unsigned component_j = this->mFeSystem.system_to_component_index(j).first;
                        
//                        if ((component_i<this->PRESSURE_COMPONENT_INDEX) &&(component_j<this->PRESSURE_COMPONENT_INDEX) )
//                        {
//                            for (unsigned N=0; N<DIM; N++)
//                            {
//                                for (unsigned M=0; M<DIM; M++)
//                                {
                                    elementMatrix(i,j) +=   T                                  //  T[M][N]
                                                          * fe_values.shape_grad(j,q_point)[0]
                                                          * fe_values.shape_grad(i,q_point)[0]
                                                          //* identity[component_i][component_j]
                                                          * fe_values.JxW(q_point);
                                                            
//                                    for (unsigned P=0; P<DIM; P++)
//                                    {
//                                        for (unsigned Q=0; Q<DIM; Q++)
//                                        {
                                    elementMatrix(i,j) +=   0.5
                                                          * dTdE //this->dTdE(M,N,P,Q)
                                                          * (
                                                              fe_values.shape_grad(j,q_point)[0]
                                                            * F                  //F[component_j][P]
                                                            +
                                                              fe_values.shape_grad(j,q_point)[0]
                                                            * F                  //F[component_j][Q]
                                                            )
                                                            *F                   //F[component_i][M]
                                                            * fe_values.shape_grad(i,q_point)[0]
                                                            * fe_values.JxW(q_point);
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                        else if ((component_i<this->PRESSURE_COMPONENT_INDEX) &&(component_j==this->PRESSURE_COMPONENT_INDEX) )
//                        {
//                            for (unsigned M=0; M<DIM; M++)
//                            {
//                                for (unsigned N=0; N<DIM; N++)
//                                {
//                                    elementMatrix(i,j) +=  - F[component_i][M]
//                                                           * inv_C[M][N]
//                                                           * fe_values.shape_grad(i,q_point)[N]
//                                                           * fe_values.shape_value(j,q_point)
//                                                           * fe_values.JxW(q_point);
//                                }
//                            }
//                        }
//                        else if ((component_i==this->PRESSURE_COMPONENT_INDEX) &&(component_j<this->PRESSURE_COMPONENT_INDEX) )
//                        {
//                            for (unsigned M=0; M<DIM; M++)
//                            {
//                                elementMatrix(i,j) +=    fe_values.shape_value(i,q_point)
//                                                       * detF
//                                                       * inv_F[M][component_j]
//                                                       * fe_values.shape_grad(j,q_point)[M]
//                                                       * fe_values.JxW(q_point);
//                            }
//                        }
//                        //else
//                        //{
//                            // do nothing, ie elementMatrix(i,j)  +=  0 * fe_values.JxW(q_point);;
//                        //}
                    }
                }
                
                if (assembleResidual)
                {
//                    if (component_i<this->PRESSURE_COMPONENT_INDEX)
//                    {
                        // body force so nothing here
                                         
//                        for (unsigned N=0; N<DIM; N++)
//                        {
//                            for (unsigned M=0; M<DIM; M++)
//                            {
                                elementRhs(i) +=   T                    //T[M][N]
                                                 * F                    //F[component_i][M]
                                                 * fe_values.shape_grad(i,q_point)[0]
                                                 * fe_values.JxW(q_point);
//                            }
//                        }
//                    }
//                    else
//                    {
//                        elementRhs(i) +=   fe_values.shape_value(i,q_point)
//                                         * (detF - 1)
//                                         * fe_values.JxW(q_point);
//                    }
                }
            }
    
            mCurrentQuadPointGlobalIndex++;
        }
        
        first = false;
    }


public: 
    std::vector<Vector<double> >& rGetDeformedPosition()
    {
        mDeformedPosition.clear();
        mDeformedPosition.resize(1);
        mDeformedPosition[0].reinit(this->mpMesh->n_vertices());
    
        DofVertexIterator<1> vertex_iter(this->mpMesh, &this->mDofHandler);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<1> old_posn = vertex_iter.GetVertex();
            
            mDeformedPosition[0](vertex_index) =   old_posn(0)
                                                 + mCurrentSolution(vertex_iter.GetDof(0));
            vertex_iter.Next();
        }
    
        return mDeformedPosition;
    }

    std::vector<Vector<double> >& rGetUndeformedPosition()
    {
        mUndeformedPosition.clear();
        mUndeformedPosition.resize(1);
        mUndeformedPosition[0].reinit(this->mpMesh->n_vertices());
    
        TriangulationVertexIterator<1> vertex_iter(this->mpMesh);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<1> old_posn = vertex_iter.GetVertex();
            mUndeformedPosition[0](vertex_index) = old_posn(0);
            vertex_iter.Next();
        }
    
        return mUndeformedPosition;
    }
};

#endif /*_1DCARDIACMECHANICSASSEMBLER_HPP_*/
