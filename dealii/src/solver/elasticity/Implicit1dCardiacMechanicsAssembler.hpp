#ifndef _IMPLICIT1DCARDIACMECHANICSASSEMBLER_HPP_
#define _IMPLICIT1DCARDIACMECHANICSASSEMBLER_HPP_

#include "AbstractElasticityAssembler.hpp"
#include "PoleZero3dIn1dLaw.hpp"
#include "Abstract1dCardiacMechanicsAssembler.hpp"
#include "NhsSystemWithImplicitSolver.hpp"

///\todo: 
// more doxy
// test of analytic jacobian


/**
 *  An Implicit 1d cardiac mechanics assembler.
 *  
 *  See Abstract1dCardiacMechanicsAssembler
 * 
 *  This class takes in the current [Ca]_i at each quad point and solves the elasticity
 *  equations implicit, which involves solving the NHS cell models implicitly.
 */
class Implicit1dCardiacMechanicsAssembler : public Abstract1dCardiacMechanicsAssembler
{
private:
    std::vector<NhsSystemWithImplicitSolver> mCellMechSystems; 
    std::vector<double> mLambdaLastTimeStep;

    double mCurrentTime;
    double mNextTime;
    double mDt;


public:
    Implicit1dCardiacMechanicsAssembler(Triangulation<1>* pMesh)
        :  Abstract1dCardiacMechanicsAssembler(pMesh)
    {
        mCellMechSystems.resize(mTotalQuadPoints);
        mLambdaLastTimeStep.resize(mTotalQuadPoints, 1.0);
    }    

    void SetForcingQuantity(std::vector<double>& caI)
    {
        assert(caI.size() == mTotalQuadPoints);
        for(unsigned i=0; i<caI.size(); i++)
        {
            mCellMechSystems[i].SetIntracellularCalciumConcentration(caI[i]);
        } 
    }

    void Solve(double currentTime, double nextTime, double timestep)
    {
        assert(currentTime < nextTime);
        mCurrentTime = currentTime;
        mNextTime = nextTime;
        mDt = timestep;
        
        Abstract1dCardiacMechanicsAssembler::Solve(currentTime,nextTime,timestep);
    
        for(unsigned i=0; i<mCellMechSystems.size(); i++)
        {
             mCellMechSystems[i].UpdateStateVariables();
             mLambdaLastTimeStep[i] = mCellMechSystems[i].GetLambda();
        }
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
            
            double F = 1 + grad_u[0][0];
            double C = F*F;
            double T;
            double dTdE;            

            // get the active tension at this quad point
            //double active_tension = mActiveTension[mCurrentQuadPointGlobalIndex];
            
            double lam = F;
            double dlam_dt = (lam-mLambdaLastTimeStep[mCurrentQuadPointGlobalIndex])/(mNextTime-mCurrentTime);

            // get active tension for (lam+h,dlamdt)
            double h1 = std::max(1e-8, lam/100);
            mCellMechSystems[mCurrentQuadPointGlobalIndex].SetLambdaAndDerivative(lam+h1, dlam_dt);
            mCellMechSystems[mCurrentQuadPointGlobalIndex].SolveDoNotUpdate(mCurrentTime,mNextTime,mDt);
            double active_tension_at_lam_plus_h = mCellMechSystems[mCurrentQuadPointGlobalIndex].GetActiveTensionAtNextTime();        

            // get active tension for (lam,dlamdt+h)
            double h2 = std::max(1e-8, dlam_dt/100);
            mCellMechSystems[mCurrentQuadPointGlobalIndex].SetLambdaAndDerivative(lam, dlam_dt+h2);
            mCellMechSystems[mCurrentQuadPointGlobalIndex].SolveDoNotUpdate(mCurrentTime,mNextTime,mDt);
            double active_tension_at_dlamdt_plus_h = mCellMechSystems[mCurrentQuadPointGlobalIndex].GetActiveTensionAtNextTime();        

            // get proper active tension
            mCellMechSystems[mCurrentQuadPointGlobalIndex].SetLambdaAndDerivative(lam, dlam_dt);
            mCellMechSystems[mCurrentQuadPointGlobalIndex].SolveDoNotUpdate(mCurrentTime,mNextTime,mDt);
            double active_tension = mCellMechSystems[mCurrentQuadPointGlobalIndex].GetActiveTensionAtNextTime();        

            double d_act_tension_dlam = (active_tension_at_lam_plus_h - active_tension)/h1;
            double d_act_tension_d_dlamdt = (active_tension_at_dlamdt_plus_h - active_tension)/h2;


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
 
                        // differentiate Ta wrt lam and dlamdt
                        elementMatrix(i,j)  +=    (
                                                       d_act_tension_dlam
                                                     +
                                                       d_act_tension_d_dlamdt/mDt
                                                  )
                                                * (fe_values.shape_grad(j,q_point)[0]/C)                 
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

#endif /*_IMPLICIT1DCARDIACMECHANICSASSEMBLER_HPP_*/
