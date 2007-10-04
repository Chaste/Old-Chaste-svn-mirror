#ifndef _EXPLICIT1DCARDIACMECHANICSASSEMBLER_HPP_
#define _EXPLICIT1DCARDIACMECHANICSASSEMBLER_HPP_

#include "Abstract1dCardiacMechanicsAssembler.hpp"

/**
 *  An Explicit 1d cardiac mechanics assembler.
 *  
 *  See Abstract1dCardiacMechanicsAssembler
 * 
 *  This class takes in the current active tension at each quad point and uses that
 *  to compute the stress there.
 */
class Explicit1dCardiacMechanicsAssembler : public Abstract1dCardiacMechanicsAssembler
{
protected :
    std::vector<double> mActiveTension;

public:
    Explicit1dCardiacMechanicsAssembler(Triangulation<1>* pMesh)
        : Abstract1dCardiacMechanicsAssembler(pMesh)
    {
        assert(mTotalQuadPoints > 0);
        mActiveTension.resize(mTotalQuadPoints, 1.0);
    }    

    void SetForcingQuantity(std::vector<double>& activeTension)
    {
        assert(activeTension.size() == mTotalQuadPoints);
        mActiveTension = activeTension;
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

            // Compute T(C), dTdE(C), using the active tension
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

#endif /*_EXPLICIT1DCARDIACMECHANICSASSEMBLER_HPP_*/
