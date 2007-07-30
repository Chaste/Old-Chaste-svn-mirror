#ifndef CARDIACMECHASSEMBLER_CPP_
#define CARDIACMECHASSEMBLER_CPP_

#include "CardiacMechAssembler.hpp"
#include "MooneyRivlinMaterialLaw.hpp"

template<unsigned DIM>
CardiacMechAssembler<DIM>::CardiacMechAssembler(Triangulation<DIM>* pMesh,
                                                          std::string outputDirectory,
                                                          AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw)
    : FiniteElasticityAssembler<DIM>(pMesh, pMaterialLaw, Vector<double>(DIM), 1.0, outputDirectory, 2, 1)
{
    mAllocatedMaterialLawMemory = false;

    if(pMaterialLaw==NULL)
    {
        this->mMaterialLaws.resize(1);
        this->mMaterialLaws[0] = new MooneyRivlinMaterialLaw<DIM>(0.02);
        this->mHeterogeneous = false; // as FiniteElasticityAssembler would have set this to be true as NULL passed in to construtor
        mAllocatedMaterialLawMemory = true;
    }
    
    // set up quad point info
    QGauss<DIM>   quadrature_formula(mNumQuadPointsInEachDimension);
    mTotalQuadPoints = quadrature_formula.n_quadrature_points *this->mpMesh->n_active_cells();
    mCurrentQuadPointGlobalIndex = 0;
    
    mLambda.resize(mTotalQuadPoints, 1.0);
}

template<unsigned DIM>
CardiacMechAssembler<DIM>::~CardiacMechAssembler()
{
    if(mAllocatedMaterialLawMemory)
    {    
        delete this->mMaterialLaws[0];
    }
}



template<unsigned DIM>
unsigned CardiacMechAssembler<DIM>::GetTotalNumQuadPoints()
{
    return mTotalQuadPoints;
}


template<unsigned DIM>
void CardiacMechAssembler<DIM>::SetActiveTension(std::vector<double> activeTension)
{
    assert(activeTension.size() == mTotalQuadPoints);
    mActiveTension = activeTension;
}

template<unsigned DIM>
std::vector<double>& CardiacMechAssembler<DIM>::GetLambda()
{
    return mLambda;
}




/*************************************
 *  AssembleOnElement
 * 
 *  Differs FiniteElasticityAssembler::AssembleOnElement in a few ways:
 *   1. It sets up the active tension at the gauss point by interpolating linearly 
 *      using the nodal values of the active tension
 *  
 *   2. The extra term in the stress (see below) arising from the active tension is incorporated 
 *      into the formulation. It is NOT added to the stress (the variable 'T'), but
 *      the extra terms in elementRhs and elementMatrix are computed and added directly
 * 
 *   3. Since the body force and Neumann tractions will be zero, the corresponding loops
 *      have been removed. 
 *      
 *  The fibre direction is currently assumed to be the x-direction, so the active tension
 *  term in the stress is T_a invF_{M0} invF_{N0}
 * 
 *************************************/
template<unsigned DIM>
void CardiacMechAssembler<DIM>::AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                                                       Vector<double>&       elementRhs,
                                                       FullMatrix<double>&   elementMatrix,
                                                       bool                  assembleResidual,
                                                       bool                  assembleJacobian
                                                      )
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


    static QGauss<DIM>   quadrature_formula(mNumQuadPointsInEachDimension);
    //static QGauss<DIM-1> face_quadrature_formula(3);
    
    const unsigned n_q_points    = quadrature_formula.n_quadrature_points;
    //const unsigned n_face_q_points = face_quadrature_formula.n_quadrature_points;
    
    
    // would want this to be static too (slight speed up), but causes errors
    // in debug mode (upon destruction of the class, in 2d, or something)
    FEValues<DIM> fe_values(this->mFeSystem, quadrature_formula,
                            UpdateFlags(update_values    |
                                        update_gradients |
                                        update_q_points  |     // needed for interpolating u and u' on the quad point
                                        update_JxW_values));
                                        
    //    // would want this to be static too (slight speed up), but causes errors
    //    // in debug mode (upon destruction of the class, in 2d, or something)
    //    FEFaceValues<DIM> fe_face_values(this->mFeSystem, face_quadrature_formula,
    //                                     UpdateFlags(update_values         |
    //                                                 update_q_points       |
    //                                                 update_normal_vectors |
    //                                                 update_JxW_values));
                                                 
                                                 
    const unsigned dofs_per_element = this->mFeSystem.dofs_per_cell;
    
    static std::vector< Vector<double> >                  local_solution_values(n_q_points);
    static std::vector< std::vector< Tensor<1,DIM> > >    local_solution_gradients(n_q_points);
    
    static Tensor<2,DIM> identity;
    
    static bool first = true;
    
    
    if (first)
    {
        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            local_solution_values[q_point].reinit(DIM+1);
            local_solution_gradients[q_point].resize(DIM+1);
        }
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned j=0; j<DIM; j++)
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
    
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law = GetMaterialLawForElement(elementIter);
    
    for (unsigned q_point=0; q_point<n_q_points; q_point++)
    {
        const std::vector< Tensor<1,DIM> >& grad_u_p = local_solution_gradients[q_point];
        
        double p = local_solution_values[q_point](this->PRESSURE_COMPONENT_INDEX);
        
        static Tensor<2,DIM> F;
        static Tensor<2,DIM> C;
        static Tensor<2,DIM> inv_C;
        static Tensor<2,DIM> inv_F;
        static SymmetricTensor<2,DIM> T;
        
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                F[i][j] = identity[i][j] + grad_u_p[i][j];
            }
        }
        
        C = transpose(F) * F;
        inv_C = invert(C);
        inv_F = invert(F);
        
        ///////////////////////////////////////////////////////////
        // Get the active tension at gauss point, and store lambda
        ///////////////////////////////////////////////////////////
        double active_tension = mActiveTension[mCurrentQuadPointGlobalIndex];
        mLambda[mCurrentQuadPointGlobalIndex] = sqrt(C[0][0]);

        
        double detF = determinant(F);
        
        p_material_law->ComputeStressAndStressDerivative(C,inv_C,p,T,this->dTdE,assembleJacobian);

                
        for (unsigned i=0; i<dofs_per_element; i++)
        {
            const unsigned component_i = this->mFeSystem.system_to_component_index(i).first;
            
            if (assembleJacobian)
            {
                for (unsigned j=0; j<dofs_per_element; j++)
                {
                    const unsigned component_j = this->mFeSystem.system_to_component_index(j).first;
                    
                    if ((component_i<this->PRESSURE_COMPONENT_INDEX) &&(component_j<this->PRESSURE_COMPONENT_INDEX) )
                    {
                        for (unsigned N=0; N<DIM; N++)
                        {
                            for (unsigned M=0; M<DIM; M++)
                            {
                                elementMatrix(i,j) +=   T[M][N]
                                                      * fe_values.shape_grad(j,q_point)[M]
                                                      * fe_values.shape_grad(i,q_point)[N]
                                                      * identity[component_i][component_j]
                                                      * fe_values.JxW(q_point);
                                                        
                                for (unsigned P=0; P<DIM; P++)
                                {
                                    for (unsigned Q=0; Q<DIM; Q++)
                                    {
                                        elementMatrix(i,j) +=   0.5
                                                              * this->dTdE[M][N][P][Q]
                                                              * (
                                                                  fe_values.shape_grad(j,q_point)[Q]
                                                                  * F[component_j][P]
                                                                  +
                                                                  fe_values.shape_grad(j,q_point)[P]
                                                                  * F[component_j][Q]
                                                              )
                                                              * F[component_i][M]
                                                              * fe_values.shape_grad(i,q_point)[N]
                                                              * fe_values.JxW(q_point);
                                    }
                                }
                                
                                ///////////////////////////////////////////////////////////
                                // The extra part of the element stiffness matrix 
                                // arising from the active tension part of the stress
                                ///////////////////////////////////////////////////////////
                                elementMatrix(i,j) +=  -active_tension
                                                      * identity[component_i][0]
                                                      * detF
                                                      * inv_F[N][component_j]
                                                      * inv_F[M][0]
                                                      * fe_values.shape_grad(j,q_point)[M]
                                                      * fe_values.shape_grad(i,q_point)[N]
                                                      * fe_values.JxW(q_point);

                            }
                        }
                    }
                    else if ((component_i<this->PRESSURE_COMPONENT_INDEX) &&(component_j==this->PRESSURE_COMPONENT_INDEX) )
                    {
                        for (unsigned M=0; M<DIM; M++)
                        {
                            for (unsigned N=0; N<DIM; N++)
                            {
                                elementMatrix(i,j) +=  - F[component_i][M]
                                                       * inv_C[M][N]
                                                       * fe_values.shape_grad(i,q_point)[N]
                                                       * fe_values.shape_value(j,q_point)
                                                       * fe_values.JxW(q_point);
                            }
                        }
                    }
                    else if ((component_i==this->PRESSURE_COMPONENT_INDEX) &&(component_j<this->PRESSURE_COMPONENT_INDEX) )
                    {
                        for (unsigned M=0; M<DIM; M++)
                        {
                            elementMatrix(i,j) +=    fe_values.shape_value(i,q_point)
                                                   * detF
                                                   * inv_F[M][component_j]
                                                   * fe_values.shape_grad(j,q_point)[M]
                                                   * fe_values.JxW(q_point);
                        }
                    }
                    //else
                    //{
                        // do nothing, ie elementMatrix(i,j)  +=  0 * fe_values.JxW(q_point);;
                    //}
                }
            }
            
            if (assembleResidual)
            {
                if (component_i<this->PRESSURE_COMPONENT_INDEX)
                {
                    /* body force is zero so do not do this: */
                    //elementRhs(i) += - mDensity * this->mBodyForce(component_i)
                    //                 * fe_values.shape_value(i,q_point)
                    //                 * fe_values.JxW(q_point);
                                     
                    for (unsigned N=0; N<DIM; N++)
                    {
                        for (unsigned M=0; M<DIM; M++)
                        {
                            elementRhs(i) +=   T[M][N]
                                             * F[component_i][M]
                                             * fe_values.shape_grad(i,q_point)[N]
                                             * fe_values.JxW(q_point);
                        }

                        ///////////////////////////////////////////////////////////
                        // The extra part of the element stiffness matrix 
                        // arising from the active tension part of the stress
                        ///////////////////////////////////////////////////////////
                        elementRhs(i) +=   active_tension
                                         * identity[component_i][0]
                                         * inv_F[N][0]
                                         * fe_values.shape_grad(i,q_point)[N]
                                         * fe_values.JxW(q_point);
                    }
                }
                else
                {
                    elementRhs(i) +=   fe_values.shape_value(i,q_point)
                                     * (detF - 1)
                                     * fe_values.JxW(q_point);
                }
            }
        }

        mCurrentQuadPointGlobalIndex++;
    }
    
    /* zero applied stress, so do not do this */
    //    
    //    ////////////////////////////
    //    // loop over faces
    //    ////////////////////////////
    //    if (assembleResidual)
    //    {
    //        for (unsigned face_index=0; face_index<GeometryInfo<DIM>::faces_per_cell; face_index++)
    //        {
    //            if (elementIter->face(face_index)->boundary_indicator() == NEUMANN_BOUNDARY)
    //            {
    //                fe_face_values.reinit(elementIter, face_index);
    //                
    //                for (unsigned q_point=0; q_point<n_face_q_points; q_point++)
    //                {
    //                    Vector<double> neumann_traction(DIM); // zeros
    //                    
    //                    neumann_traction(1)=0.0;
    //                    
    //                    for (unsigned i=0; i<dofs_per_element; i++)
    //                    {
    //                        const unsigned component_i = this->mFeSystem.system_to_component_index(i).first;
    //                        
    //                        if (component_i < this->PRESSURE_COMPONENT_INDEX)
    //                        {
    //                            elementRhs(i) +=   neumann_traction(component_i)
    //                                             * fe_face_values.shape_value(i,q_point)
    //                                             * fe_face_values.JxW(q_point);
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    
    first = false;
}

#endif /*CARDIACMECHASSEMBLER_CPP_*/
