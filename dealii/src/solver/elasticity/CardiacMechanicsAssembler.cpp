#ifndef CARDIACMECHANICSASSEMBLER_CPP_
#define CARDIACMECHANICSASSEMBLER_CPP_
/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "CardiacMechanicsAssembler.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "NashHunterPoleZeroLaw.hpp"

template<unsigned DIM>
CardiacMechanicsAssembler<DIM>::CardiacMechanicsAssembler(Triangulation<DIM>* pMesh,
                                                          std::string outputDirectory,
                                                          AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw)
    : FiniteElasticityAssembler<DIM>(pMesh, pMaterialLaw, Vector<double>(DIM), 1.0, outputDirectory, 2, 1),
      AbstractCardiacMechanicsAssembler<DIM>(pMesh)
{
    mAllocatedMaterialLawMemory = false;

    if(pMaterialLaw==NULL)
    {
        this->mMaterialLaws.resize(1);
        this->mMaterialLaws[0] = new NashHunterPoleZeroLaw<DIM>;
        this->mHeterogeneous = false; // as FiniteElasticityAssembler would have set this to be true as NULL passed in to construtor
        mAllocatedMaterialLawMemory = true;
    }

    for(unsigned i=0; i<DIM; i++)
    {
        for(unsigned j=0; j<DIM; j++)
        {
            mFibreSheetMat[i][j] = i==j ? 1.0 : 0.0;
            
        }
    }
    mTransFibreSheetMat = transpose(mFibreSheetMat);
    mScaleFactor = 1.0;
}

template<unsigned DIM>
CardiacMechanicsAssembler<DIM>::~CardiacMechanicsAssembler()
{
    if(mAllocatedMaterialLawMemory)
    {    
        delete this->mMaterialLaws[0];
    }
}


template<unsigned DIM>
void CardiacMechanicsAssembler<DIM>::SetForcingQuantity(std::vector<double>& activeTension)
{
    assert(activeTension.size() == this->mTotalQuadPoints);
    mActiveTension = activeTension;
}

template<unsigned DIM>
void CardiacMechanicsAssembler<DIM>::Solve(double currentTime, double nextTime, double timestep)
{
    // do nothing with the times (as explicit) and call Solve on the base class
    FiniteElasticityAssembler<DIM>::StaticSolve(false);
}

template<unsigned DIM>
void CardiacMechanicsAssembler<DIM>::SetScaling(double scaleFactor)
{
    assert(mScaleFactor > 0.0);
    mScaleFactor = scaleFactor;
    for(unsigned i=0; i<this->mMaterialLaws.size(); i++)
    {
        std::cout << "setting scaling on " << i << "\n";
        this->mMaterialLaws[i]->ScaleMaterialParameters(mScaleFactor);
    }
}


template<unsigned DIM>
void CardiacMechanicsAssembler<DIM>::SetFibreSheetMatrix(Tensor<2,DIM> fibreSheetMat)
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


/*************************************
 *  AssembleOnElement
 * 
 *  Differs FiniteElasticityAssembler::AssembleOnElement in a few ways:  
 *   1. The active tension at a quad point is used, and the stretch at the quad point is set.
 *   
 *   2. The extra term in the stress (see below) arising from the active tension is incorporated 
 *      into the formulation. It is added to the stress (the variable 'T'), and dTdE is amended 
 *      as well. The fibre direction is taken into account using a rotation matrix.
 * 
 *   3. Since the body force and Neumann tractions will be zero, the corresponding loops
 *      have been removed. 
 * 
 *  The active tension term in the stress is T_a/C_fibre[0][0], when C_fibre = P^T C P,
 *  where P is a rotation matrix rotating the axes onto the fibre-sheet axes.
 *      
 *************************************/
template<unsigned DIM>
void CardiacMechanicsAssembler<DIM>::AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                                                       Vector<double>&       elementRhs,
                                                       FullMatrix<double>&   elementMatrix,
                                                       bool                  assembleResidual,
                                                       bool                  assembleJacobian
                                                      )
{
    // if mCurrentQuadPointGlobalIndex is greater than the total num of quad points something
    // very bad has happened. 
    assert(this->mCurrentQuadPointGlobalIndex <= this->mTotalQuadPoints);
    
    if(this->mCurrentQuadPointGlobalIndex==this->mTotalQuadPoints)
    {
        // if we are not back to the first cell something bad has happened
        assert( elementIter == this->mDofHandler.begin_active() );
        
        this->mCurrentQuadPointGlobalIndex = 0;
    }


    static QGauss<DIM>   quadrature_formula(this->mNumQuadPointsInEachDimension);
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
    

//// for a varying fibre-direction
//    assert(DIM==2);
//    double   theta = 0.785398163/5 * elementIter->vertex(0)[0]; //0->pi/20
//    mFibreSheetMat[0][0] =  cos(theta);
//    mFibreSheetMat[0][1] =  sin(theta);
//    mFibreSheetMat[1][0] = -sin(theta);
//    mFibreSheetMat[1][1] =  cos(theta);
//    mTransFibreSheetMat = transpose(mFibreSheetMat);
    
    
    for (unsigned q_point=0; q_point<n_q_points; q_point++)
    {
        const std::vector< Tensor<1,DIM> >& grad_u_p = local_solution_gradients[q_point];
        
        double p = local_solution_values[q_point](this->PRESSURE_COMPONENT_INDEX);

        static Tensor<2,DIM> F;
        static Tensor<2,DIM> C;
        static Tensor<2,DIM> inv_C;
        static Tensor<2,DIM> inv_F;
        static Tensor<2,DIM> T;
        
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

        double detF = determinant(F);
        
        /*************************************
         *  The cardiac-specific code
         ************************************/


        static Tensor<2,DIM> C_fibre;          // C when transformed to fibre-sheet axes
        static Tensor<2,DIM> inv_C_fibre;      // C^{-1} transformed to fibre-sheet axes
        static SymmetricTensor<2,DIM> T_fibre; // T when transformed to fibre-sheet axes
        
        // transform C and invC
        C_fibre = mTransFibreSheetMat * C * mFibreSheetMat;
        inv_C_fibre = mTransFibreSheetMat * inv_C * mFibreSheetMat;

        // store the stretch in the fibre direction
        this->mLambda[this->mCurrentQuadPointGlobalIndex] = sqrt(C_fibre[0][0]);

        // get the active tension at this quad point
        double active_tension = mActiveTension[this->mCurrentQuadPointGlobalIndex]/mScaleFactor;


        //mDTdE_fibre.Zero();

        // compute the transformed tension. The material law should law be a cardiac-
        // specific law which assumes the x-axes in the fibre, the z-axes the sheet normal
        p_material_law->ComputeStressAndStressDerivative(C_fibre,inv_C_fibre,p,T_fibre,mDTdE_fibre,assembleJacobian);

        // amend the stress and dTdE using the active tension
        T_fibre[0][0] += active_tension/C_fibre[0][0];
        mDTdE_fibre(0,0,0,0) -= 2*active_tension/(C_fibre[0][0]*C_fibre[0][0]);  
        
        // transform T back to real coordinates
        // Note we explicitly do the multiplication as can't multiply
        // deal.II SymmetricTensor with a Tensor

///\todo: make efficient
        for(unsigned M=0; M<DIM; M++) 
        {
            for(unsigned N=0; N<DIM; N++)
            {
                T[M][N] = 0;        
                for(unsigned al=0; al<DIM; al++) 
                {
                    for(unsigned be=0; be<DIM; be++)
                    {
                        T[M][N] +=                T_fibre [al][be]
                                    *      mFibreSheetMat [M][al]
                                    * mTransFibreSheetMat [be][N];
                    }
                }
            }
        }            


        
//        // transform dTdE back to real coords (ie dT_{albe}dE_{gam de} to dT_{MN}dE_{PQ})
//        for(unsigned M=0; M<DIM; M++) 
//        {
//            for(unsigned N=0; N<DIM; N++)
//            {
//                for(unsigned P=0; P<DIM; P++) 
//                {
//                    for(unsigned Q=0; Q<DIM; Q++)
//                    {
//                        this->dTdE(M,N,P,Q) = 0;
//                        for(unsigned al=0; al<DIM; al++) 
//                        {
//                            for(unsigned be=0; be<DIM; be++)
//                            {
//                                for(unsigned gam=0; gam<DIM; gam++)
//                                {
//                                    for(unsigned de=0; de<DIM; de++)
//                                    {
//                                        this->dTdE(M,N,P,Q) +=            mDTdE_fibre (al,be,gam,de)
//                                                                *      mFibreSheetMat [M][al]
//                                                                * mTransFibreSheetMat [be][N]
//                                                                * mTransFibreSheetMat [gam][P]
//                                                                *      mFibreSheetMat [Q][de];
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }            

        static FourthOrderTensor<DIM> temp1;
        static FourthOrderTensor<DIM> temp2;
        static FourthOrderTensor<DIM> temp3;

        temp1.SetAsProduct(mDTdE_fibre, mFibreSheetMat, 0);
        temp2.SetAsProduct(temp1,       mFibreSheetMat, 1);
        temp3.SetAsProduct(temp2,       mFibreSheetMat, 2);
        
        this->dTdE.SetAsProduct(temp3, mFibreSheetMat, 3);

        
        /********************************
         * end of cardiac specific code
         ********************************/
///\todo: refactor somehow 

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
                                                              * this->dTdE(M,N,P,Q)
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

//// implementation of old (wrong) equation, where T = .. + T_a invF_{0M} invF_{0N}
//// Note T is now directly altered, so no need to add anything new to elementMatrix
//                                ///////////////////////////////////////////////////////////
//                                // The extra part of the element stiffness matrix 
//                                // arising from the active tension part of the stress
//                                ///////////////////////////////////////////////////////////
//                                elementMatrix(i,j) +=  -active_tension
//                                                      * identity[component_i][0]
//                                                      * detF
//                                                      * inv_F[N][component_j]
//                                                      * inv_F[M][0]
//                                                      * fe_values.shape_grad(j,q_point)[M]
//                                                      * fe_values.shape_grad(i,q_point)[N]
//                                                      * fe_values.JxW(q_point);

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

//// implementation of old (wrong) equation, where T = .. + T_a invF_{0M} invF_{0N}
//// Note T is now directly altered, so no need to add anything new to elementRhs
//                        ///////////////////////////////////////////////////////////
//                        // The extra part of the element stiffness matrix 
//                        // arising from the active tension part of the stress
//                        ///////////////////////////////////////////////////////////
//                        elementRhs(i) +=   active_tension
//                                         * identity[component_i][0]
//                                         * inv_F[N][0]
//                                         * fe_values.shape_grad(i,q_point)[N]
//                                         * fe_values.JxW(q_point);
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

        this->mCurrentQuadPointGlobalIndex++;
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

#endif /*CARDIACMECHANICSASSEMBLER_CPP_*/
