/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "ImplicitCardiacMechanicsAssembler.hpp"

template<unsigned DIM>
ImplicitCardiacMechanicsAssembler<DIM>::ImplicitCardiacMechanicsAssembler(
                                  QuadraticMesh<DIM>* pQuadMesh,
                                  std::string outputDirectory,
                                  std::vector<unsigned>& rFixedNodes,
                                  AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw)
    : NonlinearElasticityAssembler<DIM>(pQuadMesh,
                                        pMaterialLaw!=NULL ? pMaterialLaw : new NashHunterPoleZeroLaw<DIM>,
                                        zero_vector<double>(DIM),
                                        DOUBLE_UNSET,
                                        outputDirectory,
                                        rFixedNodes),
      mCurrentTime(DBL_MAX),
      mNextTime(DBL_MAX),
      mOdeTimestep(DBL_MAX)
{
    // compute total num quad points
    mTotalQuadPoints = pQuadMesh->GetNumElements()*this->mpQuadratureRule->GetNumQuadPoints();

    // initialise stores
    mLambda.resize(mTotalQuadPoints, 1.0);
    mLambdaLastTimeStep.resize(mTotalQuadPoints, 1.0);
    mCellMechSystems.resize(mTotalQuadPoints);

    // note that if pMaterialLaw is NULL a new NashHunter law was sent to the
    // NonlinElas constuctor (see above)
    mAllocatedMaterialLawMemory = (pMaterialLaw==NULL);
}

template<unsigned DIM>
ImplicitCardiacMechanicsAssembler<DIM>::~ImplicitCardiacMechanicsAssembler()
{
    if(mAllocatedMaterialLawMemory)
    {
        assert(this->mMaterialLaws.size()==1); // haven't implemented heterogeniety yet
        delete this->mMaterialLaws[0];
    }
}

template<unsigned DIM>
unsigned ImplicitCardiacMechanicsAssembler<DIM>::GetTotalNumQuadPoints()
{
    return mTotalQuadPoints;
}

template<unsigned DIM>
GaussianQuadratureRule<DIM>* ImplicitCardiacMechanicsAssembler<DIM>::GetQuadratureRule()
{
    return this->mpQuadratureRule;
}

template<unsigned DIM>
void ImplicitCardiacMechanicsAssembler<DIM>::SetIntracellularCalciumConcentrations(std::vector<double>& caI)
{
    assert(caI.size() == mTotalQuadPoints);
    for(unsigned i=0; i<caI.size(); i++)
    {
        mCellMechSystems[i].SetIntracellularCalciumConcentration(caI[i]);
    }
}

template<unsigned DIM>
std::vector<double>& ImplicitCardiacMechanicsAssembler<DIM>::rGetLambda()
{
    return mLambda;
}


template<unsigned DIM>
void ImplicitCardiacMechanicsAssembler<DIM>::Solve(double currentTime, double nextTime, double odeTimestep)
{
    // set the times, which are used in AssembleOnElement
    assert(currentTime < nextTime);
    mCurrentTime = currentTime;
    mNextTime = nextTime;
    mOdeTimestep = odeTimestep;

    // solve
    NonlinearElasticityAssembler<DIM>::Solve();

    // assemble residual again (to solve the cell models implicitly again
    // using the correct value of the deformation x (in case this wasn't the
    // last thing that was done
    this->AssembleSystem(true,false);

    // now update state variables, and set lambda at last timestep. Note
    // lambda was set in AssembleOnElement
    for(unsigned i=0; i<mCellMechSystems.size(); i++)
    {
         mCellMechSystems[i].UpdateStateVariables();
         mLambdaLastTimeStep[i] = mCellMechSystems[i].GetLambda();
    }
}



template<unsigned DIM>
void ImplicitCardiacMechanicsAssembler<DIM>::AssembleOnElement(Element<DIM, DIM>& rElement,
                       c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElem,
                       c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElemPrecond,
                       c_vector<double,STENCIL_SIZE>& rBElem,
                       bool assembleResidual,
                       bool assembleJacobian)
{
    // check these have been set
    assert(mCurrentTime != DBL_MAX);
    assert(mNextTime != DBL_MAX);
    assert(mOdeTimestep != DBL_MAX);

    c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
    double jacobian_determinant;
    this->mpQuadMesh->GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

    if (assembleJacobian)
    {
        rAElem.clear();
        rAElemPrecond.clear();
    }

    if (assembleResidual)
    {
        rBElem.clear();
    }

    ///////////////////////////////////////////////
    // Get the current displacement at the nodes
    ///////////////////////////////////////////////
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> element_current_displacements;
    static c_vector<double,NUM_VERTICES_PER_ELEMENT> element_current_pressures;
    for(unsigned II=0; II<NUM_NODES_PER_ELEMENT; II++)
    {
        for(unsigned JJ=0; JJ<DIM; JJ++)
        {
            element_current_displacements(JJ,II) = this->mCurrentSolution[DIM*rElement.GetNodeGlobalIndex(II) + JJ];
        }
    }

    ///////////////////////////////////////////////
    // Get the current pressure at the vertices
    ///////////////////////////////////////////////
    for(unsigned II=0; II<NUM_VERTICES_PER_ELEMENT; II++)
    {
        element_current_pressures(II) = this->mCurrentSolution[DIM*this->mpQuadMesh->GetNumNodes() + rElement.GetNodeGlobalIndex(II)];
    }

    // allocate memory for the basis functions values and derivative values
    c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
    c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
    c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;

    // get the material law
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law;
    if(this->mMaterialLaws.size()==1)
    {
        // homogeneous
        p_material_law = this->mMaterialLaws[0];
    }
    else
    {
        // heterogeneous
        #define COVERAGE_IGNORE // not going to have tests in cts for everything
        p_material_law = this->mMaterialLaws[rElement.GetIndex()];
        #undef COVERAGE_IGNORE
    }

//// for a varying fibre-direction
//    assert(DIM==2);
//    double   theta = 0.785398163/5 * elementIter->vertex(0)[0]; //0->pi/20
//    this->mFibreSheetMat[0][0] =  cos(theta);
//    this->mFibreSheetMat[0][1] =  sin(theta);
//    this->mFibreSheetMat[1][0] = -sin(theta);
//    this->mFibreSheetMat[1][1] =  cos(theta);
//    this->mTransFibreSheetMat = transpose(this->mFibreSheetMat);


    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    //// loop over Gauss points
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    for (unsigned quadrature_index=0; quadrature_index < this->mpQuadratureRule->GetNumQuadPoints(); quadrature_index++)
    {
        unsigned current_quad_point_global_index =   rElement.GetIndex()*this->mpQuadratureRule->GetNumQuadPoints()
                                                   + quadrature_index;

        double wJ = jacobian_determinant * this->mpQuadratureRule->GetWeight(quadrature_index);

        const ChastePoint<DIM>& quadrature_point = this->mpQuadratureRule->rGetQuadPoint(quadrature_index);

        //////////////////////////////////////
        // set up basis function info
        //////////////////////////////////////
        LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
        QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
        QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_quad_phi);

        ////////////////////////////////////////////////////
        // (dont get the body force)
        ////////////////////////////////////////////////////

        //////////////////////////////////////
        // interpolate grad_u and p
        //////////////////////////////////////
        static c_matrix<double,DIM,DIM> grad_u; // grad_u = (du_i/dX_M)
        grad_u = zero_matrix<double>(DIM,DIM);  // must be on new line!!

        for(unsigned node_index=0;
            node_index<NUM_NODES_PER_ELEMENT;
            node_index++)
        {
            for (unsigned i=0; i<DIM; i++)
            {
                for(unsigned M=0; M<DIM; M++)
                {
                    grad_u(i,M) += grad_quad_phi(M,node_index)*element_current_displacements(i,node_index);
                }
            }
        }

        double pressure = 0;
        for(unsigned vertex_index=0;
            vertex_index<NUM_VERTICES_PER_ELEMENT;
            vertex_index++)
        {
            pressure += linear_phi(vertex_index)*element_current_pressures(vertex_index);
        }

        ///////////////////////////////////////////////
        // calculate C, inv(C) and T
        ///////////////////////////////////////////////
        static c_matrix<double,DIM,DIM> F;
        static c_matrix<double,DIM,DIM> C;
        static c_matrix<double,DIM,DIM> inv_C;
        static c_matrix<double,DIM,DIM> inv_F;
        static c_matrix<double,DIM,DIM> T;

        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned M=0; M<DIM; M++)
            {
                F(i,M) = (i==M?1:0) + grad_u(i,M);
            }
        }

        C = prod(trans(F),F);
        inv_C = Inverse(C);
        inv_F = Inverse(F);

        double detF = Determinant(F);

        /************************************
         *  The cardiac-specific code PART 1
         ************************************/
        //static c_matrix<double,DIM,DIM> C_fibre;          // C when transformed to fibre-sheet axes
        //static c_matrix<double,DIM,DIM> inv_C_fibre;      // C^{-1} transformed to fibre-sheet axes
        //static c_matrix<double,DIM,DIM> T_fibre; // T when transformed to fibre-sheet axes
        //
        //// transform C and invC
        //C_fibre = this->mTransFibreSheetMat * C * this->mFibreSheetMat;
        //inv_C_fibre = this->mTransFibreSheetMat * inv_C * this->mFibreSheetMat;

        // store the stretch in the fibre direction
        this->mLambda[current_quad_point_global_index] = sqrt(C(0,0));

        double lam = sqrt(C(0,0));
        double dlam_dt = (lam-mLambdaLastTimeStep[current_quad_point_global_index])/(mNextTime-mCurrentTime);

        NhsSystemWithImplicitSolver& system = mCellMechSystems[current_quad_point_global_index];

        // get proper active tension
        // see NOTE below
        system.SetLambdaAndDerivative(lam, dlam_dt);

        double active_tension=0.0;
        try
        {
            system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
            active_tension = system.GetActiveTensionAtNextTime();
        }
        catch (Exception& e)
        {
            #define COVERAGE_IGNORE
            if(assembleJacobian)
            {
                EXCEPTION("Failed in solving NHS systems when assembling Jacobian");
            }
            #undef COVERAGE_IGNORE
        }

        // compute the derivative of the active tension wrt lam and dlam_dt
        double d_act_tension_dlam = 0.0; //Set and used if assembleJacobian==true
        double d_act_tension_d_dlamdt = 0.0; //Set and used if assembleJacobian==true

        if(assembleJacobian)
        {
            // get active tension for (lam+h,dlamdt)
            double h1 = std::max(1e-6, lam/100);
            system.SetLambdaAndDerivative(lam+h1, dlam_dt);
            system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
            double active_tension_at_lam_plus_h = system.GetActiveTensionAtNextTime();

            // get active tension for (lam,dlamdt+h)
            double h2 = std::max(1e-6, dlam_dt/100);
            system.SetLambdaAndDerivative(lam, dlam_dt+h2);
            system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
            double active_tension_at_dlamdt_plus_h = system.GetActiveTensionAtNextTime();

            d_act_tension_dlam = (active_tension_at_lam_plus_h - active_tension)/h1;
            d_act_tension_d_dlamdt = (active_tension_at_dlamdt_plus_h - active_tension)/h2;
        }

        // NOTE - have to get the active tension again, this must be done last!!
        // As if this turns out to be the correct solution, the state vars will be updated!
        // TODO: sort out this inefficiency
        system.SetLambdaAndDerivative(lam, dlam_dt);
        system.SetActiveTensionInitialGuess(active_tension);

        try
        {
            system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
            assert( fabs(system.GetActiveTensionAtNextTime()-active_tension)<1e-8);
        }
        catch (Exception& e)
        {
            #define COVERAGE_IGNORE
            LOG(2, "WARNING in ImplicitCardiacMechanicsAssembler!\n");
            active_tension = 1e10;
            // should have done something above..
            #undef COVERAGE_IGNORE
        }

        //this->mDTdE_fibre.Zero();

        // compute the transformed tension. The material law should law be a cardiac-
        // specific law which assumes the x-axes in the fibre, the z-axes the sheet normal
        p_material_law->ComputeStressAndStressDerivative(C,inv_C,pressure,T,this->dTdE,assembleJacobian);

        // amend the stress and dTdE using the active tension
        T(0,0) += active_tension/C(0,0);
        this->dTdE(0,0,0,0) -= 2*active_tension/(C(0,0)*C(0,0));

//            //// could do this without the loop now
//            // transform T back to real coordinates
//            for(unsigned M=0; M<DIM; M++)
//            {
//                for(unsigned N=0; N<DIM; N++)
//                {
//                    T[M][N] = 0;
//                    for(unsigned al=0; al<DIM; al++)
//                    {
//                        for(unsigned be=0; be<DIM; be++)
//                        {
//                            T[M][N] +=                      T_fibre [al][be]
//                                        *      this->mFibreSheetMat [M] [al]
//                                        * this->mTransFibreSheetMat [be][N];
//                        }
//                    }
//                }
//            }
//            static FourthOrderTensor<DIM> temp1;
//            static FourthOrderTensor<DIM> temp2;
//            static FourthOrderTensor<DIM> temp3;
//
//            temp1.SetAsProduct(this->mDTdE_fibre, this->mFibreSheetMat, 0);
//            temp2.SetAsProduct(temp1,             this->mFibreSheetMat, 1);
//            temp3.SetAsProduct(temp2,             this->mFibreSheetMat, 2);
//
//            this->dTdE.SetAsProduct(temp3, this->mFibreSheetMat, 3);



        static FourthOrderTensor<DIM> dTdE_F;
        static FourthOrderTensor<DIM> dTdE_FF1;
        static FourthOrderTensor<DIM> dTdE_FF2;
  
        dTdE_F.SetAsProduct(this->dTdE, F, 0);  // B^{aNPQ}  = F^a_M * dTdE^{MNPQ}
        dTdE_FF1.SetAsProduct(dTdE_F, F, 3);    // B1^{aNPb} = F^a_M * F^b_Q * dTdE^{MNPQ} 
        dTdE_FF2.SetAsProduct(dTdE_F, F, 2);    // B2^{aNbQ} = F^a_M * F^b_P * dTdE^{MNPQ}


        /*************************************
         * end of cardiac specific code PART 1
         *************************************/

        /////////////////////////////////////////
        // residual vector
        /////////////////////////////////////////
        if (assembleResidual)
        {
            for(unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
            {
                unsigned spatial_dim = index%DIM;
                unsigned node_index = (index-spatial_dim)/DIM;

                assert(node_index < NUM_NODES_PER_ELEMENT);

                // no body force bit as body force = 0

                for (unsigned M=0; M<DIM; M++)
                {
                    for (unsigned N=0; N<DIM; N++)
                    {
                        rBElem(index) +=   T(M,N)
                                         * F(spatial_dim,M)
                                         * grad_quad_phi(N,node_index)
                                         * wJ;
                    }
                }
            }

            for(unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                rBElem( NUM_NODES_PER_ELEMENT*DIM + vertex_index ) +=   linear_phi(vertex_index)
                                                                      * (detF - 1)
                                                                      * wJ;
            }
        }

        /////////////////////////////////////////
        // Jacobian matrix
        /////////////////////////////////////////
        if(assembleJacobian)
        {
            for(unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
            {
                unsigned spatial_dim1 = index1%DIM;
                unsigned node_index1 = (index1-spatial_dim1)/DIM;


                for(unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    for (unsigned M=0; M<DIM; M++)
                    {
                        for (unsigned N=0; N<DIM; N++)
                        {
                            rAElem(index1,index2) +=   T(M,N)
                                                     * grad_quad_phi(N,node_index1)
                                                     * grad_quad_phi(M,node_index2)
                                                     * (spatial_dim1==spatial_dim2?1:0)
                                                     * wJ;

//                            for (unsigned P=0; P<DIM; P++)
//                            {
//                                for (unsigned Q=0; Q<DIM; Q++)
//                                {
//                                    rAElem(index1,index2)  +=   0.5
//                                                              * this->dTdE(M,N,P,Q)
//                                                              * (
//                                                                  grad_quad_phi(P,node_index2)
//                                                                * F(spatial_dim2,Q)
//                                                                   +
//                                                                  grad_quad_phi(Q,node_index2)
//                                                                * F(spatial_dim2,P)
//                                                                 )
//                                                              * F(spatial_dim1,M)
//                                                              * grad_quad_phi(N,node_index1)
//                                                              * wJ;
//                                }
//                            }
                        }
                    }
                    
                    for (unsigned N=0; N<DIM; N++)
                    {
                        for (unsigned P=0; P<DIM; P++)
                        {
                            rAElem(index1,index2)  +=   0.5
                                                      * dTdE_FF1(spatial_dim1,N,P,spatial_dim2)
                                                      * grad_quad_phi(P,node_index2)
                                                      * grad_quad_phi(N,node_index1)
                                                      * wJ;
                        }
                        
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                           rAElem(index1,index2)  +=   0.5
                                                     * dTdE_FF2(spatial_dim1,N,spatial_dim2,Q)
                                                     * grad_quad_phi(Q,node_index2)
                                                     * grad_quad_phi(N,node_index1)
                                                     * wJ;
                        }
                    }





                    /************************************
                     *  The cardiac-specific code PART 2
                     ************************************/
                    rAElem(index1,index2) +=  (
                                                   d_act_tension_dlam
                                                 +
                                                   d_act_tension_d_dlamdt/(mNextTime-mCurrentTime)
                                                )
                                                * (F(spatial_dim2,0)/lam)
                                                * grad_quad_phi(0,node_index2)
                                                * F(spatial_dim1,0)
                                                * grad_quad_phi(0,node_index1)
                                                * wJ;
                   /************************************
                    *  End cardiac-specific code PART 2
                    ************************************/
                }


                for(unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                    for (unsigned M=0; M<DIM; M++)
                    {
                        for (unsigned N=0; N<DIM; N++)
                        {
                            rAElem(index1,index2)  +=  - F(spatial_dim1,M)
                                                       * inv_C(M,N)
                                                       * grad_quad_phi(N,node_index1)
                                                       * linear_phi(vertex_index)
                                                       * wJ;
                        }
                    }
                }
            }

            for(unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                unsigned index1 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                for(unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    for (unsigned M=0; M<DIM; M++)
                    {
                        rAElem(index1,index2) +=   linear_phi(vertex_index)
                                                 * detF
                                                 * inv_F(M,spatial_dim2)
                                                 * grad_quad_phi(M,node_index2)
                                                 * wJ;
                    }
                }

                /////////////////////////////////////////////////////
                // Preconditioner matrix
                // Fill the mass matrix (ie \intgl phi_i phi_j) in the
                // pressure-pressure block. Note, the rest of the
                // entries are filled in at the end
                /////////////////////////////////////////////////////
                for(unsigned vertex_index2=0; vertex_index2< NUM_VERTICES_PER_ELEMENT; vertex_index2++)
                {
                    unsigned index2 =  NUM_NODES_PER_ELEMENT*DIM + vertex_index2;
                    rAElemPrecond(index1,index2) +=   linear_phi(vertex_index)
                                                    * linear_phi(vertex_index2)
                                                    * wJ;
                }
            }
        }
    }


    if (assembleJacobian)
    {
        // Fill in the other blocks of the preconditioner matrix. (This doesn't
        // effect the pressure-pressure block of the rAElemPrecond but the
        // pressure-pressure block of rAElem is zero
        rAElemPrecond = rAElemPrecond + rAElem;
    }
}

template class ImplicitCardiacMechanicsAssembler<2>;
template class ImplicitCardiacMechanicsAssembler<3>;


