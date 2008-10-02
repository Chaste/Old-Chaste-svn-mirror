/*

Copyright (C) University of Oxford, 2008

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


#ifndef IMPLICITCARDIACMECHANICSASSEMBLER2_HPP_
#define IMPLICITCARDIACMECHANICSASSEMBLER2_HPP_

#include "NonlinearElasticityAssembler.hpp"
#include "NhsSystemWithImplicitSolver.hpp"
#include "NashHunterPoleZeroLaw2.hpp"
#include "LogFile.hpp"
#include <cfloat>

// NOTE: with zero body force and quasi-static elasticity the density is taken in but
// not used (as multiplied by zero body force), hence we just use the value 1.0.
const double CARDIAC_TISSUE_DENSITY = 1.0;  


/**
 *  Implicit Cardiac Mechanics Assembler
 * 
 *  Solves cardiac mechanics implicitly (together with the NHS cell
 *  models for determining the active tension), taking in the intracellular
 *  Calcium concentration. See CardiacElectroMechanicsProblem2 documentation
 *  for more detail. 
 */
template<unsigned DIM>
class ImplicitCardiacMechanicsAssembler2 : public NonlinearElasticityAssembler<DIM>
{
friend class TestImplicitCardiacMechanicsAssembler2;

private:
    /** 
     *  The NHS cell systems (with their own implicit solvers, which take in 
     *  [Ca]_i and return Ta. Note the indexing: the i-th entry corresponds to
     *  the i-th global quad point, when looping over elements and then
     *  quad points */
    std::vector<NhsSystemWithImplicitSolver> mCellMechSystems;
    
    /** The stretch ratio (in the fibre direction) at the last timestep.
     *  Note the indexing: the i-th entry corresponds to the i-th global 
     *  quad point, when looping over elements and then quad points 
     */ 
    std::vector<double> mLambdaLastTimeStep;
    
    /** The current stretch ratio (in the fibre direction). Note the indexing: 
     *  the i-th entry corresponds to the i-th global quad point, when looping 
     *  over elements and then quad points 
     */ 
    std::vector<double> mLambda;

    /*< Current time */
    double mCurrentTime;
    /*< Time to which the solver has been asked to solve to */
    double mNextTime;
    /*< Time used to integrate the NHS model */
    double mOdeTimestep;

    /*< Whether the material law was passed in or the default used */
    bool mAllocatedMaterialLawMemory;
    
    /*< Total number of quad points in the (mechanics) mesh */
    unsigned mTotalQuadPoints;
    
public:
    /**
     *  Constructor
     *
     *  @param pMesh. A pointer to the mesh. Should have a surface set as the fixed surface
     *  @param outputDirectory. The output directory, relative to TEST_OUTPUT
     *  @param pMaterialLaw. The material law for the tissue. Defaults to NULL, in which case
     *   a default material law is used.
     */
    ImplicitCardiacMechanicsAssembler2(QuadraticMesh<DIM>* pQuadMesh,
                                       std::string outputDirectory,
                                       std::vector<unsigned>& rFixedNodes,
                                       AbstractIncompressibleMaterialLaw2<DIM>* pMaterialLaw = NULL)
        : NonlinearElasticityAssembler<DIM>(pQuadMesh, 
                                            pMaterialLaw!=NULL ? pMaterialLaw : new NashHunterPoleZeroLaw2<DIM>,
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
    
    /** 
     *  Destructor just deletes memory if it was allocated
     */
    ~ImplicitCardiacMechanicsAssembler2()
    {
        if(mAllocatedMaterialLawMemory)
        {
            assert(this->mMaterialLaws.size()==1); // haven't implemented heterogeniety yet
            delete this->mMaterialLaws[0];
        }
    }
    
    /*< Get the total number of quad points in the mesh */
    unsigned GetTotalNumQuadPoints()
    {
        return mTotalQuadPoints;
    }
    
    /*< Get the quadrature rule used in the elements */
    GaussianQuadratureRule<DIM>* GetQuadratureRule()
    {
        return this->mpQuadratureRule;
    }

    /**
     *  Set the intracellular Calcium concentrations (note: in an explicit algorithm we 
     *  would set the active tension as the forcing quantity; the implicit algorithm
     *  takes in the Calcium concentration and solves for the active tension implicitly
     *  together with the mechanics.
     */
    void SetIntracellularCalciumConcentrations(std::vector<double>& caI)
    {
        assert(caI.size() == mTotalQuadPoints);
        for(unsigned i=0; i<caI.size(); i++)
        {
            mCellMechSystems[i].SetIntracellularCalciumConcentration(caI[i]);
        }
    }

    /**
     *  Solve for the deformation using quasi-static nonlinear elasticity.
     *  (not dynamic nonlinear elasticity, despite the times taken in - just ONE
     *  deformation is solved for. The cell models are integrated implicitly
     *  over the time range using the ODE timestep provided, as part of the solve,
     *  and updated at the end once the solution has been found, as is lambda.
     */
    void Solve(double currentTime, double nextTime, double odeTimestep)
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


private:

    /**
     *  Overloaded AssembleOnElement. Apart from a tiny bit of initial set up and 
     *  the lack of the body force term in the residual, the bits were this is
     *  different to the base class AssembleOnElement are restricted to two bits
     *  (see code): calculating Ta implicitly and using it to compute the stress,
     *  and the addition of a corresponding extra term to the Jacobian
     */
    void AssembleOnElement(Element<DIM, DIM>& rElement,
                           c_matrix<double, NonlinearElasticityAssembler<DIM>::STENCIL_SIZE, NonlinearElasticityAssembler<DIM>::STENCIL_SIZE >& rAElem,
                           c_vector<double, NonlinearElasticityAssembler<DIM>::STENCIL_SIZE>& rBElem,
                           bool assembleResidual,
                           bool assembleJacobian)
    {
        // check these have been set
        assert(mCurrentTime != DBL_MAX);
        assert(mNextTime != DBL_MAX);
        assert(mOdeTimestep != DBL_MAX);
        
        const c_matrix<double, DIM, DIM>* p_inverse_jacobian = rElement.GetInverseJacobian();
        double jacobian_determinant = rElement.GetJacobianDeterminant();

        if (assembleJacobian)
        {
            rAElem.clear();
        }

        if (assembleResidual)
        {
            rBElem.clear();
        }

        ///////////////////////////////////////////////
        // Get the current displacement at the nodes
        ///////////////////////////////////////////////
        static c_matrix<double,DIM,NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT> element_current_displacements;
        static c_vector<double,NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT> element_current_pressures;
        for(unsigned II=0; II<NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT; II++)
        {
            for(unsigned JJ=0; JJ<DIM; JJ++)
            {
                element_current_displacements(JJ,II) = this->mCurrentSolution[DIM*rElement.GetNodeGlobalIndex(II) + JJ];
            }
        }

        ///////////////////////////////////////////////
        // Get the current pressure at the vertices
        ///////////////////////////////////////////////
        for(unsigned II=0; II<NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT; II++)
        {
            element_current_pressures(II) = this->mCurrentSolution[DIM*this->mpQuadMesh->GetNumNodes() + rElement.GetNodeGlobalIndex(II)];
        }

        // allocate memory for the basis functions values and derivative values
        c_vector<double, NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT> linear_phi;
        c_vector<double, NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT> quad_phi;
        c_matrix<double, DIM, NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT> grad_quad_phi;

        // get the material law
        AbstractIncompressibleMaterialLaw2<DIM>* p_material_law;
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
            QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, *p_inverse_jacobian, grad_quad_phi);
            
            ////////////////////////////////////////////////////
            // (dont get the body force)
            ////////////////////////////////////////////////////

            //////////////////////////////////////
            // interpolate grad_u and p
            //////////////////////////////////////
            static c_matrix<double,DIM,DIM> grad_u; // grad_u = (du_i/dX_M)
            grad_u = zero_matrix<double>(DIM,DIM);  // must be on new line!!

            for(unsigned node_index=0; 
                node_index<NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT; 
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
                vertex_index<NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT;
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

            double active_tension;
            try
            {
                system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
                active_tension = system.GetActiveTensionAtNextTime(); ///this->mScaleFactor;
            }
            catch (Exception& e)
            {
                active_tension = 1e10;
                if(assembleJacobian) 
                {
                    EXCEPTION("Failed in solving NHS systems when assembling Jacobian");
                }
            }

            // compute the derivative of the active tension wrt lam and dlam_dt
            double d_act_tension_dlam;
            double d_act_tension_d_dlamdt;

            if(assembleJacobian)
            {
                // get active tension for (lam+h,dlamdt)
                double h1 = std::max(1e-6, lam/100);
                system.SetLambdaAndDerivative(lam+h1, dlam_dt);
                system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
                double active_tension_at_lam_plus_h = system.GetActiveTensionAtNextTime();       ///this->mScaleFactor;

                // get active tension for (lam,dlamdt+h)
                double h2 = std::max(1e-6, dlam_dt/100);
                system.SetLambdaAndDerivative(lam, dlam_dt+h2);
                system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
                double active_tension_at_dlamdt_plus_h = system.GetActiveTensionAtNextTime(); //this->mScaleFactor;

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
                assert( fabs(system.GetActiveTensionAtNextTime()/*/this->mScaleFactor*/-active_tension)<1e-8);
            }
            catch (Exception& e)
            {
                LOG(2, "WARNING in ImplicitCardiacMechanicsAssembler!\n");
                active_tension = 1e10;
                // should have done something above..
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


            /*************************************
             * end of cardiac specific code PART 1
             *************************************/

            /////////////////////////////////////////
            // residual vector
            /////////////////////////////////////////
            if (assembleResidual)
            {
                for(unsigned index=0; index<NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT*DIM; index++)
                {
                    unsigned spatial_dim = index%DIM;
                    unsigned node_index = (index-spatial_dim)/DIM;

                    assert(node_index < NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT);

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
                
                for(unsigned vertex_index=0; vertex_index<NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    rBElem( NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT*DIM + vertex_index ) +=   linear_phi(vertex_index)
                                                                          * (detF - 1)
                                                                          * wJ;
                }
            }

            /////////////////////////////////////////
            // Jacobian matrix
            /////////////////////////////////////////
            if(assembleJacobian)
            {
                for(unsigned index1=0; index1<NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT*DIM; index1++)
                {
                    unsigned spatial_dim1 = index1%DIM;
                    unsigned node_index1 = (index1-spatial_dim1)/DIM;
                    
                    
                    for(unsigned index2=0; index2<NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT*DIM; index2++)
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

                                for (unsigned P=0; P<DIM; P++)
                                {
                                    for (unsigned Q=0; Q<DIM; Q++)
                                    {
                                        rAElem(index1,index2)  +=   0.5
                                                                  * this->dTdE(M,N,P,Q)
                                                                  * (
                                                                      grad_quad_phi(P,node_index2)
                                                                    * F(spatial_dim2,Q)
                                                                       +
                                                                      grad_quad_phi(Q,node_index2)
                                                                    * F(spatial_dim2,P)
                                                                     )
                                                                  * F(spatial_dim1,M)
                                                                  * grad_quad_phi(N,node_index1)
                                                                  * wJ;
                                    }
                                }
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
                                                    //* (fe_values.shape_grad(j,q_point)[0]/C[0][0])
                                                    * F(spatial_dim1,0)
                                                    * grad_quad_phi(0,node_index1)
                                                    //* fe_values.shape_grad(i,q_point)[0]
                                                    * wJ;                        
                       /************************************
                        *  End cardiac-specific code PART 2
                        ************************************/
                    }
                    
                    
                    for(unsigned vertex_index=0; vertex_index<NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT; vertex_index++)
                    {
                        unsigned index2 = NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT*DIM + vertex_index;
                        
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

                for(unsigned vertex_index=0; vertex_index<NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    unsigned index1 = NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                    for(unsigned index2=0; index2<NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT*DIM; index2++)
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
                }
            }
        }
    }



//// THE FOLLOWING IS IF THE MONODOMAIN EQUATIONS ARE ADJUSTED TO USE inverse(C)
//// (still the dealii version though)
//// DO NOT DELETE
////
////
//public:
//    std::vector<std::vector<unsigned> > mNodesContainedInElement;
//
//    void ComputeElementsContainingNodes(ConformingTetrahedralMesh<DIM,DIM>* pOtherMesh)
//    {
//        assert(DIM==2);
//
//        mNodesContainedInElement.resize(this->mpMesh->n_active_cells());
//
//        unsigned element_number = 0;
//        typename DoFHandler<DIM>::active_cell_iterator  element_iter = this->mDofHandler.begin_active();
//
//        while (element_iter!=this->mDofHandler.end())
//        {
//            double xmin = element_iter->vertex(0)(0);
//            double xmax = element_iter->vertex(1)(0);
//            double ymin = element_iter->vertex(0)(1);
//            double ymax = element_iter->vertex(3)(1);
//
//            assert(element_iter->vertex(2)(0)==xmax);
//            assert(element_iter->vertex(2)(1)==ymax);
//
//            for(unsigned i=0; i<pOtherMesh->GetNumNodes(); i++)
//            {
//                double x = pOtherMesh->GetNode(i)->rGetLocation()[0];
//                double y = pOtherMesh->GetNode(i)->rGetLocation()[1];
//                if((x>=xmin) && (x<=xmax) && (y>=ymin) && (y<=ymax))
//                {
//                    mNodesContainedInElement[element_number].push_back(i);
//                }
//            }
//
//            element_iter++;
//            element_number++;
//        }
//    }
//
//    void WriteLambda(std::string directory, std::string fileName)
//    {
//        OutputFileHandler handler(directory,false);
//        out_stream p_file = handler.OpenOutputFile(fileName);
//
//        std::vector<std::vector<double> > quad_point_posns
//           = FiniteElasticityTools<DIM>::GetQuadPointPositions(*(this->mpMesh), this->GetNumQuadPointsInEachDimension());
//
//
//        for(unsigned i=0; i<quad_point_posns.size(); i++)
//        {
//            (*p_file) << quad_point_posns[i][0] << " " << quad_point_posns[i][1] << " "
//                      << mCellMechSystems[i].GetLambda() << "\n";
//        }
//    }
//
//
//    void CalculateCinverseAtNodes(ConformingTetrahedralMesh<DIM,DIM>* pOtherMesh, std::vector<std::vector<double> >& rValuesAtNodes)
//    {
//        assert(DIM==2);
//        rValuesAtNodes.resize(pOtherMesh->GetNumNodes());
//
//        unsigned element_number = 0;
//
//        static QTrapez<DIM>   trapezoid_quadrature_formula; //trapeziod rule - values at NODES
//        const unsigned n_q_points = trapezoid_quadrature_formula.n_quadrature_points;
//
//        FEValues<DIM> fe_values(this->mFeSystem, trapezoid_quadrature_formula,
//                                UpdateFlags(update_values    |
//                                            update_gradients |
//                                            update_q_points  |     // needed for interpolating u and u' on the quad point
//                                            update_JxW_values));
//
//        std::vector< Vector<double> >                  local_solution_values(n_q_points);
//        std::vector< std::vector< Tensor<1,DIM> > >    local_solution_gradients(n_q_points);
//
//        for (unsigned q_point=0; q_point<n_q_points; q_point++)
//        {
//            local_solution_values[q_point].reinit(DIM+1);
//            local_solution_gradients[q_point].resize(DIM+1);
//        }
//
//
//        Tensor<2,DIM> identity;
//        for (unsigned i=0; i<DIM; i++)
//        {
//            for (unsigned j=0; j<DIM; j++)
//            {
//                identity[i][j] = i==j ? 1.0 : 0.0;
//            }
//        }
//
//        typename DoFHandler<DIM>::active_cell_iterator  element_iter = this->mDofHandler.begin_active();
//
//        while (element_iter!=this->mDofHandler.end())
//        {
//            double xmin = element_iter->vertex(0)(0);
//            double xmax = element_iter->vertex(1)(0);
//            double ymin = element_iter->vertex(0)(1);
//            double ymax = element_iter->vertex(3)(1);
//            assert(element_iter->vertex(2)(0)==xmax);
//            assert(element_iter->vertex(2)(1)==ymax);
//
//            fe_values.reinit(element_iter); // compute fe values for this element
//            fe_values.get_function_values(this->mCurrentSolution, local_solution_values);
//            fe_values.get_function_grads(this->mCurrentSolution, local_solution_gradients);
//
//            std::vector<Point<DIM> > quad_points =fe_values.get_quadrature_points();
//
//
//            AbstractIncompressibleMaterialLaw<DIM>* p_material_law = this->GetMaterialLawForElement(element_iter);
//
//            std::vector<Tensor<2,DIM> > inv_C_at_nodes(4);// 4=2^DIM
//
//            for (unsigned q_point=0; q_point<n_q_points; q_point++)
//            {
//                const std::vector< Tensor<1,DIM> >& grad_u_p = local_solution_gradients[q_point];
//                static Tensor<2,DIM> F;
//                static Tensor<2,DIM> C;
//
//                for (unsigned i=0; i<DIM; i++)
//                {
//                    for (unsigned j=0; j<DIM; j++)
//                    {
//                        F[i][j] = identity[i][j] + grad_u_p[i][j];
//                    }
//                }
//
//                C = transpose(F) * F;
//                inv_C_at_nodes[q_point] = invert(C);
//            }
//
//
///// QUAD POINT ORDER: (0,0), (1,0), (0,1), (1,1)
////            std::cout << quad_points[0](0) << " " << quad_points[0](1) << "\n";
////            std::cout << quad_points[1](0) << " " << quad_points[1](1) << "\n";
////            std::cout << quad_points[2](0) << " " << quad_points[2](1) << "\n";
////            std::cout << quad_points[3](0) << " " << quad_points[3](1) << "\n";
////            std::cout << xmin << " " << ymin << " " << local_solution_values[0](0) << "\n";
////            std::cout << xmin << " " << ymax << " " << local_solution_values[1](0) << "\n";
////            std::cout << xmax << " " << ymin << " " << local_solution_values[2](0) << "\n";
////            std::cout << xmax << " " << ymax << " " << local_solution_values[3](0) << "\n";
//
//
//
//            for(unsigned j=0; j<mNodesContainedInElement[element_number].size(); j++)
//            {
//                unsigned node_num = mNodesContainedInElement[element_number][j];
//                double x = pOtherMesh->GetNode(node_num)->rGetLocation()[0];
//                double y = pOtherMesh->GetNode(node_num)->rGetLocation()[1];
//
//                assert((x>=xmin) && (x<=xmax) && (y>=ymin) && (y<=ymax));
//                double xi  = (x-xmin)/(xmax-xmin);
//                double eta = (y-ymin)/(ymax-ymin);
//                assert((0<=xi) && (x<=1) && (0<=eta) && (eta<=1));
//
//                rValuesAtNodes[node_num][0] = InterpolateCinverse(xi,eta,inv_C_at_nodes,0,0);
//                rValuesAtNodes[node_num][1] = InterpolateCinverse(xi,eta,inv_C_at_nodes,0,1);
//                rValuesAtNodes[node_num][2] = InterpolateCinverse(xi,eta,inv_C_at_nodes,1,1);
//            }
//
//
//            element_iter++;
//            element_number++;
//        }
//    }
//
//
//    double InterpolateCinverse(const double xi, const double eta,
//                               const std::vector<Tensor<2,DIM> >& inverseCAtNodes,
//                               unsigned i, unsigned j)
//    {
//        return    inverseCAtNodes[0][i][j] * (1-xi) * (1-eta)
//                + inverseCAtNodes[1][i][j] * (1-xi) *   eta
//                + inverseCAtNodes[2][i][j] *   xi   * (1-eta)
//                + inverseCAtNodes[3][i][j] *   xi   *   eta;
//    }
 };

#endif /*IMPLICITCARDIACMECHANICSASSEMBLER2_HPP_*/
