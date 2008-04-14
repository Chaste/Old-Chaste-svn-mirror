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

#ifndef IMPLICITCARDIACMECHANICSASSEMBLER_HPP_
#define IMPLICITCARDIACMECHANICSASSEMBLER_HPP_

#include "CardiacMechanicsAssembler.cpp"
#include "NhsSystemWithImplicitSolver.hpp"
#include "LogFile.hpp"
#include "FiniteElasticityTools.hpp"
#include <cfloat>


template<unsigned DIM> 
class ImplicitCardiacMechanicsAssembler : public CardiacMechanicsAssembler<DIM>
{
friend class TestImplicitCardiacMechanicsAssembler;
    
private:
    std::vector<NhsSystemWithImplicitSolver> mCellMechSystems; 
    std::vector<double> mLambdaLastTimeStep;

    double mCurrentTime;
    double mNextTime;
    double mOdeTimestep;
    
public:
    /**
     *  Constructor
     *  
     *  @param pMesh. A pointer to the mesh. Should have a surface set as the fixed surface
     *  @param outputDirectory. The output directory, relative to TEST_OUTPUT
     *  @param pMaterialLaw. The material law for the tissue. Defaults to NULL, in which case
     *   a default material law is used.
     */
    ImplicitCardiacMechanicsAssembler(Triangulation<DIM>* pMesh, 
                                      std::string outputDirectory,
                                      AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw = NULL)
        : CardiacMechanicsAssembler<DIM>(pMesh, outputDirectory, pMaterialLaw),
          mCurrentTime(DBL_MAX),
          mNextTime(DBL_MAX),
          mOdeTimestep(DBL_MAX)
    {
        mCellMechSystems.resize(this->mTotalQuadPoints);
        mLambdaLastTimeStep.resize(this->mTotalQuadPoints, 1.0);        

        // we don't need this active tension in the implicit method. clear to
        // free memory and more importantly to force an error if someone accidentally 
        // uses it in this class
        this->mActiveTension.resize(mCellMechSystems.size());
    }

    ~ImplicitCardiacMechanicsAssembler()
    {
    }
    
    /**
     *  Overloaded SetForcingQuantity(), expecting Calcium concentrations not
     *  active tensions
     */
    void SetForcingQuantity(std::vector<double>& caI)
    {
        assert(caI.size() == this->mTotalQuadPoints);
        for(unsigned i=0; i<caI.size(); i++)
        {
            mCellMechSystems[i].SetIntracellularCalciumConcentration(caI[i]);
        } 
    }

    /** 
     *  Overloaded Solve, which stores the time info, calls the base Solve(),
     *  then updates cell mechanics systems and lambda
     */
    void Solve(double currentTime, double nextTime, double odeTimestep)
    {
        this->mActiveTension.clear();
        
        assert(currentTime < nextTime);
        mCurrentTime = currentTime;
        mNextTime = nextTime;
        mOdeTimestep = odeTimestep;
        
        CardiacMechanicsAssembler<DIM>::Solve(currentTime,nextTime,odeTimestep);

        this->AssembleSystem(true,false);

        for(unsigned i=0; i<mCellMechSystems.size(); i++)
        {
             mCellMechSystems[i].UpdateStateVariables();
             mLambdaLastTimeStep[i] = mCellMechSystems[i].GetLambda();
        }
    }

    void SetFibreSheetMatrix(Tensor<2,DIM> fibreSheetMat)
    {
        EXCEPTION("ImplicitCardiacMechanicsAssembler can't do different fibre directions yet");
    }
    
private:

    /**
     *  AssembleOnElement
     * 
     *  FILL IN     
     */
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                           Vector<double>&       elementRhs,
                           FullMatrix<double>&   elementMatrix,
                           bool                  assembleResidual,
                           bool                  assembleJacobian
                          )
    {
        // check these have been set
        assert(mCurrentTime != DBL_MAX);
        assert(mNextTime != DBL_MAX);
        assert(mOdeTimestep != DBL_MAX);
        
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
        
        const unsigned n_q_points = quadrature_formula.n_quadrature_points;
        
        
        // would want this to be static too (slight speed up), but causes errors
        // in debug mode (upon destruction of the class, in 2d, or something)
        FEValues<DIM> fe_values(this->mFeSystem, quadrature_formula,
                                UpdateFlags(update_values    |
                                            update_gradients |
                                            update_q_points  |     // needed for interpolating u and u' on the quad point
                                            update_JxW_values));
                                                                                             
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
    //    this->mFibreSheetMat[0][0] =  cos(theta);
    //    this->mFibreSheetMat[0][1] =  sin(theta);
    //    this->mFibreSheetMat[1][0] = -sin(theta);
    //    this->mFibreSheetMat[1][1] =  cos(theta);
    //    this->mTransFibreSheetMat = transpose(this->mFibreSheetMat);
        
        
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
    
            // get the active tension at this quad point
            //double active_tension = mActiveTension[this->mCurrentQuadPointGlobalIndex];
    
            static Tensor<2,DIM> C_fibre;          // C when transformed to fibre-sheet axes
            static Tensor<2,DIM> inv_C_fibre;      // C^{-1} transformed to fibre-sheet axes
            static SymmetricTensor<2,DIM> T_fibre; // T when transformed to fibre-sheet axes
            
            // transform C and invC
            C_fibre = this->mTransFibreSheetMat * C * this->mFibreSheetMat;
            inv_C_fibre = this->mTransFibreSheetMat * inv_C * this->mFibreSheetMat;
    
            // store the stretch in the fibre direction
            this->mLambda[this->mCurrentQuadPointGlobalIndex] = sqrt(C_fibre[0][0]);
    

            double lam = sqrt(C_fibre[0][0]);
            double dlam_dt = (lam-mLambdaLastTimeStep[this->mCurrentQuadPointGlobalIndex])/(mNextTime-mCurrentTime);

            NhsSystemWithImplicitSolver& system = mCellMechSystems[this->mCurrentQuadPointGlobalIndex];

            // get proper active tension
            // see NOTE below
            system.SetLambdaAndDerivative(lam, dlam_dt);

            double active_tension;        
            try
            {
                system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
                active_tension = system.GetActiveTensionAtNextTime()/this->mScaleFactor;        
            }
            catch (Exception& e)
            {
                LOG(1, "\nCAUGHT EXCEPTION!!\n");
                active_tension = 1e10;
                if(assembleJacobian) EXCEPTION("Failed");
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
                double active_tension_at_lam_plus_h = system.GetActiveTensionAtNextTime()/this->mScaleFactor;        

                // get active tension for (lam,dlamdt+h)
                double h2 = std::max(1e-6, dlam_dt/100);
                system.SetLambdaAndDerivative(lam, dlam_dt+h2);
                system.SolveDoNotUpdate(mCurrentTime,mNextTime,mOdeTimestep);
                double active_tension_at_dlamdt_plus_h = system.GetActiveTensionAtNextTime()/this->mScaleFactor;        

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
                assert( fabs(system.GetActiveTensionAtNextTime()/this->mScaleFactor-active_tension)<1e-8);
            }
            catch (Exception& e)
            {        
                LOG(1, "WARNING!\n");
                active_tension = 1e10;
                // should have done something abouve..
            }

            //this->mDTdE_fibre.Zero();
    
            // compute the transformed tension. The material law should law be a cardiac-
            // specific law which assumes the x-axes in the fibre, the z-axes the sheet normal
            p_material_law->ComputeStressAndStressDerivative(C_fibre,inv_C_fibre,p,T_fibre,this->mDTdE_fibre,assembleJacobian);

            // amend the stress and dTdE using the active tension
            T_fibre[0][0] += active_tension/C_fibre[0][0];
            this->mDTdE_fibre(0,0,0,0) -= 2*active_tension/(C_fibre[0][0]*C_fibre[0][0]);  

            
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
                            T[M][N] +=                      T_fibre [al][be]
                                        *      this->mFibreSheetMat [M] [al]
                                        * this->mTransFibreSheetMat [be][N];
                        }
                    }
                }
            }            
    
    
            static FourthOrderTensor<DIM> temp1;
            static FourthOrderTensor<DIM> temp2;
            static FourthOrderTensor<DIM> temp3;
    
            temp1.SetAsProduct(this->mDTdE_fibre, this->mFibreSheetMat, 0);
            temp2.SetAsProduct(temp1,             this->mFibreSheetMat, 1);
            temp3.SetAsProduct(temp2,             this->mFibreSheetMat, 2);
            
            this->dTdE.SetAsProduct(temp3, this->mFibreSheetMat, 3);
    
            
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
                                }
                            }
                            
                            // extra bit to the matrix coming from differentiating
                            // Ta wrt U_I (Ta is dependent of new nodal positions 
                            // U_I as Ta is dependent on x through lam=C00 and lam_dot
                            // We don't have to use the term coming from differentiated
                            // the (1/C00) bit, that is accounted for in dTdE
                            elementMatrix(i,j) +=   (
                                                       d_act_tension_dlam
                                                     +
                                                       d_act_tension_d_dlamdt/(mNextTime-mCurrentTime)
                                                    )
                                                    * (F[component_j][0]/lam)
                                                    * (fe_values.shape_grad(j,q_point)[0]/C[0][0])
                                                    * F[component_i][0]
                                                    * fe_values.shape_grad(i,q_point)[0]
                                                    * fe_values.JxW(q_point);
                                    

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
        
        first = false;
    }



//// DONT DELETE
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

#endif /*IMPLICITCARDIACMECHANICSASSEMBLER_HPP_*/
