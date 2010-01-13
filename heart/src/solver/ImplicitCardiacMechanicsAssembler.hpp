/*

Copyright (C) University of Oxford, 2005-2010

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


#ifndef IMPLICITCARDIACMECHANICSASSEMBLER_HPP_
#define IMPLICITCARDIACMECHANICSASSEMBLER_HPP_

#include "AbstractCardiacMechanicsAssembler.hpp"
#include "QuadraticBasisFunction.hpp"
#include "LinearBasisFunction.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "AbstractContractionModel.hpp"
#include "LogFile.hpp"
#include <cfloat>


/**
 *  Implicit Cardiac Mechanics Assembler
 *
 *  Solves cardiac mechanics implicitly (together with the contraction
 *  models for determining the active tension), taking in the intracellular
 *  calcium concentration. See CardiacElectroMechanicsProblem documentation
 *  for more detail.
 */
template<unsigned DIM>
class ImplicitCardiacMechanicsAssembler : public AbstractCardiacMechanicsAssembler<DIM>
{
friend class TestImplicitCardiacMechanicsAssembler;

private:
    /** The stretch in the fibre direction at the last timestep, in order
     *  to compute the stretch rate.
     *  Note the indexing: the i-th entry corresponds to the i-th global
     *  quad point, when looping over elements and then quad points
     */
    std::vector<double> mStretchesLastTimeStep;

    /** The current stretch ratio (in the fibre direction). Note the indexing:
     *  the i-th entry corresponds to the i-th global quad point, when looping
     *  over elements and then quad points
     */

    /** This solver is an implicit solver (overloaded pure method) */
    bool IsImplicitSolver()
    {
        return true;
    }
    /**
     *  A method called by AbstractCardiacMechanicsAssembler::AssembleOnElement(), providing
     *  the active tension (and other info) at a particular quadrature point. This version uses C to 
     *  determine the current stretch and stretch rate, and integrates the contraction model ODEs to determine
     *  the active tension, and derivatives of the active tension with respect to stretch and
     *  stretch rate.
     * 
     *  @param currentFibreStretch The stretch in the fibre direction
     *  @param currentQuadPointGlobalIndex Quadrature point integrand currently being evaluated at in AssembleOnElement.
     *  @param assembleJacobian  A bool stating whether to assemble the Jacobian matrix.
     *  @param rActiveTension The returned active tension
     *  @param rDerivActiveTensionWrtLambda The returned dT_dLam, derivative of active tension wrt stretch
     *  @param rDerivActiveTensionWrtDLambdaDt The returned dT_dLamDot, derivative of active tension wrt stretch rate
     */
    void GetActiveTensionAndTensionDerivs(double currentFibreStretch,  
                                          unsigned currentQuadPointGlobalIndex,
                                          bool assembleJacobian,
                                          double& rActiveTension,
                                          double& rDerivActiveTensionWrtLambda,
                                          double& rDerivActiveTensionWrtDLambdaDt);

public:
    /**
     * Constructor
     *
     * @param contractionModel The contraction model.
     * @param pQuadMesh A pointer to the mesh.
     * @param outputDirectory The output directory, relative to TEST_OUTPUT
     * @param rFixedNodes The fixed nodes
     * @param pMaterialLaw The material law for the tissue. Defaults to NULL, in which case
     *   a default material law is used.
     */
    ImplicitCardiacMechanicsAssembler(ContractionModel contractionModel,
                                      QuadraticMesh<DIM>* pQuadMesh,
                                      std::string outputDirectory,
                                      std::vector<unsigned>& rFixedNodes,
                                      AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw = NULL);

    /**
     *  Destructor
     */
    virtual ~ImplicitCardiacMechanicsAssembler();


    /**
     *  Get lambda (the stretch ratio).
     *  NOTE: the i-th entry of this vector is assumed to be the i-th quad point
     *  obtained by looping over cells in the obvious way and then looping over
     *  quad points. These quad points, in the same order, can be obtained by
     *  using the QuadraturePointsGroup class.
     */
    std::vector<double>& rGetFibreStretches();

    /**
     *  Solve for the deformation using quasi-static nonlinear elasticity.
     *  (not dynamic nonlinear elasticity, despite the times taken in - just ONE
     *  deformation is solved for. The cell models are integrated implicitly
     *  over the time range using the ODE timestep provided, as part of the solve,
     *  and updated at the end once the solution has been found, as is lambda.
     * 
     *  @param time the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    void Solve(double time, double nextTime, double odeTimestep);
};


//// THE FOLLOWING IS IF THE MONODOMAIN EQUATIONS ARE ADJUSTED TO USE inverse(C)
//// (still the dealii version though)
////
////  *** DO NOT DELETE ***
////
////
//public:
//    std::vector<std::vector<unsigned> > mNodesContainedInElement;
//
//    void ComputeElementsContainingNodes(TetrahedralMesh<DIM,DIM>* pOtherMesh)
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
//    void CalculateCinverseAtNodes(TetrahedralMesh<DIM,DIM>* pOtherMesh, std::vector<std::vector<double> >& rValuesAtNodes)
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


#endif /*IMPLICITCARDIACMECHANICSASSEMBLER_HPP_*/
