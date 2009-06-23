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


#ifndef IMPLICITCARDIACMECHANICSASSEMBLER_HPP_
#define IMPLICITCARDIACMECHANICSASSEMBLER_HPP_

#include "NonlinearElasticityAssembler.hpp"
#include "QuadraticBasisFunction.hpp"
#include "LinearBasisFunction.hpp"
#include "NhsSystemWithImplicitSolver.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "LogFile.hpp"
#include <cfloat>


/**
 *  Implicit Cardiac Mechanics Assembler
 *
 *  Solves cardiac mechanics implicitly (together with the NHS cell
 *  models for determining the active tension), taking in the intracellular
 *  Calcium concentration. See CardiacElectroMechanicsProblem documentation
 *  for more detail.
 */
template<unsigned DIM>
class ImplicitCardiacMechanicsAssembler : public NonlinearElasticityAssembler<DIM>
{
friend class TestImplicitCardiacMechanicsAssembler;

private:
    static const unsigned STENCIL_SIZE = NonlinearElasticityAssembler<DIM>::STENCIL_SIZE;
    static const unsigned NUM_NODES_PER_ELEMENT = NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT;
    static const unsigned NUM_VERTICES_PER_ELEMENT = NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT;
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

    /** Current time */
    double mCurrentTime;
    /** Time to which the solver has been asked to solve to */
    double mNextTime;
    /** Time used to integrate the NHS model */
    double mOdeTimestep;

    /** Whether the material law was passed in or the default used */
    bool mAllocatedMaterialLawMemory;

    /** Total number of quad points in the (mechanics) mesh */
    unsigned mTotalQuadPoints;

public:
    /**
     * Constructor
     *
     * @param pQuadMesh A pointer to the mesh. Should have a surface set as the fixed surface
     * @param outputDirectory The output directory, relative to TEST_OUTPUT
     * @param rFixedNodes The fixed nodes
     * @param pMaterialLaw The material law for the tissue. Defaults to NULL, in which case
     *   a default material law is used.
     */
    ImplicitCardiacMechanicsAssembler(QuadraticMesh<DIM>* pQuadMesh,
                                      std::string outputDirectory,
                                      std::vector<unsigned>& rFixedNodes,
                                      AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw = NULL);

    /**
     *  Destructor just deletes memory if it was allocated
     */
    ~ImplicitCardiacMechanicsAssembler();

    /** Get the total number of quad points in the mesh */
    unsigned GetTotalNumQuadPoints();

    /** Get the quadrature rule used in the elements */
    GaussianQuadratureRule<DIM>* GetQuadratureRule();

    /**
     *  Set the intracellular Calcium concentrations (note: in an explicit algorithm we
     *  would set the active tension as the forcing quantity; the implicit algorithm
     *  takes in the Calcium concentration and solves for the active tension implicitly
     *  together with the mechanics).
     * 
     *  @param caI the intracellular calcium concentrations
     */
    void SetIntracellularCalciumConcentrations(std::vector<double>& caI);

    /**
     *  Get lambda (the stretch ratio).
     *  NOTE: the i-th entry of this vector is assumed to be the i-th quad point
     *  obtained by looping over cells in the obvious way and then looping over
     *  quad points. These quad points, in the same order, can be obtained by
     *  using the QuadraturePointsGroup class.
     */
    std::vector<double>& rGetLambda();

    /**
     *  Solve for the deformation using quasi-static nonlinear elasticity.
     *  (not dynamic nonlinear elasticity, despite the times taken in - just ONE
     *  deformation is solved for. The cell models are integrated implicitly
     *  over the time range using the ODE timestep provided, as part of the solve,
     *  and updated at the end once the solution has been found, as is lambda.
     * 
     *  @param currentTime the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    void Solve(double currentTime, double nextTime, double odeTimestep);


private:

    /**
     * Overloaded AssembleOnElement. Apart from a tiny bit of initial set up and
     * the lack of the body force term in the residual, the bits where this is
     * different to the base class AssembleOnElement are restricted to two bits
     * (see code): calculating Ta implicitly and using it to compute the stress,
     * and the addition of a corresponding extra term to the Jacobian.
     * 
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     * @param rAElemPrecond \todo Document this parameter
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    void AssembleOnElement(Element<DIM, DIM>& rElement,
                           c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElem,
                           c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElemPrecond,
                           c_vector<double,STENCIL_SIZE>& rBElem,
                           bool assembleResidual,
                           bool assembleJacobian);
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
