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


#ifndef FINITEELASTICITYASSEMBLER_HPP_
#define FINITEELASTICITYASSEMBLER_HPP_


#include "AbstractElasticityAssembler.hpp"

/* Some speed notes:
 *  1. Don't use Tensor<4,DIM>, use double[][][][] - factor of about 2 difference
 *  2. Compile using OPTIMISATION!! - factor of about 200 difference when
 *     assembling system!
 *  3. Make reused variables static rather then repeatedly creating them? or member
 *     variables. static may cause problems if a 2d sim is run then a 3d sim is run?
 *     or probably not.
 */


// TODO: better tests against other code. esp 2d or against other code using quadratics

// fix heterogeneity: !form Initial guess just works with one mat law!
// then cover heterogeniety

// choose newton tolerances better.

// refactor (and test) WriteStresses

// nonzero neumann
// refactor out the newton solver?
// chaste style output
// change quad rule if linears


#include <fstream>
#include <iostream>
#include <sstream>
#include <base/symmetric_tensor.h>

const unsigned FIXED_BOUNDARY = 10;
const unsigned NEUMANN_BOUNDARY = 11;
const unsigned DIRICHLET_BOUNDARY = 12;

const double NEWTON_ABS_TOL = 1e-8;
const double NEWTON_REL_TOL = 1e-6;

#include "AbstractIncompressibleMaterialLaw.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
#include "DofVertexIterator.hpp"



/**
 *  FiniteElasticityAssembler
 *
 *  Solve a static incompressible finite elasticity problem. In the Lagrangian
 *  coordinates, this is
 *  \f[
 *  \partial{T_{MN} F^i_M}{X^N} + \rho \g^i = 0
 *  \f]
 *  where
 *  \f$ T^{MN} \f$ is the second Piola-Kirchoff stress tensor,
 *  \f$ F^i_M \f$ is the deformation gradient tensor
 *  \f$ \rho \f$ is the body mass density,
 *  and \f$ g^i \f$ is the body force per unit volume (eg. gravitational acceleration).
 *
 *  See DynamicFiniteElasticityAssembler for time-dependent problems.
 *
 *  The mesh, material law for the material, body force and density are passed in in the
 *  constructor. Alternative methods are available in order to specify material laws
 *  on different regions of the mesh.
 *
 *  Note that the mesh must have some surface elements with their boundary indicator
 *  set to FIXED_BOUNDARY. Zero-valued dirichlet boundary conditions will be specified
 *  on this region.
 *
 *  Neumann boundary conditions are not implemented yet.
 *
 *  Call Solve() to compute the deformed shape.
 *
 *  The Newton method is used to solve the nonlinear set of finite element equations.
 *  The default degree of the basis functions is quadratic for displacement and linear
 *  for pressure.
 *
 *  Note: Calling Solve() repeatedly will use the previous solution as the starting guess.
 */
template<unsigned DIM>
class FiniteElasticityAssembler : public AbstractElasticityAssembler<DIM>
{
protected:
    // note that this must be defined before mDofHandler. except this doesn't
    // matter now that the dof handler is in the abstract class
    /*< The dealii finite element object */
    FESystem<DIM>        mFeSystem;

    /** The derivative of stress (ie second derivative of the strain energy).
     *  dTdE[M][N][P][Q] = d(T^{MN})/d(E_{PQ}) */
    FourthOrderTensor<DIM> dTdE;

    /*< Whether the material is heterogeneous or not */
    bool                 mHeterogeneous;
    /*< The material laws at each region of the mesh */
    std::vector<AbstractIncompressibleMaterialLaw<DIM>*>  mMaterialLaws;
    /*< Map from region number of material law (needed for heterogeneous simulations */
    std::vector<int>     mMaterialIdToMaterialLawIndexMap;
    /*< Helper function for getting material law given element */
    AbstractIncompressibleMaterialLaw<DIM>* GetMaterialLawForElement(typename DoFHandler<DIM>::active_cell_iterator elementIter);

    /*< Body force per unit volume */
    Vector<double>       mBodyForce;
    /*< Mass density of the material (currently as if homogeneous material) */
    double               mDensity;

    /*< just set to be DIM, ie if DIM==2 the spatial indices are 0 and 1, the pressure index is 2 */
    const unsigned       PRESSURE_COMPONENT_INDEX;

    /*< Map from degree of freedom index to boundary value on that degree of freedom */
    std::map<unsigned,double> mBoundaryValues;
    /*< Storage for a numerical approximation of the Jacobian (mainly used in testing) */
    SparseMatrix<double> mNumericalJacobianMatrix;


    /**
     *  Since this assembler will be used repeatedly in quasi-static simulations (eg cardiac,
     *  growth), we want to only call FormInitialGuess() (which guess the zero deformation
     *  solution) the first time, and use the current solution as the guess the rest of the time.
     *  This bool is used for this.
     */
    bool mADeformedHasBeenSolved;

    virtual void WriteStresses(unsigned counter);

    /**
     *  AssembleOnElement
     *
     *  Assemble of the element matrix and/or element vector for the given element. Called
     *  by AssembleSystem in the base class
     *
     *  @elementIter Iterator pointing at current element
     *  @elementRhs Small vector to be filled in. Should be of size AbstractDealiiAssembler::mDofsPerElement
     *  @elementMatrix Small matrix to be filled in. Should be of square, of size AbstractDealiiAssembler::mDofsPerElement
     *  @assembleResidual Whether to assemble the small vector
     *  @assembleJacobian Whether to assemble the small matrix
     */
    virtual void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                                   Vector<double>&                                 elementRhs,
                                   FullMatrix<double>&                             elementMatrix,
                                   bool                                            assembleResidual,
                                   bool                                            assembleJacobian);

    /**
     *  Apply boundary using the mBoundaryValues map to the system matrix and system
     *  vector. Takes into account the fact that this is a nonlinear problem, not a
     *  linear problem, so calculates the appropriate boundary conditions on the linear
     *  sub-problem for the update vector.
     */
    void ApplyDirichletBoundaryConditions();


    /**
     *  Gets the material law corresponding to the required region of the mesh.
     *  Reads mMaterialIdToMaterialLawIndexMap but with error checking
     */
    unsigned GetMaterialLawIndexFromMaterialId(unsigned materialId);

    /**
     *  Set up mCurrentSolution as the solution of zero body force problem. This
     *  is obviously zero displacement, but has non-zero pressure.
     */
    void FormInitialGuess();

    /**
     *  Set up a numerical approximation to the jacobian. Won't generally be needed.
     */
    void ComputeNumericalJacobian();

    /**
     *  Simple method called by the base class which needs to be implemented in
     *  this concrete class
     */
    void DistributeDofs();


public:
    /**
     *  Constructor
     *
     *  @param pMesh A pointer to a dealii mesh. Note, the mesh must have some surface
     *   elements which have had their boundary indicator set to FIXED_BOUNDARY
     *  @param pMaterialLaw A pointer to an incompressible material law. If this is null
     *   SetMaterialLawsForHeterogeneousProblem() must be called before Solve()
     *  @bodyForce A vector of size DIM represents the body force (force per unit volume)
     *  @density The mass density. Must be strictly positive
     *  @outputDirectory The output directory, relative the the testoutput directory. If
     *   empty no output is written
     *  @degreeOfBasesForPosition Degree of the polynomials used for interpolating positions.
     *   Defaults to 2, ie quadratic interpolation
     *  @degreeOfBasesForPressure Degree of the polynomials used for interpolating pressue.
     *   Defaults to 2, ie linear interpolation
     */
    FiniteElasticityAssembler(Triangulation<DIM>* pMesh,
                              AbstractIncompressibleMaterialLaw<DIM>*  pMaterialLaw,
                              Vector<double> bodyForce,
                              double density,
                              std::string outputDirectory,
                              unsigned degreeOfBasesForPosition=2,
                              unsigned degreeOfBasesForPressure=1);

    virtual ~FiniteElasticityAssembler();



    // Note: this type of function doesn't really work
    // void SetDisplacementBoundaryConditions(std::vector<unsigned> nodes,
    //                                        std::vector<unsigned> coordinates,
    //                                        std::vector<double>   values);

    /**
     *  Setting boundary conditions is a hassle. Currently, assuming the
     *  default boundary conditions are not required, the user has to set
     *  up the dof->value map themselves and pass it in using this method.
     *  Note: call rGetDofHandler() to get the dof handler first.
     */
    void SetBoundaryValues(std::map<unsigned, double> boundary_values);


    /** Set the material laws
     *
     *  @materialLaws The material laws of the different regions of the mesh
     *  @materialIds The values of the labels of the different region. Each element in
     *   the mesh should have it's material id set to a value in this vector
     */
    void SetMaterialLawsForHeterogeneousProblem(std::vector<AbstractIncompressibleMaterialLaw<DIM>*> materialLaws,
                                                std::vector<unsigned> materialIds);

    /**
     *  Solve the static finite elasticity problem
     *
     *  @param whether to write output (which will be the solution at the end of every
     *  Newton iteration) to the output directory (if one exists). Defaults to true.
     */
    virtual void StaticSolve(bool writeOutput=true, double tol = -1.0);


    /**
     *  Verify that the analytic jacobian is the same as the numerical jacobian
     *
     *  @param tol. The tolerance with which to compare the absolute component-wise
     *  difference between the two jacobians.
     */
    void CompareJacobians(double tol=1e-8);
};




#endif /*FINITEELASTICITYASSEMBLER_HPP_*/
