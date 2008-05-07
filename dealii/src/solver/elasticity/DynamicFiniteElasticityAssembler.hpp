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


#ifndef DYNAMICFINITEELASTICITYASSEMBLER_HPP_
#define DYNAMICFINITEELASTICITYASSEMBLER_HPP_


// speed notes:
//  - don't use Tensor<4,DIM>, use double[][][][] - factor of about 2 difference
//  - compile using OPTIMISATION!!!!!! - factor of about 200 difference when
//     assembling system!!
//  - make reused variables static rather then repeatedly creating them? or member
//     variables. static may cause problems if a 2d sim is run then a 3d sim is run?
//     or probably not.


// todos: see finite elas assembler

#include "FiniteElasticityAssembler.cpp"   //cpp!



/**
 *  DynamicFiniteElasticityAssembler
 *
 *  Solve a time-dependent incompressible finite elasticity problem. In the Lagrangian
 *  coordinates, this is
 *  \f[
 *  \frac{\partial T_{MN} F^i_M}{\partial X^N} + \rho \g^i = \rho \frac{\partial{v^i}}{\partial t}
 *  \f]
 *  where
 *  \f$ T^{MN} \f$ is the second Piola-Kirchoff stress tensor,
 *  \f$ F^i_M \f$ is the deformation gradient tensor
 *  \f$ \rho \f$ is the body mass density,
 *  \f$ g^i \f$ is the body force per unit volume (eg. gravitational acceleration),
 *  and \f$ v^i \f$ is the velocity
 *
 *  See FiniteElasticityAssembler for time-independent problems.
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
 *  Call SetTimes(), and then Solve() to compute the deformed shape, and rGetDeformedPosition()
 *  (in the base class) to get the solution.
 *
 *  The Newton method is used to solve the nonlinear set of finite element equations.
 *  The default degree of the basis functions is quadratic for displacement and linear
 *  for pressure.
 */
template<unsigned DIM>
class DynamicFiniteElasticityAssembler : public FiniteElasticityAssembler<DIM>
{
private:

    /* Inherited 
        Triangulation<DIM>*  mpMesh;
        FESystem<DIM>        mFeSystem;
        DoFHandler<DIM>      mDofHandler;
    
        SparsityPattern      mSparsityPattern;
    
        SparseMatrix<double> mJacobianMatrix;
        Vector<double>       mCurrentSolution;
        Vector<double>       mResidual;
    
        std::string          mOutputDirectoryFullPath;
 
        FourthOrderTensor<DIM> dTdE;
    
        AbstractIncompressibleMaterialLaw<DIM>*  mpMaterialLaw;
        Vector<double>       mBodyForce;
        double               mDensity;
        const unsigned       PRESSURE_COMPONENT_INDEX;
    
        std::map<unsigned,double> mBoundaryValues;
    
        void AssembleSystem(bool assembleResidual, bool assembleJacobian);
        void ApplyDirichletBoundaryConditions(bool assembleResidual, bool assembleJacobian);
    
        void OutputResults(unsigned newtonIteration);
        SparseMatrix<double> mNumericalJacobianMatrix;
     */
    
    /*< Stored solution at the previous time step */
    Vector<double> mSolutionAtLastTimestep;
    
    /*< Start time */
    double mTstart;
    /*< End time */
    double mTend;
    /*< Timestep */
    double mDt;
    /*< 1.0/timestep */
    double mDtInverse;
    /*< Whether the start time, end time and timestep have been set */
    bool   mTimesSet;
    
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
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                           Vector<double>&                                 elementRhs,
                           FullMatrix<double>&                             elementMatrix,
                           bool                                            assembleResidual,
                           bool                                            assembleJacobian);
                           
                           
public:
    /**
     *  Constructor
     *  
     *  @param pMesh A pointer to a dealii mesh. Note, the mesh must have some surface
     *   elements which have had their boundary indicator set to FIXED_BOUNDARY
     *  @param pMaterialLaw A pointer to an incompressible material law
     *  @bodyForce A vector of size DIM represents the body force (force per unit volume)
     *  @density The mass density. Must be strictly positive
     *  @outputDirectory The output directory, relative the the testoutput directory. If 
     *   empty no output is written
     *  @degreeOfBasesForPosition Degree of the polynomials used for interpolating positions.
     *   Defaults to 2, ie quadratic interpolation 
     *  @degreeOfBasesForPressure Degree of the polynomials used for interpolating pressue.
     *   Defaults to 2, ie linear interpolation
     */
    DynamicFiniteElasticityAssembler(Triangulation<DIM>* pMesh,
                                     AbstractIncompressibleMaterialLaw<DIM>*  pMaterialLaw,
                                     Vector<double> bodyForce,
                                     double density,
                                     std::string outputDirectory,
                                     unsigned degreeOfBasesForPosition=2,
                                     unsigned degreeOfBasesForPressure=1);
    virtual ~DynamicFiniteElasticityAssembler();
    
    /**
     *  Set the start and end times, and dt, for the simulation. Must be called 
     *  before Solve()
     */
    void SetTimes(double Tstart, double Tend, double dt);
    
    /**
     *  Solve the dynamic finite elasticity problem. Note, SetTimes() must be called 
     *  before this. rGetDeformedPosition() in the base class can be used to get the 
     *  deformed positions post-solve.
     */
    void Solve();
    
    
    /* Inherited
        void SetBoundaryValues(std::map<unsigned, double> boundary_values);
        void ComputeNumericalJacobian();
        void CompareJacobians();
    */
};




#endif /*DYNAMICFINITEELASTICITYASSEMBLER_HPP_*/
