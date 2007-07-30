#ifndef FINITEELASTICITYASSEMBLER_HPP_
#define FINITEELASTICITYASSEMBLER_HPP_


#include "AbstractDealiiAssembler.hpp"

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

const double NEWTON_ABS_TOL = 1e-13;
const double NEWTON_REL_TOL = 1e-7;

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
 */
template<unsigned DIM>
class FiniteElasticityAssembler : public AbstractDealiiAssembler<DIM>
{
protected:
    // note that this must be defined before mDofHandler. except this doesn't
    // matter now that the dof handler is in the abstract class
    /*< The dealii finite element object */
    FESystem<DIM>        mFeSystem;
    
    /*< Full path of output directory, including chaste testoutput */
    std::string          mOutputDirectoryFullPath;
    /*< Whether to write any output or not */
    bool                 mWriteOutput;
    
    /** The derivative of stress (ie second derivative of the strain energy).
     *  dTdE[M][N][P][Q] = d(T^{MN})/d(E_{PQ}) */
    double               dTdE[DIM][DIM][DIM][DIM];
    
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
    
    /*< Number of newton iterations needed to solve the problem, for interest and testing */
    unsigned mNumNewtonIterations;
    
    /*< Map from degree of freedom index to boundary value on that degree of freedom */
    std::map<unsigned,double> mBoundaryValues;
    /*< Storage for a numerical approximation of the Jacobian (mainly used in testing) */
    SparseMatrix<double> mNumericalJacobianMatrix;
    
    /*< Data structure containing the deformed position, by vertex index, in easily
     * accessable form. Only created if asked for */
    std::vector<Vector<double> > mDeformedPosition;
    /*< Data structure containing the undeformed position, by vertex index, in easily
     * accessable form. Only created if asked for */
    std::vector<Vector<double> > mUndeformedPosition;
    
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
     *  Compute the L2 norm of the current residual vector divided by it's length.
     */
    double CalculateResidualNorm();
    
    /**
     *  Take one Newton step, ie assemble the jacobian matrix and residual, solve for
     *  the update, and determine best damping value.
     */
    void TakeNewtonStep();
    
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
    
    
    /**
     *  Output current deformed position to file (or the undeformed mesh, if the 
     *  second parameter is set to false)
     *  @counter A number to suffix the file. The output file will be 
     *   <out_dir>/finiteelas_solution_<counter.[nodes/elem/undefnodes/undefelem]
     *  @writeDeformed whether to write the deformed position or the undeformed
     *   position, defaults to deformed
     */
    void WriteOutput(unsigned counter, bool writeDeformed=true);
    
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
     */
    virtual void Solve();
    
    /**
     *  Get the deformed position. rGetDeformedPosition()[i][j] is the x_i value at node j
     */
    std::vector<Vector<double> >& rGetDeformedPosition();
    /**
     *  Get the undeformed position. rGetUndeformedPosition()[i][j] is the X_i value at node j
     *  Obviously this data is accessible from the mesh as well, this method is more useful
     *  in some situations are. Note, this data structure is not set up unless this method 
     *  is called.
     */
    std::vector<Vector<double> >& rGetUndeformedPosition();
       
    /**
     *  Get the number of newton iterations that had been required to solve the problem
     */
    unsigned GetNumNewtonIterations();
    
    /**
     *  Verify that the analytic jacobian is the same as the numerical jacobian
     *  
     *  @param tol. The tolerance with which to compare the absolute component-wise
     *  difference between the two jacobians. 
     */
    void CompareJacobians(double tol=1e-8);
};




#endif /*FINITEELASTICITYASSEMBLER_HPP_*/
