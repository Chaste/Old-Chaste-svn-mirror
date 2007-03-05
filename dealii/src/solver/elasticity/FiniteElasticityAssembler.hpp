#ifndef FINITEELASTICITYASSEMBLER_HPP_
#define FINITEELASTICITYASSEMBLER_HPP_


// speed notes:
//  - don't use Tensor<4,DIM>, use double[][][][] - factor of about 2 difference
//  - compile using OPTIMISATION!!!!!! - factor of about 200 difference when 
//     assembling system!!
//  - make reused variables static rather then repeatedly creating them? or member 
//     variables. static may cause problems if a 2d sim is run then a 3d sim is run? 
//     or probably not.



// TODO: TEST AGAINST RESULTS FROM SOMEWHERE ELSE

// todo: possibly important: initial guess is the zero displacement solution,
//  which isn't p=0
// refactor the assembler stuff into a dealii abstract assembler class?
// refactor out the newton solver?
// change quad rule if linears
// tobefixed: numerical jacobian: works but the boundary condition is wrong for some
// reason

// other todos: doxygen,chaste style output, nonzero neumann, heterogeneity, nondim. 


#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>

#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>

#include <dofs/dof_tools.h>
 
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>

#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>

#include <lac/vector_memory.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <sstream>
#include <base/symmetric_tensor.h>

// for dealing with hanging nodes..
#include <dofs/dof_constraints.h>

#define FIXED_BOUNDARY 10
#define NEUMANN_BOUNDARY 11
#define DIRICHLET_BOUNDARY 12

#define NEWTON_ABS_TOL 1e-11
#define NEWTON_REL_TOL 1e-7

#include "AbstractIncompressibleMaterialLaw.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
#include "DofVertexIterator.hpp"

template<int DIM>
class FiniteElasticityAssembler
{
protected:
    Triangulation<DIM>*   mpMesh;

    // an FE_Q object seems to be equivalent to our basis functions
    // It is templated over dimension, with the order of the bases taken 
    // in the contructor
    // note that this must be defined before mDofHandler!
    /*FE_Q<DIM>            mFe;*/            

    // an Fe system seems to be, loosely, a set of fe_q objects
    FESystem<DIM>        mFeSystem;          // note that this must be defined before mDofHandler!
    DoFHandler<DIM>      mDofHandler;

    ConstraintMatrix     mHangingNodeConstraints;

    SparsityPattern      mSparsityPattern;

    SparseMatrix<double> mJacobianMatrix;
    Vector<double>       mCurrentSolution;
    Vector<double>       mResidual;
    
    std::string          mOutputDirectoryFullPath;
    bool                 mWriteOutput;

    double dTdE[DIM][DIM][DIM][DIM]; 


    bool mHeterogeneous;
    std::vector<AbstractIncompressibleMaterialLaw<DIM>*>  mMaterialLaws;
    std::vector<int> mMaterialIdToMaterialLawIndexMap;
    
    Vector<double>       mBodyForce;
    double               mDensity;

    const unsigned       PRESSURE_COMPONENT_INDEX; // just set to be DIM, ie if DIM==2, 
                                                   // the spatial indices are 0 and 1,
                                                   // the pressure index is 2


    std::map<unsigned,double> mBoundaryValues;
    SparseMatrix<double> mNumericalJacobianMatrix;

    virtual void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter, 
                                   Vector<double>&                                 elementRhs,
                                   FullMatrix<double>&                             elementMatrix,
                                   bool                                            assembleResidual,
                                   bool                                            assembleJacobian);

    void AssembleSystem(bool assembleResidual, bool assembleJacobian);
    void ApplyDirichletBoundaryConditions(bool assembleResidual, bool assembleJacobian);
    
    void OutputResults(unsigned counter);
    double CalculateResidualNorm();
    
    void TakeNewtonStep();
    
    unsigned GetMaterialLawIndexFromMaterialId(unsigned materialId);
    
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

    
    
    //// this type of function doesn't really work
    //void SetDisplacementBoundaryConditions(std::vector<unsigned> nodes, 
    //                                       std::vector<unsigned> coordinates, 
    //                                       std::vector<double>   values);

    //// this type of function doesn't really work
    //void SetFixedNodes(std::vector<unsigned> nodes);
    
    // Setting boundary conditions is a hassle. Currently, assuming the 
    // default boundary conditions are not required, the user has to set 
    // up the dof->value map themselves and pass it in using this method. 
    // Note: call GetDofHandler() to get the dof handler first.
    void SetBoundaryValues(std::map<unsigned, double> boundary_values);


    void SetMaterialLawsForHeterogeneousProblem(std::vector<AbstractIncompressibleMaterialLaw<DIM>*> materialLaws,
                                                std::vector<unsigned> materialIds);

    virtual void Solve();
    
    Vector<double>& GetSolutionVector();
    DoFHandler<DIM>& GetDofHandler();
    Triangulation<DIM>* GetMesh();

    void ComputeNumericalJacobian();
    void CompareJacobians();    
};




#endif /*FINITEELASTICITYASSEMBLER_HPP_*/
