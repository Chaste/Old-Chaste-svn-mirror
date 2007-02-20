#ifndef FINITEELASTICITYASSEMBLER_HPP_
#define FINITEELASTICITYASSEMBLER_HPP_


// speed notes:
//  - don't use Tensor<4,DIM>, use double[][][][] - factor of about 2 difference
//  - compile using OPTIMISATION!!!!!! - factor of about 200 difference when 
//     assembling system!!
//  - make reused variables static rather then repeatedly creating them? or member 
//     variables. static may cause problems if a 2d sim is run then a 3d sim is run? 
//     or probably not.


// todos: useful output, bcc object, read in dirichlet/neumann?, heterogeneity, nondim?


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
private:
    Triangulation<DIM>*   mpMesh;

    // an FE_Q object seems to be equivalent to our basis functions
    // It is templated over dimension, with the order of the bases taken 
    // in the contructor
    /*FE_Q<DIM>            mFe;*/            // note that this must be defined before mDofHandler!

    // an Fe system seems to be, loosely, a set of fe_q objects
    FESystem<DIM>        mFeSystem;          // note that this must be defined before mDofHandler!
    DoFHandler<DIM>      mDofHandler;

    SparsityPattern      mSparsityPattern;

    SparseMatrix<double> mJacobianMatrix;
    Vector<double>       mCurrentSolution;
    Vector<double>       mResidual;
    
    std::string          mOutputDirectoryFullPath;

    double dTdE[DIM][DIM][DIM][DIM]; 

    AbstractIncompressibleMaterialLaw<DIM>*  mpMaterialLaw;
    Vector<double>       mBodyForce;
    double               mDensity;

    const unsigned unsigned   PRESSURE_COMPONENT_INDEX; // just set to be DIM, ie if DIM==2, 
                                                   // the spatial indices are 0 and 1,
                                                   // the pressure index is 2


    std::map<unsigned,double> mBoundaryValues;


    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter, 
                           Vector<double>&                                 elementRhs,
                           FullMatrix<double>&                             elementMatrix,
                           bool                                            assembleResidual,
                           bool                                            assembleJacobian);

    void AssembleSystem(bool assembleResidual, bool assembleJacobian);
    void ApplyDirichletBoundaryConditions(bool assembleResidual, bool assembleJacobian);
    
    void OutputResults(unsigned newtonIteration);
    double CalculateResidualNorm();
    
    
public:
    FiniteElasticityAssembler(Triangulation<DIM>* pMesh,
                              AbstractIncompressibleMaterialLaw<DIM>*  pMaterialLaw,
                              Vector<double> bodyForce,
                              double density,
                              std::string outputDirectory,
                              unsigned orderOfBasesForPosition=2,
                              unsigned orderOfBasesForPressure=1);
    ~FiniteElasticityAssembler();
    
    
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

    void Solve();
    
    Vector<double>& GetSolutionVector();
    DoFHandler<DIM>& GetDofHandler();
    
};




#endif /*FINITEELASTICITYASSEMBLER_HPP_*/
