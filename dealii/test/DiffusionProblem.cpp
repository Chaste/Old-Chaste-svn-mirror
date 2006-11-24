#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>

#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>

#include <dofs/dof_tools.h>

#include <fe/fe_values.h> // *** PROBLEM with boost
#include <base/quadrature_lib.h>

#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <sstream>


class DiffusionProblem 
{
public:
    DiffusionProblem();
    void Run();
      
private:
    void MakeGridAndDofs();
    void AssembleSystem();
    void Solve();
    void OutputResults(int timeStepCounter) const;
  
    Triangulation<2>     mMesh;
    FE_Q<2>              mFe;
    DoFHandler<2>        mDofHandler;
  
//    SparsityPattern      mSparsityPattern;
//    SparseMatrix<double> mSystemMatrix;
//  
    Vector<double>       mSolution;
//    Vector<double>       mSystemRhs;
      
    bool                 mMatrixIsAssembled;
    double               mDt;
    double               mDiffusionCoeff;
};
  
  
DiffusionProblem::DiffusionProblem() :
                  mFe(1),
                  mDofHandler(mMesh)
{
    mMatrixIsAssembled = false;
    mDt = 0.01;
    mDiffusionCoeff = 0.1;
}
  
  
void DiffusionProblem::MakeGridAndDofs()
{
    GridGenerator::hyper_cube(mMesh, -1, 1);
    mMesh.refine_global(4);
    std::cout << "Number of active cells: "
              << mMesh.n_active_cells()
              << std::endl;
    std::cout << "Total number of cells: "
              << mMesh.n_cells()
              << std::endl;
    
    mDofHandler.distribute_dofs(mFe);
  
    std::cout << "Number of degrees of freedom: "
              << mDofHandler.n_dofs()
              << std::endl;
  
//    mSparsityPattern.reinit(mDofHandler.n_dofs(),
//                            mDofHandler.n_dofs(),
//                            mDofHandler.max_couplings_between_dofs());
//    DoFTools::make_sparsity_pattern(mDofHandler, mSparsityPattern);
//    mSparsityPattern.compress();
//  
//    mSystemMatrix.reinit(mSparsityPattern);
//  
    mSolution.reinit(mDofHandler.n_dofs());
//    mSystemRhs.reinit(mDofHandler.n_dofs());
//    
//    // initialise solution
//    for(unsigned dof=0; dof<mDofHandler.n_dofs(); dof++)
//    {
//        mSolution(dof)=0;
//    }
}
  
  
void DiffusionProblem::AssembleSystem() 
{
//    QGauss<2>  quadrature_formula(2);
//    FEValues<2> fe_values(mFe, quadrature_formula, 
//                          UpdateFlags(update_values    |
//                                      update_gradients |
//                                      update_JxW_values));
//  
//    const unsigned int   dofs_per_cell = mFe.dofs_per_cell;
//    const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;
//  
//    FullMatrix<double>   cell_matrix(dofs_per_cell, dofs_per_cell);
//    Vector<double>       cell_rhs(dofs_per_cell);
//  
//    std::vector<unsigned int> local_dof_indices(dofs_per_cell);
//  
//    DoFHandler<2>::active_cell_iterator cell = mDofHandler.begin_active(),
//                                        endc = mDofHandler.end();
//                                        
//    for(unsigned dof=0; dof<mDofHandler.n_dofs(); dof++)
//    {
//        mSystemRhs(dof)=0.0;
//    } 
//                                        
//                                        
//    for(; cell!=endc; cell++)
//    {
//        fe_values.reinit(cell);
//  
//        cell_matrix = 0;
//        cell_rhs = 0;
//  
//        cell->get_dof_indices(local_dof_indices);
//  
//        if(!mMatrixIsAssembled)
//        {
//            for(unsigned int i=0; i<dofs_per_cell; i++)
//            {
//                for(unsigned int j=0; j<dofs_per_cell; j++)
//                {
//                    for(unsigned int q_point=0; q_point<n_q_points; q_point++)
//                    {
//                        cell_matrix(i,j) +=(fe_values.shape_grad(i, q_point) *
//                                            fe_values.shape_grad(j, q_point) *
//                                            fe_values.JxW(q_point))*mDt*mDiffusionCoeff;
// 
//                        cell_matrix(i,j) += fe_values.shape_value(i, q_point) *
//                                            fe_values.shape_value(j, q_point) *
//                                            fe_values.JxW(q_point);
//                    }
//                }
//            }
//        }                           
//  
//        for(unsigned int i=0; i<dofs_per_cell; i++)
//        {
//            for(unsigned int q_point=0; q_point<n_q_points; q_point++)
//            {
//                cell_rhs(i) += (fe_values.shape_value(i, q_point) *
//                                1 * // source term
//                                fe_values.JxW(q_point)) *mDt;
//             
//                double U=0;
//                for(unsigned int j=0; j<dofs_per_cell; j++)
//                {
//                    U += mSolution(local_dof_indices[j])*fe_values.shape_value(j,q_point); 
//                }
//             
//                cell_rhs(i) += fe_values.shape_value(i, q_point) *
//                               fe_values.JxW(q_point) *
//                               U;
//            }
//        }
//  
//        if(!mMatrixIsAssembled)
//        {
//            for(unsigned int i=0; i<dofs_per_cell; i++)
//            {
//                for(unsigned int j=0; j<dofs_per_cell; j++)
//                {
//                    mSystemMatrix.add(local_dof_indices[i],
//                                      local_dof_indices[j],
//                                      cell_matrix(i,j));
//                }
//            }
//        }
//  
//        for(unsigned int i=0; i<dofs_per_cell; i++)
//        {
//            mSystemRhs(local_dof_indices[i]) += cell_rhs(i);
//        }
//    }
//  
//  
//    std::map<unsigned int,double> boundary_values;
//    VectorTools::interpolate_boundary_values(mDofHandler,
//                                             0,
//                                             ZeroFunction<2>(),
//                                             boundary_values);
//    
//    MatrixTools::apply_boundary_values(boundary_values,
//                                        mSystemMatrix,
//                                        mSolution,
//                                        mSystemRhs);
//                                        
//    mMatrixIsAssembled = true;
}
  
  
void DiffusionProblem::Solve() 
{
//    SolverControl           solver_control(1000, 1e-12);
//    PrimitiveVectorMemory<> vector_memory;
//    SolverCG<>              cg(solver_control, vector_memory);
//  
//    cg.solve(mSystemMatrix, mSolution, mSystemRhs,
//             PreconditionIdentity());
}
  
  
void DiffusionProblem::OutputResults(int timeStepCounter) const
{
    DataOut<2> data_out;
    data_out.attach_dof_handler(mDofHandler);
    std::stringstream filename;
    filename << "solution_" << timeStepCounter << ".gpl";
    data_out.add_data_vector(mSolution, "solution");
    data_out.build_patches();

    std::ofstream output(filename.str().c_str());
    data_out.write_gnuplot(output);
}
  
  
void DiffusionProblem::Run() 
{
    MakeGridAndDofs();
    
    double t=0;
    int time_step_counter=0;
    OutputResults(time_step_counter);
  
    while(t < 2)
    {
//        AssembleSystem();
//        Solve();
        t += mDt;
//        time_step_counter++;
//        std::cout << "Time = " << t << std::endl;
//
//        OutputResults(time_step_counter);
    }
}
