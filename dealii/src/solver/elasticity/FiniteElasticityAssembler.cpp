#ifndef FINITEELASTICITYASSEMBLER_CPP_
#define FINITEELASTICITYASSEMBLER_CPP_

#include "FiniteElasticityAssembler.hpp"

#include <dofs/dof_tools.h>

template<int DIM>
FiniteElasticityAssembler<DIM>::FiniteElasticityAssembler(Triangulation<DIM>* pMesh,
                                                          AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                                          Vector<double> bodyForce,
                                                          double density,
                                                          std::string outputDirectory,
                                                          unsigned orderOfBasesForPosition,
                                                          unsigned orderOfBasesForPressure
                                                          )  :
            // DIM bases for position, 1 for pressure
            mFeSystem(FE_Q<DIM>(orderOfBasesForPosition), DIM, FE_Q<DIM>(1), orderOfBasesForPressure),
            mDofHandler(*pMesh),                            // associate the mesh with the dof handler
            mBodyForce(bodyForce),
            mDensity(density),
            PRESSURE_COMPONENT_INDEX(DIM) // ie if DIM=2, the space indices are 0 and 1, pressure index is 2
{
    assert(pMesh!=NULL); // probably will fail in the mDofHandler(*pMesh) line above before here if pMesh==NULL
    mpMesh = pMesh;
    
    assert(pMaterialLaw != NULL);
    mpMaterialLaw = pMaterialLaw;

    assert(bodyForce.size()==DIM);
    assert(density > 0.0);
    
    OutputFileHandler output_file_handler(outputDirectory);
    mOutputDirectoryFullPath = output_file_handler.GetTestOutputDirectory(outputDirectory);
    
    // check the mesh has a region on the surface which has been set to 
    // be the fixed boudary.

    bool found_fixed_boundary = false;
    
    //loop over surface elements and set indicator as dirichlet or neumman
    typename Triangulation<DIM>::cell_iterator element_iter = mpMesh->begin();
    while(element_iter!=mpMesh->end())
    {
        for(unsigned face_index=0; face_index<GeometryInfo<DIM>::faces_per_cell; face_index++)
        {
            // note: the boundary_indicator is set to be 255 for internal faces, at_boundary()
            // essentially checks whether face->boundary_indicator()==255.
            if(element_iter->face(face_index)->at_boundary()) 
            {
                if(element_iter->face(face_index)->boundary_indicator()==FIXED_BOUNDARY)
                {
                    found_fixed_boundary = true;
                }
            }
       }
        element_iter++;
    }
    
    if(!found_fixed_boundary)
    {
        EXCEPTION("No fixed surface found. (no surface elements in the mesh have had their boundary indicator set to be FIXED_BOUNDARY");
    }         

    // distribute dofs
    /* mDofHandler.distribute_dofs(mFe);*/
    mDofHandler.distribute_dofs(mFeSystem);
    
    // form sparsity pattern
    mSparsityPattern.reinit(mDofHandler.n_dofs(), 
                            mDofHandler.n_dofs(),
                            mDofHandler.max_couplings_between_dofs());

    DoFTools::make_sparsity_pattern(mDofHandler, mSparsityPattern);
    mSparsityPattern.compress();
    
    // initialise vectors and matrices
    mJacobianMatrix.reinit(mSparsityPattern);
    mCurrentSolution.reinit(mDofHandler.n_dofs());
    mResidual.reinit(mDofHandler.n_dofs());
    
    std::cerr << "Number of active cells: " << mpMesh->n_active_cells() << std::endl;
    std::cerr << "Total number of cells: "  << mpMesh->n_cells() << std::endl;
    std::cerr << "Number of degrees of freedom: " << mDofHandler.n_dofs() << std::endl;    
    
    std::vector<bool> component_mask(DIM+1);

    for(int i=0; i<DIM; i++)
    {
        component_mask[i] = true;
    }
    component_mask[DIM] = false;

    VectorTools::interpolate_boundary_values(mDofHandler,
                                             FIXED_BOUNDARY,
                                             ZeroFunction<DIM>(DIM+1),  // note the "+1" here! - number of components
                                             mBoundaryValues,
                                             component_mask);

    std::map<unsigned,double>::iterator iter = mBoundaryValues.begin();


// random inputing code
//    GridIn<DIM> grid_in;
//    grid_in.attach_triangulation(mMesh);
//    std::ifstream input(meshFile.c_str());
//    grid_in.read_ucd(input);
    
// random outputing code
//    std::ofstream out("mesh.inp");
//    GridOut grid_out;
//    grid_out.write_ucd(mMesh, out);

}

template<int DIM>
FiniteElasticityAssembler<DIM>::~FiniteElasticityAssembler()
{
}


//template<int DIM>
//void FiniteElasticityAssembler<DIM>::SetDisplacementBoundaryConditions(std::vector<unsigned> node,
//                                                                       std::vector<unsigned> coordinate,
//                                                                       std::vector<double> value)
//{
//    mBoundaryValues.clear();
//
//    assert(node.size()==coordinate.size());
//    assert(node.size()==value.size());
//    
//    unsigned num_bcs = node.size();
//    
//    for(unsigned i=0; i<num_bcs; i++) 
//    {
//        assert(coordinate[i] < DIM);
//    }
//    
//    DofVertexIterator<DIM> vertex_iter(mpMesh, &mDofHandler);
//    
//    while(!vertex_iter.ReachedEnd())
//    {
//        unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
//        
//        for(unsigned i=0; i<num_bcs; i++)
//        {
//            if( node[i]==vertex_index )
//            {
//                unsigned dof = vertex_iter.GetDof( coordinate[i] );
//                mBoundaryValues[dof] = value[i];
//            }
//        }
//        
//        vertex_iter.Next();
//    }
//}    
//
//template<int DIM>
//void FiniteElasticityAssembler<DIM>::SetFixedNodes(std::vector<unsigned> nodes)                                        
//{
//    mBoundaryValues.clear();
//    
//    DofVertexIterator<DIM> vertex_iter(mpMesh, &mDofHandler);
//    
//    while(!vertex_iter.ReachedEnd())
//    {
//        unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
//        
//        for(unsigned i=0; i<nodes.size(); i++)
//        {
//            if( nodes[i]==vertex_index )
//            {
//                for(unsigned j=0; j<DIM; j++)
//                {
//                    unsigned dof = vertex_iter.GetDof(j);
//                    mBoundaryValues[dof] = 0.0;
//                }
//            }
//        }
//        
//        vertex_iter.Next();
//    }
//    
//    unsigned quads[5] = {41,53,56,65,83};
//    
//    for(unsigned i=0; i<5; i++)
//    {
//        mBoundaryValues[quads[i]] = 0.0;
//        mBoundaryValues[quads[i]+1] = 0.0;
//        mBoundaryValues[quads[i]+2] = 0.0;
//    }
//    
//    
//    
//    std::map<unsigned,double>::iterator iter = mBoundaryValues.begin();
//
//    while(iter!=mBoundaryValues.end())
//    {
//        std::cout << iter->first << " " << iter->second << "\n";
//        iter++;
//    }    
//}    

template<int DIM>
void FiniteElasticityAssembler<DIM>::SetBoundaryValues(std::map<unsigned,double> boundaryValues)
{
    assert(!boundaryValues.empty());
    
    mBoundaryValues.clear();
    std::map<unsigned,double>::iterator iter = boundaryValues.begin();
    while(iter!=boundaryValues.end())
    {
        mBoundaryValues[ iter->first ] = iter->second;
        iter++;
    }   

    assert(!mBoundaryValues.empty());
}



//////////////////////////////////////////////////////////////////////////////////////////
// AssembleOnElement
//////////////////////////////////////////////////////////////////////////////////////////
template<int DIM>
void FiniteElasticityAssembler<DIM>::AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter, 
                                                       Vector<double>&       elementRhs,
                                                       FullMatrix<double>&   elementMatrix,
                                                       bool                  assembleResidual,
                                                       bool                  assembleJacobian
                                                       )
{    
    static QGauss<DIM>   quadrature_formula(3);
    

    static QGauss<DIM-1> face_quadrature_formula(3);
    
    const unsigned n_q_points    = quadrature_formula.n_quadrature_points;
    const unsigned n_face_q_points = face_quadrature_formula.n_quadrature_points;
    
    
    // would want this to be static too (slight speed up), but causes errors
    // in debug mode (upon destruction of the class, in 2d, or something)
    FEValues<DIM> fe_values(mFeSystem, quadrature_formula, 
                            UpdateFlags(update_values    |
                                        update_gradients |
                                        update_q_points  |     // needed for interpolating u and u' on the quad point
                                        update_JxW_values));

    // would want this to be static too (slight speed up), but causes errors
    // in debug mode (upon destruction of the class, in 2d, or something)
    FEFaceValues<DIM> fe_face_values(mFeSystem, face_quadrature_formula, 
                                     UpdateFlags(update_values         |
                                                 update_q_points       |
                                                 update_normal_vectors |
                                                 update_JxW_values));

    
    const unsigned dofs_per_element = mFeSystem.dofs_per_cell;


    static std::vector< unsigned >                    local_dof_indices(dofs_per_element);
    static std::vector< Vector<double> >                  local_solution_values(n_q_points);
    static std::vector< std::vector< Tensor<1,DIM> > >    local_solution_gradients(n_q_points);

    static Tensor<2,DIM> identity;         // how do you do this properly??

    static bool first = true;


    if(first)
    {
        for(unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            local_solution_values[q_point].reinit(DIM+1);
            local_solution_gradients[q_point].resize(DIM+1);
        }
        for(unsigned i=0; i<DIM; i++)
        {
            for(unsigned j=0; j<DIM; j++)
            {
                identity[i][j] = i==j ? 1.0 : 0.0;
            }
        }
    }
    

    elementMatrix = 0;
    elementRhs = 0;
      
    elementIter->get_dof_indices(local_dof_indices);

    fe_values.reinit(elementIter); // compute fe values for this element
    fe_values.get_function_values(mCurrentSolution, local_solution_values);
    fe_values.get_function_grads(mCurrentSolution, local_solution_gradients);        


    for(unsigned q_point=0; q_point<n_q_points; q_point++)
    {       
        const std::vector< Tensor<1,DIM> >& grad_u_p = local_solution_gradients[q_point];
        
        double p = local_solution_values[q_point](PRESSURE_COMPONENT_INDEX);
        
        static Tensor<2,DIM> F;
        static Tensor<2,DIM> C;
        static Tensor<2,DIM> inv_C;
        static Tensor<2,DIM> inv_F;
        static SymmetricTensor<2,DIM> T;

        for(unsigned i=0; i<DIM; i++)
        {
            for(unsigned j=0; j<DIM; j++)
            {
                F[i][j] = identity[i][j] + grad_u_p[i][j];
            }
        }
        
        C = transpose(F) * F;
        inv_C = invert(C);
        inv_F = invert(F);

        double detF = determinant(F);

        static SymmetricTensor<2,DIM> T2;
        mpMaterialLaw->ComputeStressAndStressDerivative(C,inv_C,p,T,dTdE,assembleJacobian);

        for(unsigned i=0; i<dofs_per_element; i++)
        {
            const unsigned component_i = mFeSystem.system_to_component_index(i).first;
            
            if(assembleJacobian)
            {
                for(unsigned j=0; j<dofs_per_element; j++)
                {
                    const unsigned component_j = mFeSystem.system_to_component_index(j).first;
        
                    if((component_i<PRESSURE_COMPONENT_INDEX) &&(component_j<PRESSURE_COMPONENT_INDEX) )
                    {
                        for(unsigned M=0; M<DIM; M++)
                        {
                            for(unsigned N=0; N<DIM; N++)
                            {
                                elementMatrix(i,j) +=   T[M][N] 
                                                      * fe_values.shape_grad(j,q_point)[M]
                                                      * fe_values.shape_grad(i,q_point)[N]
                                                      * identity[component_i][component_j]
                                                      * fe_values.JxW(q_point);
    
                                for(unsigned P=0; P<DIM; P++)
                                {
                                    for(unsigned Q=0; Q<DIM; Q++)
                                    {
                                        elementMatrix(i,j) +=   0.5 
                                                              * dTdE[M][N][P][Q]
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
                    }
                    else if((component_i<PRESSURE_COMPONENT_INDEX) &&(component_j==PRESSURE_COMPONENT_INDEX) )
                    {
                        for(unsigned M=0; M<DIM; M++)
                        {
                            for(unsigned N=0; N<DIM; N++)
                            {
                                elementMatrix(i,j) +=  - F[component_i][M]
                                                       * inv_C[M][N]
                                                       * fe_values.shape_grad(i,q_point)[N] 
                                                       * fe_values.shape_value(j,q_point)
                                                       * fe_values.JxW(q_point);
                            }
                        }
                    }
                    else if((component_i==PRESSURE_COMPONENT_INDEX) &&(component_j<PRESSURE_COMPONENT_INDEX) )
                    {
                        for(unsigned M=0; M<DIM; M++)
                        {
                            elementMatrix(i,j) +=    fe_values.shape_value(i,q_point)
                                                   * detF
                                                   * inv_F[M][component_j]
                                                   * fe_values.shape_grad(j,q_point)[M] 
                                                   * fe_values.JxW(q_point);
                        }
                    }
                    else
                    {
                        // do nothing, ie elementMatrix(i,j)  +=  0 * fe_values.JxW(q_point);;
                    }
                }                
            }
            
            if(assembleResidual)
            {
                if(component_i<PRESSURE_COMPONENT_INDEX)
                {
                    elementRhs(i) += - mDensity * mBodyForce(component_i)
                                     * fe_values.shape_value(i,q_point)
                                     * fe_values.JxW(q_point);
                                
                    for(unsigned M=0; M<DIM; M++)
                    {
                        for(unsigned N=0; N<DIM; N++)
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
    }


    ////////////////////////////
    // loop over faces
    ////////////////////////////
    if(assembleResidual)
    {
        for (unsigned face_index=0; face_index<GeometryInfo<DIM>::faces_per_cell; face_index++)
        {
            if (elementIter->face(face_index)->boundary_indicator() == NEUMANN_BOUNDARY)
            {
                fe_face_values.reinit(elementIter, face_index);
  
                for (unsigned q_point=0; q_point<n_face_q_points; q_point++)
                {
                    Vector<double> neumann_traction(DIM); // zeros
    
                    neumann_traction(1)=0.0;
    
                    for (unsigned i=0; i<dofs_per_element; i++)
                    {
                        const unsigned component_i = mFeSystem.system_to_component_index(i).first;

                        if(component_i < PRESSURE_COMPONENT_INDEX)
                        {
                            elementRhs(i) +=   neumann_traction(component_i)
                                             * fe_face_values.shape_value(i,q_point)
                                             * fe_face_values.JxW(q_point);
                        }
                    }
                }
            }
        }
    }

    first = false;
}



//////////////////////////////////////////////////////////////////////////////////////////
// AssembleSystem
//////////////////////////////////////////////////////////////////////////////////////////
template<int DIM>
void FiniteElasticityAssembler<DIM>::AssembleSystem(bool assembleResidual, 
                                                    bool assembleJacobian)
{
    const unsigned       dofs_per_element = mFeSystem.dofs_per_cell;
  
    FullMatrix<double>   element_matrix(dofs_per_element, dofs_per_element);
    Vector<double>       element_rhs(dofs_per_element);

    // the dofs associated with the nodes of an element
    std::vector<unsigned> local_dof_indices(dofs_per_element);

    typename DoFHandler<DIM>::active_cell_iterator  element_iter = mDofHandler.begin_active();
    
    
    if(assembleResidual)
    {
        mResidual = 0;
    }
    
    if(assembleJacobian)
    {
        mJacobianMatrix = 0;
    }

    unsigned elem_counter = 0;
    
    while(element_iter!=mDofHandler.end())   // huh? mDof.end() returns an element iterator?
    {
        // zero the small matrix and vector
        element_matrix = 0;
        element_rhs = 0;
      
        element_iter->get_dof_indices(local_dof_indices);

        AssembleOnElement(element_iter,
                          element_rhs, 
                          element_matrix,
                          assembleResidual,                     
                          assembleJacobian);                    

        if(assembleJacobian)
        {
            std::cout << elem_counter++ << " " << std::flush;
        }

        for(unsigned i=0; i<dofs_per_element; i++)
        {
            if(assembleJacobian)
            {
                for(unsigned j=0; j<dofs_per_element; j++)
                {
                    mJacobianMatrix.add(local_dof_indices[i],
                                        local_dof_indices[j], 
                                        element_matrix(i,j));
                }
            }
            
            if(assembleResidual)
            {
                mResidual(local_dof_indices[i]) += element_rhs(i);
            }
        }

        element_iter++;
    }
    if(assembleJacobian) { std::cout << "\n"; }
    
    ApplyDirichletBoundaryConditions(assembleResidual,assembleJacobian);
}


template<int DIM>
void FiniteElasticityAssembler<DIM>::ApplyDirichletBoundaryConditions(bool assembleResidual,
                                                                      bool assembleJacobian)
{
    // The boundary conditions on the NONLINEAR SYSTEM are x=boundary_values 
    // on the boundary nodes. However:
    // The boundary conditions on the LINEAR SYSTEM  Ju=f, where J is the 
    // u the negative update vector and f is the residual is 
    // u=current_soln-boundary_values on the boundary nodes
    std::map<unsigned,double>  applied_boundary_values;

    std::map<unsigned,double>::iterator iter = mBoundaryValues.begin();
    while(iter!=mBoundaryValues.end())
    {
        unsigned dof = iter->first;
        double value = iter->second;
   
        applied_boundary_values[dof] = mCurrentSolution(dof)-mBoundaryValues[dof];
        iter++;
    }
    
    // don't have access to u (the solution of the linear system) at the moment,
    // so pass in a dummy vector 
    Vector<double> dummy(mResidual.size());

    MatrixTools::apply_boundary_values(applied_boundary_values,
                                       mJacobianMatrix,
                                       dummy,
                                       mResidual);
}


template<int DIM>
void FiniteElasticityAssembler<DIM>::OutputResults(unsigned newtonIteration)
{
    std::stringstream ss;
    ss << mOutputDirectoryFullPath << "/finiteelas_solution_" << newtonIteration << ".gmv";
    std::string filename = ss.str();
    std::ofstream output(filename.c_str());
  
    DataOut<DIM> data_out;
    data_out.attach_dof_handler(mDofHandler);
  
    std::vector<std::string> solution_names;

    solution_names.push_back("x_displacement");
    if(DIM>1)
    {
        solution_names.push_back("y_displacement");
    }
    if(DIM>2)
    {
        solution_names.push_back("z_displacement");
    }
        
    solution_names.push_back("pressure");

    data_out.add_data_vector(mCurrentSolution, solution_names);
    data_out.build_patches();
    data_out.write_gmv(output);
}



template<int DIM>
void FiniteElasticityAssembler<DIM>::Solve()
{
    OutputResults(0);
    
    // compute residual
    AssembleSystem(true, false);
    double norm_resid = CalculateResidualNorm();
    std::cout << "\nNorm of residual is " << norm_resid << "\n";
    
    unsigned counter = 1;
    
    // use the larger of the tolerances formed from the absolute or 
    // relative possibilities
    double tol = NEWTON_ABS_TOL;
    if( tol < NEWTON_REL_TOL*norm_resid )
    {
        tol = NEWTON_REL_TOL*norm_resid;
    }
    std::cout << "Solving with tolerance " << tol << "\n";                                            
    
     
    while(norm_resid > tol && counter < 10)
    {
        std::cout <<  "\n-------------------\n"
                  <<   "Newton iteration " << counter
                  << ":\n-------------------\n";

        // compute Jacobian
        AssembleSystem(false, true);

        // solve the linear system
        SolverControl  solver_control(20000, 1e-6, false, false);
        PrimitiveVectorMemory<> vector_memory;

        Vector<double> update;
        update.reinit(mDofHandler.n_dofs());
    
        SolverGMRES<>::AdditionalData gmres_additional_data(200);
        SolverGMRES<>  gmres(solver_control, vector_memory, gmres_additional_data);    
        gmres.solve(mJacobianMatrix, update, mResidual, PreconditionIdentity());
        

        // save the old current solution
        Vector<double> old_solution = mCurrentSolution;

        double best_norm_resid = norm_resid;
        double best_damping_value = 0.0;
        
        std::vector<double> damping_values;
        damping_values.push_back(0.0);
        damping_values.push_back(0.05);
        for(unsigned i=1; i<=10; i++)
        {
            damping_values.push_back((double)i/10.0);
        }

        for(unsigned i=0; i<damping_values.size(); i++)
        {
            mCurrentSolution.equ(1.0, old_solution, -damping_values[i], update);
            
            // compute residual
            AssembleSystem(true, false);
            norm_resid = CalculateResidualNorm();
            
            std::cout << "\tTesting s = " << damping_values[i] << ", |f| = " << norm_resid << "\n";
            
            if(norm_resid < best_norm_resid)
            {
                best_norm_resid = norm_resid;
                best_damping_value = damping_values[i];
            }
        }
        
        if(best_damping_value == 0.0)
        {
            std::cout << "\nResidual does not decrease in newton direction, quitting\n";
            assert(0);
        }
        else
        {
            std::cout << "\tBest s = " << best_damping_value << "\n";
        }
        
        // implement best update and recalculate residual
        mCurrentSolution.equ(1.0, old_solution, -best_damping_value, update);
        AssembleSystem(true, false);
        norm_resid = CalculateResidualNorm();

        std::cout << "Norm of residual is " << norm_resid << "\n";

        OutputResults(counter);

        counter++;
    }
    
    if(norm_resid > tol)
    {
        EXCEPTION("Failed to converge");
    }
    

//    std::cout << "Node    X    Y    x    y\n";
//
//    std::vector<bool> vertex_touched(mpMesh->n_vertices(),false);
//
//    typename DoFHandler<DIM>::active_cell_iterator cell = mDofHandler.begin_active();
//    while(cell != mDofHandler.end())
//    {
//        for (unsigned v=0; v<GeometryInfo<DIM>::vertices_per_cell; v++)
//        {
//            if(vertex_touched[cell->vertex_index(v)] == false)
//            {
//                vertex_touched[cell->vertex_index(v)] = true;
//
//                std::cout << cell->vertex_index(v) << "  "
//                          << cell->vertex(v)(0) << "   "  
//                          << cell->vertex(v)(1) << "   ";  
//
//                unsigned dof = cell->vertex_dof_index(v,0);
//                std::cout << mCurrentSolution(dof) << "   ";
//                                
//                dof = cell->vertex_dof_index(v,1);
//                std::cout << mCurrentSolution(dof) << "\n";
//            }
//        }
//        cell++;
//    }
}


template<int DIM>
double FiniteElasticityAssembler<DIM>::CalculateResidualNorm()
{
    return mResidual.norm_sqr()/mDofHandler.n_dofs();
}


template<int DIM>
Vector<double>& FiniteElasticityAssembler<DIM>::GetSolutionVector()
{
    return mCurrentSolution;
}


template<int DIM>
DoFHandler<DIM>& FiniteElasticityAssembler<DIM>::GetDofHandler()
{
    return mDofHandler;
}


#endif // FINITEELASTICITYASSEMBLER_CPP_


