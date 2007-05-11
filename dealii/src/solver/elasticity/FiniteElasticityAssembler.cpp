#ifndef FINITEELASTICITYASSEMBLER_CPP_
#define FINITEELASTICITYASSEMBLER_CPP_

#include "FiniteElasticityAssembler.hpp"
#include "TriangulationVertexIterator.hpp"

#include <dofs/dof_tools.h>

template<unsigned DIM>
FiniteElasticityAssembler<DIM>::FiniteElasticityAssembler(Triangulation<DIM>* pMesh,
                                                          AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                                          Vector<double> bodyForce,
                                                          double density,
                                                          std::string outputDirectory,
                                                          unsigned degreeOfBasesForPosition,
                                                          unsigned degreeOfBasesForPressure
                                                         )  :
        AbstractDealiiAssembler<DIM>(pMesh),
        // DIM bases for position, 1 for pressure
        mFeSystem(FE_Q<DIM>(degreeOfBasesForPosition), DIM, FE_Q<DIM>(1), degreeOfBasesForPressure),
        mBodyForce(bodyForce),
        mDensity(density),
        PRESSURE_COMPONENT_INDEX(DIM) // ie if DIM=2, the space indices are 0 and 1, pressure index is 2
{
    // distribute dofs
    this->mDofHandler.distribute_dofs(mFeSystem);
    this->InitialiseMatricesVectorsAndConstraints();
    this->mDofsPerElement = mFeSystem.dofs_per_cell;
    
    
    if (pMaterialLaw != NULL)
    {
        mMaterialLaws.resize(1);
        mMaterialLaws[0] = pMaterialLaw;
        mHeterogeneous = false;
    }
    else
    {
        mHeterogeneous = true;
    }
    
    if (bodyForce.size()!=DIM)
    {
        EXCEPTION("Body force dimension does not match dimension of assembler");
    }
    if (density <= 0.0)
    {
        EXCEPTION("Density must be strictly positive");
    }
    
    if (outputDirectory!="")
    {
        mWriteOutput = true;
        OutputFileHandler output_file_handler(outputDirectory);
        mOutputDirectoryFullPath = output_file_handler.GetTestOutputDirectory(outputDirectory);
    }
    else
    {
        mWriteOutput = false;
    }
    
    // check the mesh has a region on the surface which has been set to
    // be the fixed boudary.
    
    bool found_fixed_boundary = false;
    
    //loop over surface elements and set indicator as dirichlet or neumman
    typename Triangulation<DIM>::cell_iterator element_iter = this->mpMesh->begin();
    while (element_iter!=this->mpMesh->end())
    {
        for (unsigned face_index=0; face_index<GeometryInfo<DIM>::faces_per_cell; face_index++)
        {
            // note: the boundary_indicator is set to be 255 for internal faces, at_boundary()
            // essentially checks whether face->boundary_indicator()==255.
            if (element_iter->face(face_index)->at_boundary())
            {
                if (element_iter->face(face_index)->boundary_indicator()==FIXED_BOUNDARY)
                {
                    found_fixed_boundary = true;
                }
            }
        }
        element_iter++;
    }
    
    if (!found_fixed_boundary)
    {
        EXCEPTION("No fixed surface found. (no surface elements in the mesh have had their boundary indicator set to be FIXED_BOUNDARY");
    }
    
    // initialise the deformed positions structure. It will
    // be a std::vector of x_posns, y_posns, z_posns, which
    // x_posns being a Vector<double> etc
    mDeformedPosition.resize(DIM);
    for (unsigned i=0; i<DIM; i++)
    {
        mDeformedPosition[i].reinit(this->mpMesh->n_vertices());
    }

    TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
    while(!vertex_iter.ReachedEnd())
    {
        unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
        Point<DIM> old_posn = vertex_iter.GetVertex();

        for(unsigned i=0; i<DIM; i++)
        {
            mDeformedPosition[i](vertex_index) =  old_posn(i);
        }
        
        vertex_iter.Next();
    }

//    std::cerr << "Number of active cells: " << this->mpMesh->n_active_cells() << std::endl;
//    std::cerr << "Total number of cells: "  << this->mpMesh->n_cells() << std::endl;
//    std::cerr << "Number of degrees of freedom: " << this->mDofHandler.n_dofs() << std::endl;

    std::vector<bool> component_mask(DIM+1);
    
    for (unsigned i=0; i<DIM; i++)
    {
        component_mask[i] = true;
    }
    component_mask[DIM] = false;
    
    VectorTools::interpolate_boundary_values(this->mDofHandler,
                                             FIXED_BOUNDARY,
                                             ZeroFunction<DIM>(DIM+1),  // note the "+1" here! - number of components
                                             mBoundaryValues,
                                             component_mask);
                                             
    mNumNewtonIterations = 0;
    
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

template<unsigned DIM>
FiniteElasticityAssembler<DIM>::~FiniteElasticityAssembler()
{}

template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::SetMaterialLawsForHeterogeneousProblem(std::vector<AbstractIncompressibleMaterialLaw<DIM>*> materialLaws,
        std::vector<unsigned> materialIds)
{
    // check sizes match
    if (materialLaws.size()!=materialIds.size())
    {
        EXCEPTION("materialLaws and materialIds must be the same size");
    }
    
    // set as heterogeneous (no checking that the sizes of these
    // vectors is greater than one at the moment
    mHeterogeneous = true;
    
    // copy the material laws
    assert(materialLaws.size()>0);
    mMaterialLaws = materialLaws;
    
    // check that every element in the mesh has material id which is in
    // the materialIds vector
    typename Triangulation<DIM>::cell_iterator element_iter = this->mpMesh->begin_active();
    while (element_iter!=this->mpMesh->end())
    {
        bool found_id = false;
        for (unsigned i=0; i<materialIds.size(); i++)
        {
            if (element_iter->material_id()==materialIds[i])
            {
                found_id = true;
            }
        }
        if (!found_id)
        {
            EXCEPTION("Found element in mesh whose material id does not correspond to any in materialIds");
        }
        element_iter++;
    }
    
    unsigned max_material_id = 0;
    for (unsigned i=0; i<materialIds.size(); i++)
    {
        if (materialIds[i]>max_material_id)
        {
            max_material_id = materialIds[i];
        }
    }
    
    mMaterialIdToMaterialLawIndexMap.resize(max_material_id+1);
    for (unsigned i=0; i<mMaterialIdToMaterialLawIndexMap.size(); i++)
    {
        mMaterialIdToMaterialLawIndexMap[i] = -1;
    }
    
    for (unsigned i=0; i<materialIds.size(); i++)
    {
        mMaterialIdToMaterialLawIndexMap[ materialIds[i] ]  = i;
    }
}




template<unsigned DIM>
unsigned FiniteElasticityAssembler<DIM>::GetMaterialLawIndexFromMaterialId(unsigned materialId)
{
    // something gone wrong in setting up this map if the
    // following assertion fails
    assert(mMaterialIdToMaterialLawIndexMap[materialId]!=-1);
    
    // another check
    assert(mMaterialIdToMaterialLawIndexMap[materialId]>=0);
    
    // convert int to unsigned, we know it is positive from above checks
    unsigned index = abs(mMaterialIdToMaterialLawIndexMap[materialId]);
    return index;
}


template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::FormInitialGuess()
{
    double zero_strain_pressure = mMaterialLaws[0]->GetZeroStrainPressure();
    
    std::vector<unsigned> local_dof_indices(this->mDofsPerElement);
    
    typename DoFHandler<DIM>::active_cell_iterator  element_iter = this->mDofHandler.begin_active();
    while (element_iter!=this->mDofHandler.end())
    {
        element_iter->get_dof_indices(local_dof_indices);
        
        for (unsigned i=0; i<this->mDofsPerElement; i++)
        {
            const unsigned component_i = mFeSystem.system_to_component_index(i).first;
            if (component_i == PRESSURE_COMPONENT_INDEX)
            {
                this->mCurrentSolution(local_dof_indices[i]) = zero_strain_pressure;
            }
        }
        element_iter++;
    }
}

//template<unsigned DIM>
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
//    DofVertexIterator<DIM> vertex_iter(this->mpMesh, &this->mDofHandler);
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
//template<unsigned DIM>
//void FiniteElasticityAssembler<DIM>::SetFixedNodes(std::vector<unsigned> nodes)
//{
//    mBoundaryValues.clear();
//
//    DofVertexIterator<DIM> vertex_iter(this->mpMesh, &this->mDofHandler);
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

template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::SetBoundaryValues(std::map<unsigned,double> boundaryValues)
{
    assert(!boundaryValues.empty());
    
    mBoundaryValues.clear();
    std::map<unsigned,double>::iterator iter = boundaryValues.begin();
    while (iter!=boundaryValues.end())
    {
        mBoundaryValues[ iter->first ] = iter->second;
        iter++;
    }
    
    assert(!mBoundaryValues.empty());
}



//////////////////////////////////////////////////////////////////////////////////////////
// AssembleOnElement
//////////////////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>
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
    
    
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law;
    if (!mHeterogeneous)
    {
        p_material_law = mMaterialLaws[0];
    }
    else
    {
        unsigned index = GetMaterialLawIndexFromMaterialId(elementIter->material_id());
        p_material_law = mMaterialLaws[index];
    }
    
    
    
    for (unsigned q_point=0; q_point<n_q_points; q_point++)
    {
        const std::vector< Tensor<1,DIM> >& grad_u_p = local_solution_gradients[q_point];
        
        double p = local_solution_values[q_point](PRESSURE_COMPONENT_INDEX);
        
        static Tensor<2,DIM> F;
        static Tensor<2,DIM> C;
        static Tensor<2,DIM> inv_C;
        static Tensor<2,DIM> inv_F;
        static SymmetricTensor<2,DIM> T;
        
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
        
        p_material_law->ComputeStressAndStressDerivative(C,inv_C,p,T,dTdE,assembleJacobian);
        
        for (unsigned i=0; i<dofs_per_element; i++)
        {
            const unsigned component_i = mFeSystem.system_to_component_index(i).first;
            
            if (assembleJacobian)
            {
                for (unsigned j=0; j<dofs_per_element; j++)
                {
                    const unsigned component_j = mFeSystem.system_to_component_index(j).first;
                    
                    if ((component_i<PRESSURE_COMPONENT_INDEX) &&(component_j<PRESSURE_COMPONENT_INDEX) )
                    {
                        for (unsigned M=0; M<DIM; M++)
                        {
                            for (unsigned N=0; N<DIM; N++)
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
                    else if ((component_i<PRESSURE_COMPONENT_INDEX) &&(component_j==PRESSURE_COMPONENT_INDEX) )
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
                    else if ((component_i==PRESSURE_COMPONENT_INDEX) &&(component_j<PRESSURE_COMPONENT_INDEX) )
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
                    else
                    {
                        // do nothing, ie elementMatrix(i,j)  +=  0 * fe_values.JxW(q_point);;
                    }
                }
            }
            
            if (assembleResidual)
            {
                if (component_i<PRESSURE_COMPONENT_INDEX)
                {
                    elementRhs(i) += - mDensity * mBodyForce(component_i)
                                     * fe_values.shape_value(i,q_point)
                                     * fe_values.JxW(q_point);
                                     
                    for (unsigned M=0; M<DIM; M++)
                    {
                        for (unsigned N=0; N<DIM; N++)
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
    if (assembleResidual)
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
                        
                        if (component_i < PRESSURE_COMPONENT_INDEX)
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


template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::WriteStresses(unsigned counter)
{
    // only write output if the flag mWriteOutput has been set
    if (!mWriteOutput)
    {
        return;
    }

    // create stresses file
    std::stringstream ss;
    ss << this->mOutputDirectoryFullPath << "/finiteelas_solution_" << counter << ".str";
    std::string stress_filename = ss.str();
    std::ofstream stress_output(stress_filename.c_str());
    
    static QGauss<DIM>   quadrature_formula(1);
    const unsigned n_q_points = quadrature_formula.n_quadrature_points;

    FEValues<DIM> fe_values(mFeSystem, quadrature_formula,
                            UpdateFlags(update_values    |
                                        update_gradients |
                                        update_q_points  |     // needed for interpolating u and u' on the quad point
                                        update_JxW_values));
                                        
    //const unsigned dofs_per_element = mFeSystem.dofs_per_cell;

    std::vector< Vector<double> >                  local_solution_values(n_q_points);
    std::vector< std::vector< Tensor<1,DIM> > >    local_solution_gradients(n_q_points);
    
    Tensor<2,DIM> identity;
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
    
    unsigned elem_number = 0;
    typename DoFHandler<DIM>::active_cell_iterator  element_iter = this->mDofHandler.begin_active();
   
    while (element_iter!=this->mDofHandler.end())  
    {
        fe_values.reinit(element_iter); // compute fe values for this element
        fe_values.get_function_values(this->mCurrentSolution, local_solution_values);
        fe_values.get_function_grads(this->mCurrentSolution, local_solution_gradients);

        AbstractIncompressibleMaterialLaw<DIM>* p_material_law;
        if (!mHeterogeneous)
        {
            p_material_law = mMaterialLaws[0];
        }
        else
        {
            unsigned index = GetMaterialLawIndexFromMaterialId(element_iter->material_id());
            p_material_law = mMaterialLaws[index];
        }
                
        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            const std::vector< Tensor<1,DIM> >& grad_u_p = local_solution_gradients[q_point];
            
            Vector<double> x;
            x.reinit(DIM);
            for(unsigned i=0; i<DIM; i++)
            {
                x(i) = local_solution_values[q_point](i);
            }
            
            double p = local_solution_values[q_point](PRESSURE_COMPONENT_INDEX);
            
            static Tensor<2,DIM> F;
            static Tensor<2,DIM> C;
            static Tensor<2,DIM> inv_C;
            static Tensor<2,DIM> inv_F;
            static SymmetricTensor<2,DIM> T;
            
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
            
            p_material_law->ComputeStressAndStressDerivative(C,inv_C,p,T,dTdE,false);
            
            stress_output << elem_number++ << " " << T[0][0] << " " << T[1][0] << " " << T[1][1] << "\n";
        }
        
        element_iter++;
    }
    stress_output.close();
}

template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::DistributeDofs()
{
    this->mDofHandler.distribute_dofs(mFeSystem);
}


template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::ApplyDirichletBoundaryConditions()
{
    // The boundary conditions on the NONLINEAR SYSTEM are x=boundary_values
    // on the boundary nodes. However:
    // The boundary conditions on the LINEAR SYSTEM  Ju=f, where J is the
    // u the negative update vector and f is the residual is
    // u=current_soln-boundary_values on the boundary nodes
    std::map<unsigned,double>  applied_boundary_values;
    
    std::map<unsigned,double>::iterator iter = mBoundaryValues.begin();
    while (iter!=mBoundaryValues.end())
    {
        unsigned dof = iter->first;
        double value = iter->second;
        
        applied_boundary_values[dof] = this->mCurrentSolution(dof)-value;
        iter++;
    }
    
    // don't have access to u (the solution of the linear system) at the moment,
    // so pass in a dummy vector
    Vector<double> dummy(this->mRhsVector.size());
    
    MatrixTools::apply_boundary_values(applied_boundary_values,
                                       this->mSystemMatrix,
                                       dummy,
                                       this->mRhsVector);
}


template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::ComputeNumericalJacobian()
{
    if (mNumericalJacobianMatrix.empty())
    {
        mNumericalJacobianMatrix.reinit(this->mSparsityPattern);
    }
    
    unsigned size = this->mCurrentSolution.size();
    
    // save the current solution
    Vector<double> current_guess = this->mCurrentSolution;
    
    Vector<double> residual(size);
    Vector<double> residual_perturbed(size);
    
    double epsilon= 1e-6;
    
    // save the residual for the current guess
    this->AssembleSystem(true,false);
    residual = this->mRhsVector;
    
    for (unsigned global_column=0; global_column<size; global_column++)
    {
        // reset this->mCurrentSolution...
        this->mCurrentSolution = current_guess;
        //.. and then perturb
        this->mCurrentSolution(global_column) += epsilon;
        
        // compute and store the perturbed residual
        this->AssembleSystem(true,false);
        residual_perturbed = this->mRhsVector;
        
        // compute residual_perturbed - residual
        double one_over_eps=1.0/epsilon;
        for (unsigned i=0; i<size; i++)
        {
            // if value != 0 set in the matrix
            double value = one_over_eps*(residual_perturbed(i) - residual(i));
            if (fabs(value)>1e-12)
            {
                mNumericalJacobianMatrix.set(i,global_column,value);
            }
        }
    }
    
    // reset this->mCurrentSolution to what it was initially
    this->mCurrentSolution = current_guess;
    this->mRhsVector = residual;
    
    
    std::map<unsigned,double>  applied_boundary_values;
    std::map<unsigned,double>::iterator iter = mBoundaryValues.begin();
    while (iter!=mBoundaryValues.end())
    {
        unsigned dof = iter->first;
        double value = iter->second;
        applied_boundary_values[dof] = this->mCurrentSolution(dof)-value;
        iter++;
    }
    
    // don't have access to u (the solution of the linear system) at the moment,
    // so pass in a dummy vector
    Vector<double> dummy(this->mRhsVector.size());
    
    MatrixTools::apply_boundary_values(applied_boundary_values,
                                       mNumericalJacobianMatrix,
                                       dummy,
                                       this->mRhsVector);
}

template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::CompareJacobians()
{
    // compute analytic Jacobian
    this->AssembleSystem(false, true);
    
    // for some reason this has to be computed AFTER the analytic jacobian
    // for them to match (otherwise the boundary condition rows don't match..)
    ComputeNumericalJacobian();
    
    /*
    std::cout << "\nAnalytic Jacobian:\n";
    for(unsigned i=0; i<this->mSystemMatrix.m(); i++)
    {
        for(unsigned j=0; j<this->mSystemMatrix.n(); j++)
        {
            double value = this->mSystemMatrix.el(i,j);
            if(fabs(value)<1e-8)
            {
                value = 0.0;
            }
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "\nNumerical Jacobian:\n";
    for(unsigned i=0; i<this->mSystemMatrix.m(); i++)
    {
        for(unsigned j=0; j<this->mSystemMatrix.n(); j++)
        {
            double value = mNumericalJacobianMatrix.el(i,j);
            if(fabs(value)<1e-8)
            {
                value = 0.0;
            }
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
    */
    
    bool no_difference = true;
    
    std::cout << "\nDifference matrix:\n";
    for (unsigned i=0; i<this->mSystemMatrix.m(); i++)
    {
        for (unsigned j=0; j<this->mSystemMatrix.n(); j++)
        {
            double value = this->mSystemMatrix.el(i,j)-mNumericalJacobianMatrix.el(i,j);
            if (fabs(value)<1e-8)
            {
                value = 0.0;
            }
            else
            {
                #define COVERAGE_IGNORE
                no_difference = false;
                #undef COVERAGE_IGNORE
            }
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
    
    if (!no_difference)
    {
        EXCEPTION("Numerical and analytical Jacobians do not match");
    }
}


template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::OutputResults(unsigned counter)
{
    // only write output if the flag mWriteOutput has been set
    if (!mWriteOutput)
    {
        return;
    }
    
    std::stringstream ss;
    ss << mOutputDirectoryFullPath << "/finiteelas_solution_" << counter << ".gmv";
    std::string filename = ss.str();
    std::ofstream output(filename.c_str());
    
    DataOut<DIM> data_out;
    data_out.attach_dof_handler(this->mDofHandler);
    
    std::vector<std::string> solution_names;
    
    solution_names.push_back("x_displacement");
    if (DIM>1)
    {
        solution_names.push_back("y_displacement");
    }
    if (DIM>2)
    {
        solution_names.push_back("z_displacement");
    }
    
    solution_names.push_back("pressure");
    
    data_out.add_data_vector(this->mCurrentSolution, solution_names);
    data_out.build_patches();
    data_out.write_gmv(output);
}


template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::TakeNewtonStep()
{
    // compute Jacobian
    this->AssembleSystem(false, true);
    
    // solve the linear system
    SolverControl  solver_control(200000, 1e-6, false, true);
    PrimitiveVectorMemory<> vector_memory;
    
    Vector<double> update;
    update.reinit(this->mDofHandler.n_dofs());
    
    SolverGMRES<>::AdditionalData gmres_additional_data(200);
    SolverGMRES<>  gmres(solver_control, vector_memory, gmres_additional_data);
    
    gmres.solve(this->mSystemMatrix, update, this->mRhsVector, PreconditionIdentity());
    
    // deal with hanging nodes - form a continuous solutions
    this->mHangingNodeConstraints.distribute(update);
    
    // save the old current solution
    Vector<double> old_solution = this->mCurrentSolution;
    
    double best_norm_resid = 1e10;
    double best_damping_value = 0.0;
    
    std::vector<double> damping_values;
    damping_values.push_back(0.0);
    damping_values.push_back(0.05);
    for (unsigned i=1; i<=10; i++)
    {
        damping_values.push_back((double)i/10.0);
    }
    
    for (unsigned i=0; i<damping_values.size(); i++)
    {
        this->mCurrentSolution.equ(1.0, old_solution, -damping_values[i], update);
        
        // compute residual
        this->AssembleSystem(true, false);
        double norm_resid = CalculateResidualNorm();
        
        std::cout << "\tTesting s = " << damping_values[i] << ", |f| = " << norm_resid << "\n" << std::flush;
        if (norm_resid < best_norm_resid)
        {
            best_norm_resid = norm_resid;
            best_damping_value = damping_values[i];
        }
    }
    
    
    if (best_damping_value == 0.0)
    {
        std::cout << "\nResidual does not decrease in newton direction, quitting\n" << std::flush;
        assert(0);
    }
    else
    {
        std::cout << "\tBest s = " << best_damping_value << "\n"  << std::flush;
    }
    // implement best update and recalculate residual
    this->mCurrentSolution.equ(1.0, old_solution, -best_damping_value, update);
}


template<unsigned DIM>
void FiniteElasticityAssembler<DIM>::Solve()
{
    if (mMaterialLaws.size()==0)
    {
        EXCEPTION("No material laws have been set");
    }
    
    OutputResults(0);
    
    FormInitialGuess();
    
    // compute residual
    this->AssembleSystem(true, false);
    double norm_resid = CalculateResidualNorm();
    std::cout << "\nNorm of residual is " << norm_resid << "\n";
    
    mNumNewtonIterations = 0;
    unsigned counter = 1;
    
    // use the larger of the tolerances formed from the absolute or
    // relative possibilities
    double tol = NEWTON_ABS_TOL;
    if ( tol < NEWTON_REL_TOL*norm_resid )
    {
        tol = NEWTON_REL_TOL*norm_resid;
    }
    std::cout << "Solving with tolerance " << tol << "\n";
    
    
    while (norm_resid > tol)
    {
        std::cout <<  "\n-------------------\n"
                  <<   "Newton iteration " << counter
                  << ":\n-------------------\n";
        
        TakeNewtonStep();
        
        this->AssembleSystem(true, false);
        norm_resid = CalculateResidualNorm();
        
        std::cout << "Norm of residual is " << norm_resid << "\n";
        
        OutputResults(counter);
        mNumNewtonIterations = counter;
        
        counter++;
        if (counter==20)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("Not converged after 20 newton iterations, quitting");
            #undef COVERAGE_IGNORE
        }
    }
    
    if (norm_resid > tol)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("Failed to converge");
        #undef COVERAGE_IGNORE
    }
    
    // set up mDeformedPosition
    rGetDeformedPosition();
}

template<unsigned DIM>
std::vector<Vector<double> >& FiniteElasticityAssembler<DIM>::rGetDeformedPosition()
{
    for (unsigned i=0; i<DIM; i++)
    {
        mDeformedPosition[i].reinit(this->mpMesh->n_vertices());
    }
    
    DofVertexIterator<DIM> vertex_iter(this->mpMesh, &this->mDofHandler);
    while (!vertex_iter.ReachedEnd())
    {
        unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
        Point<DIM> old_posn = vertex_iter.GetVertex();
        
        for (unsigned i=0; i<DIM; i++)
        {
            mDeformedPosition[i](vertex_index) =   old_posn(i)
                                                 + mCurrentSolution(vertex_iter.GetDof(i));
        }
        
        vertex_iter.Next();
    }
    
    return mDeformedPosition;
}

template<unsigned DIM>
double FiniteElasticityAssembler<DIM>::CalculateResidualNorm()
{
    return this->mRhsVector.norm_sqr()/this->mDofHandler.n_dofs();
}


template<unsigned DIM>
std::vector<Vector<double> >& FiniteElasticityAssembler<DIM>::rGetUndeformedPosition()
{
    // initialise
    mUndeformedPosition.resize(DIM);
    for (unsigned i=0; i<DIM; i++)
    {
        mUndeformedPosition[i].reinit(this->mpMesh->n_vertices());
    }
    
    // populate
    TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
    while (!vertex_iter.ReachedEnd())
    {
        unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
        Point<DIM> old_posn = vertex_iter.GetVertex();
        
        for (unsigned i=0; i<DIM; i++)
        {
            mUndeformedPosition[i](vertex_index) = vertex_iter.GetVertex()(i);
        }
        
        vertex_iter.Next();
    }
    
    return mUndeformedPosition;
}

template<unsigned DIM>
unsigned FiniteElasticityAssembler<DIM>::GetNumNewtonIterations()
{
    return mNumNewtonIterations;
}

#endif // FINITEELASTICITYASSEMBLER_CPP_


