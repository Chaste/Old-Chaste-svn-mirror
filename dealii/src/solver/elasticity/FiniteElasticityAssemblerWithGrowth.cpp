#ifndef FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_
#define FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_

#include "FiniteElasticityAssemblerWithGrowth.hpp"
#include "TriangulationVertexIterator.hpp"

#include <dofs/dof_tools.h>



  
#include <numerics/solution_transfer.h>



template<unsigned DIM>
FiniteElasticityAssemblerWithGrowth<DIM>::FiniteElasticityAssemblerWithGrowth(Triangulation<DIM>* pMesh,
                                                                              AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                                                              Vector<double> bodyForce,
                                                                              double density,
                                                                              std::string outputDirectory,
                                                                              AbstractGrowingTumourSourceModel<DIM>* pSourceModel,  
                                                                              unsigned degreeOfBasesForPosition,
                                                                              unsigned degreeOfBasesForPressure
                                                                              )  :
        FiniteElasticityAssembler<DIM>(pMesh,
                                       pMaterialLaw,
                                       bodyForce,
                                       density, 
                                       outputDirectory,
                                       degreeOfBasesForPosition,
                                       degreeOfBasesForPressure)
{
    assert(pSourceModel!=NULL);
    mpSourceModel = pSourceModel;
    
    mTimesSet = false;
    
    ///////////////////////////////////////////////////////////
    // initialise growth variables
    ///////////////////////////////////////////////////////////
    mGrowthValuesAtVertices.reinit(this->mpMesh->n_vertices());
    mGrowthOdeSystems.resize(this->mpMesh->n_vertices());
    
    for(unsigned i=0; i<this->mpMesh->n_vertices(); i++)
    {
        mGrowthOdeSystems[i] = NULL;
        mGrowthValuesAtVertices(i) = 1.0;
    }
    




    /////////////////////////////////////////////////////////////
    // find growing region and create odes for each node in the
    // region 
    /////////////////////////////////////////////////////////////
    bool found_growing_region = false;
    unsigned eval_point_index = 0;
    
    typename Triangulation<DIM>::active_cell_iterator element_iter = this->mpMesh->begin_active();

    double total_element_size = 0;
    
    // loop over all the elements in the mesh..
    while(element_iter!=this->mpMesh->end())
    {
        total_element_size += element_iter->measure();
        
        unsigned region = element_iter->material_id();
        // check the element is set as GROWING_REGION or NON_GROWING_REGION
        if( (region!=GROWING_REGION) && (region!=NON_GROWING_REGION))
        {
            std::string err_message("Found element in mesh with is does not have it's material ");
            err_message += "id set to either GROWING_REGION or NON_GROWING_REGION";
            EXCEPTION(err_message);
        }
        
        // if the element is a growing region element (eg tumour)
        if(region == GROWING_REGION)
        {
            found_growing_region = true;
            
            // loop over all vertices..
            for(unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
            {
                unsigned vertex_index = element_iter->vertex_index(i);
                // create a growth ode system for the vertex, assuming one has not
                // been created already, and an evaluation point in the source model
                if(mGrowthOdeSystems[vertex_index]==NULL)
                {
                    mGrowthOdeSystems[vertex_index] 
                         = new GrowthByConstantMassOdeSystem<DIM>(this->mDensity,
                                                                  eval_point_index,
                                                                  mpSourceModel);

                    Point<DIM> position = element_iter->vertex(i);
                    mpSourceModel->AddEvaluationPoint(eval_point_index,
                                                      position,
                                                      vertex_index);
                    eval_point_index++;                
                }
            }
        }        

        element_iter++;    
    }

    mAverageElementVolume = total_element_size/this->mpMesh->n_active_cells();
    std::cout << "\nAverage element volume = " << mAverageElementVolume << "\n";

    // check there was at least one growing element...
    if(!found_growing_region)
    {
        EXCEPTION("No elements in the mesh was labelled as growing");
    }
    
    mUseRefinement = true;
}

template<unsigned DIM>
FiniteElasticityAssemblerWithGrowth<DIM>::~FiniteElasticityAssemblerWithGrowth()
{       
    for(unsigned i=0; i<this->mpMesh->n_vertices(); i++)
    {
        delete mGrowthOdeSystems[i];
    }
}

template<unsigned DIM>
bool FiniteElasticityAssemblerWithGrowth<DIM>::IsGrowingNode(unsigned vertexIndex)
{
    assert(vertexIndex < mGrowthOdeSystems.size());
    return (mGrowthOdeSystems[vertexIndex]!=NULL);
}

template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::DoNotUseRefinement()
{
    mUseRefinement = false;
}

template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::UseRefinement()
{
    mUseRefinement = true;
}


//////////////////////////////////////////////////////////////////////////////////////////
// AssembleOnElement
//////////////////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter, 
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
    FEValues<DIM> fe_values(this->mFeSystem, quadrature_formula, 
                            UpdateFlags(update_values    |
                                        update_gradients |
                                        update_q_points  |     // needed for interpolating u and u' on the quad point
                                        update_JxW_values));

    // would want this to be static too (slight speed up), but causes errors
    // in debug mode (upon destruction of the class, in 2d, or something)
    FEFaceValues<DIM> fe_face_values(this->mFeSystem, face_quadrature_formula, 
                                     UpdateFlags(update_values         |
                                                 update_q_points       |
                                                 update_normal_vectors |
                                                 update_JxW_values));

    
    const unsigned dofs_per_element = this->mFeSystem.dofs_per_cell;


    static std::vector< unsigned >                     local_dof_indices(dofs_per_element);
    static std::vector< Vector<double> >               local_solution_values(n_q_points);
    static std::vector< std::vector< Tensor<1,DIM> > > local_solution_gradients(n_q_points);

    static Tensor<2,DIM> identity;         

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
    fe_values.get_function_values(this->mCurrentSolution, local_solution_values);
    fe_values.get_function_grads(this->mCurrentSolution, local_solution_gradients);        

    AbstractIncompressibleMaterialLaw<DIM>* p_material_law;        
    if(!this->mHeterogeneous)
    {
        p_material_law = this->mMaterialLaws[0];
    }
    else
    {
        unsigned index = this->GetMaterialLawIndexFromMaterialId(elementIter->material_id());
        p_material_law = this->mMaterialLaws[index];
    }
    
    double element_volume = 0;

    for(unsigned q_point=0; q_point<n_q_points; q_point++)
    {       
        const std::vector< Tensor<1,DIM> >& grad_u_p = local_solution_gradients[q_point];
        
        double p = local_solution_values[q_point](this->PRESSURE_COMPONENT_INDEX);

        ////////////////////////////////////////////////////////
        //      this is the part that changes for growth      //
        ////////////////////////////////////////////////////////
        static Tensor<2,DIM> full_F;       // F_e*F_g (F for both elastic and growth parts
        static Tensor<2,DIM> F;            // just F_e, the elastic part
        static Tensor<2,DIM> C;            // everything else defined from F_e
        static Tensor<2,DIM> inv_C;
        static Tensor<2,DIM> inv_F;
        static SymmetricTensor<2,DIM> T;
        
        unsigned n0 = elementIter->vertex_index(0);
        unsigned n1 = elementIter->vertex_index(1);
        unsigned n2 = elementIter->vertex_index(2);
        unsigned n3 = elementIter->vertex_index(3);
        
        // std::cout << n0 << " " << n1 << " " << n2 << " "<< n3 << "\n";

/// TODO: proper interpolation:
  
        double growth_term_g = 0.25*(
                                        mGrowthValuesAtVertices(n0)
                                      + mGrowthValuesAtVertices(n1)
                                      + mGrowthValuesAtVertices(n2)
                                      + mGrowthValuesAtVertices(n3) 
                                    );


        for(unsigned i=0; i<DIM; i++)
        {
            for(unsigned j=0; j<DIM; j++)
            {
                full_F[i][j] = identity[i][j] + grad_u_p[i][j];
                F[i][j] = full_F[i][j]/growth_term_g;
            }
        }

        assert(DIM==2);
        element_volume +=  determinant(full_F) * fe_values.JxW(q_point);
        ////////////////////////////////////////////////////////
        // no more changes after this, until userflag setting //
        ////////////////////////////////////////////////////////
        
        
        C = transpose(F) * F;
        inv_C = invert(C);
        inv_F = invert(F);

        double detF = determinant(F);

        
        static SymmetricTensor<2,DIM> T2;
        p_material_law->ComputeStressAndStressDerivative(C,inv_C,p,T,this->dTdE,assembleJacobian);

        for(unsigned i=0; i<dofs_per_element; i++)
        {
            const unsigned component_i = this->mFeSystem.system_to_component_index(i).first;
            
            if(assembleJacobian)
            {
                for(unsigned j=0; j<dofs_per_element; j++)
                {
                    const unsigned component_j = this->mFeSystem.system_to_component_index(j).first;
        
                    if((component_i<this->PRESSURE_COMPONENT_INDEX) &&(component_j<this->PRESSURE_COMPONENT_INDEX) )
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
                                                              * this->dTdE[M][N][P][Q]
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
                    else if((component_i<this->PRESSURE_COMPONENT_INDEX) &&(component_j==this->PRESSURE_COMPONENT_INDEX) )
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
                    else if((component_i==this->PRESSURE_COMPONENT_INDEX) &&(component_j<this->PRESSURE_COMPONENT_INDEX) )
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
                if(component_i<this->PRESSURE_COMPONENT_INDEX)
                {
                    elementRhs(i) += - this->mDensity * this->mBodyForce(component_i)
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
        for(unsigned face_index=0; face_index<GeometryInfo<DIM>::faces_per_cell; face_index++)
        {
            if(elementIter->face(face_index)->boundary_indicator() == NEUMANN_BOUNDARY)
            {
                fe_face_values.reinit(elementIter, face_index);
  
                for(unsigned q_point=0; q_point<n_face_q_points; q_point++)
                {
                    Vector<double> neumann_traction(DIM); // zeros
    
                    neumann_traction(1)=0.0;
    
                    for(unsigned i=0; i<dofs_per_element; i++)
                    {
                        const unsigned component_i = this->mFeSystem.system_to_component_index(i).first;

                        if(component_i < this->PRESSURE_COMPONENT_INDEX)
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



    if(element_volume > pow(2,DIM)*mAverageElementVolume)
    {
        elementIter->set_user_flag();
    }
    else
    {
        elementIter->clear_user_flag();
    }
 
    first = false;
}



template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::WriteBasicOutput(unsigned counter)
{
    // only write output if the flag mWriteOutput has been set
    if(!this->mWriteOutput)
    {
        return;
    }
 

    ///////////////////////////////////////////////////////////////////// 
    // create an node file, by looping over vertices and writing
    //   vertex_index x y [z]
    ///////////////////////////////////////////////////////////////////// 
    std::stringstream ss;
    ss << this->mOutputDirectoryFullPath << "/finiteelas_solution_" << counter << ".nodes";
    std::string nodes_filename = ss.str();
    std::ofstream nodes_output(nodes_filename.c_str());
    
    DofVertexIterator<DIM> vertex_iter(this->mpMesh, &(this->mDofHandler));
    while(!vertex_iter.ReachedEnd())
    {
        unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
        nodes_output << vertex_index << " ";

        Point<DIM> old_posn = vertex_iter.GetVertex();
        
        Point<DIM> new_posn;
        for(unsigned i=0; i<DIM; i++)
        {
            new_posn(i) = old_posn(i) + this->mCurrentSolution(vertex_iter.GetDof(i));
            nodes_output << new_posn(i) << " "; 
        }
        nodes_output << "\n";
        vertex_iter.Next();
    }
    
    nodes_output.close();


    ///////////////////////////////////////////////////////////////////// 
    // create an element file, by looping over elements and writing
    //   node1 node2  .... nodeN tumour
    // where node_i is the vertex index and tumour = 0 or 1 indicating
    // whether in grwoing region or not
    ///////////////////////////////////////////////////////////////////// 
    std::stringstream ss2;
    ss2 << this->mOutputDirectoryFullPath << "/finiteelas_solution_" << counter << ".elem";
    std::string elem_filename = ss2.str();
    std::ofstream elem_output(elem_filename.c_str());
    
    typename Triangulation<DIM>::active_cell_iterator element_iter = this->mpMesh->begin_active();
    while(element_iter!=this->mpMesh->end())
    {
        // loop over all vertices..
        for(unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
        {
            elem_output << element_iter->vertex_index(i) << " ";
        }

        unsigned region = element_iter->material_id();
        if(region == GROWING_REGION)
        {
            elem_output << 1 << " ";
        }
        else
        {
            elem_output << 0 << " ";
        }
        elem_output << "\n";
        element_iter++;
    }
    elem_output.close();

}



template<unsigned DIM>
bool FiniteElasticityAssemblerWithGrowth<DIM>::RefineOvergrownElements(unsigned i)
{
    typename Triangulation<DIM>::active_cell_iterator element_iter = this->mpMesh->begin_active();

    // determine if there any elements to refine
    bool elements_to_refine = false;
    while(element_iter!=this->mpMesh->end())
    {
        if(element_iter->user_flag_set())
        {
            element_iter->set_refine_flag();
            elements_to_refine = true;
        }
        element_iter++;    
    }


 
    if(elements_to_refine)
    {
        this->mpMesh->prepare_coarsening_and_refinement();
  
        SolutionTransfer<DIM,double> solution_transfer(this->mDofHandler);
        solution_transfer.prepare_for_coarsening_and_refinement(this->mCurrentSolution);
    
        SolutionTransfer<DIM,double> solution_transfer2(this->mDofHandler);
        solution_transfer2.prepare_for_coarsening_and_refinement(this->mRhsVector);
    
        SolutionTransfer<DIM,double> solution_transfer3(this->mDofHandler);
        solution_transfer3.prepare_for_coarsening_and_refinement(mGrowthValuesAtVertices);
    
    
      
        std::cout << "\n\n************\nRefining\n*************\n\n" << std::flush;
        this->mpMesh->execute_coarsening_and_refinement();
    
    
    
        this->mDofHandler.distribute_dofs(this->mFeSystem);
      
        Vector<double> temp(this->mDofHandler.n_dofs());
        solution_transfer.interpolate(this->mCurrentSolution, temp);
        this->mCurrentSolution = temp;
        
        temp = 0;
        solution_transfer2.interpolate(this->mRhsVector, temp);
        this->mRhsVector = temp;
    
    
      
      
    Vector<double> temp2(this->mpMesh->n_vertices());
    assert(temp2.size() > mGrowthValuesAtVertices.size());
    for(unsigned i=0; i<mGrowthValuesAtVertices.size(); i++)
    {
        temp2(i) = mGrowthValuesAtVertices(i);
    }
    
    // this should be interpolated!!
    for(unsigned i=mGrowthValuesAtVertices.size(); i<temp2.size(); i++)
    {
        temp2(i) = 1.0;
    }
    
    mGrowthValuesAtVertices = temp2;
    
    
    
    
        this->mHangingNodeConstraints.clear();
        DoFTools::make_hanging_node_constraints(this->mDofHandler,
                                                this->mHangingNodeConstraints);
        this->mHangingNodeConstraints.close();  
        this->mHangingNodeConstraints.distribute(this->mCurrentSolution);
      
    
        // form sparsity pattern
        this->mSparsityPattern.reinit(this->mDofHandler.n_dofs(), 
                                      this->mDofHandler.n_dofs(),
                                      this->mDofHandler.max_couplings_between_dofs());
    
        DoFTools::make_sparsity_pattern(this->mDofHandler, this->mSparsityPattern);
    
        // see dealii tutorial 2
        this->mHangingNodeConstraints.condense(this->mSparsityPattern);
    
    
        this->mSparsityPattern.compress();
        
        // initialise vectors and matrices
        this->mSystemMatrix.reinit(this->mSparsityPattern);
        this->mRhsVector.reinit(this->mDofHandler.n_dofs());
    

    
    
        ///////////////////////////////////////////////////////////////
        // recalculate the boundary values
        ///////////////////////////////////////////////////////////////
        this->mBoundaryValues.clear();
        std::vector<bool> component_mask(DIM+1);
        for(unsigned i=0; i<DIM; i++)
        {
            component_mask[i] = true;
        }
        component_mask[DIM] = false;
        VectorTools::interpolate_boundary_values(this->mDofHandler,
                                                 FIXED_BOUNDARY,
                                                 ZeroFunction<DIM>(DIM+1),  // note the "+1" here! - number of components
                                                 this->mBoundaryValues,
                                                 component_mask);
    
        this->mNumericalJacobianMatrix.reinit(this->mSparsityPattern);    
        
        
        ////////////////////////////////////////////////////////////////
        // a check that all elements still have a correct material id
        ////////////////////////////////////////////////////////////////
        element_iter = this->mpMesh->begin_active();
        while(element_iter!=this->mpMesh->end())
        {
            unsigned region = element_iter->material_id();
            if( (region!=GROWING_REGION) && (region!=NON_GROWING_REGION))
            {
                assert(0);
            }     
            
            element_iter++;   
        }        
        
        return true;
    }
    else
    {
        std::cout << "\n\n***Not Refining****\n\n" << std::flush;
        return false;
    }
}
    


template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::SetTimes(double Tstart, double Tend, double odeDt)
{
    mTstart = Tstart;
    mTend   = Tend;
    mOdeDt  = odeDt;
    
    if(mTstart >= mTend)
    {
        EXCEPTION("Starting time has to less than ending time");
    }
    if(mOdeDt <= 0)
    {
        EXCEPTION("Time step has to be greater than zero");
    }
    
    assert(mOdeDt <= mTend - mTstart + 1e-10);
    
    mTimesSet = true;
}


template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::Run()
{
    if(!mTimesSet)
    {
        EXCEPTION("Start time, end time, dt have not been set. Call SetTimes() before Solve()");
    }
    

    this->OutputResults(0);
    WriteBasicOutput(0);
    
    unsigned counter=1;
    double time = mTstart;
    while(time < mTend)
    {
        // check everything is still fine
        assert(this->mpMesh->n_vertices()==mGrowthOdeSystems.size());
        assert(mGrowthValuesAtVertices.size()==mGrowthOdeSystems.size());

        
        std::cout << "===========================\n";
        std::cout << "Time = " << time << "\n"; 
        std::cout << "===========================\n";
        
        //////////////////////////////////////////////////////
        // Run the source model up to the next timestep
        //////////////////////////////////////////////////////
        mpSourceModel->Run(time, time+mOdeDt, this);

        //////////////////////////////////////////////////////
        // integrate the odes
        //////////////////////////////////////////////////////
        for(unsigned i=0; i<mGrowthOdeSystems.size(); i++)
        {
            if(mGrowthOdeSystems[i]!=NULL)
            {
                mOdeSolver.SolveAndUpdateStateVariable(mGrowthOdeSystems[i],
                                                       time,
                                                       time+mOdeDt,
                                                       mOdeDt);
            }
        }
        
        

        unsigned num_vertices_before = this->mpMesh->n_vertices();

// temporary - just to compute volumes.. - make volume calc/flag setting safe..
        this->AssembleSystem(true,false);
        bool refined = RefineOvergrownElements(counter);

        if(refined)
        {
            unsigned num_vertices_after = this->mpMesh->n_vertices();
            
            mGrowthOdeSystems.resize(num_vertices_after);
            for(unsigned i=num_vertices_before; i<num_vertices_after; i++)
            {
                mGrowthOdeSystems[i] = NULL;
            }

            
            typename Triangulation<DIM>::active_cell_iterator element_iter = this->mpMesh->begin_active();

            unsigned eval_point_index = mpSourceModel->GetNumEvaluationPoints();

            // loop over all the elements in the mesh..
            while(element_iter!=this->mpMesh->end())
            {
                unsigned region = element_iter->material_id();
                if(region == GROWING_REGION)
                {
                    // loop over all vertices..
                    for(unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
                    {
                        unsigned vertex_index = element_iter->vertex_index(i);
                        // create a growth ode system for the vertex, assuming one has not
                        // been created already, and an evaluation point in the source model
                        if(mGrowthOdeSystems[vertex_index]==NULL)
                        {
                            mGrowthOdeSystems[vertex_index] 
                                 = new GrowthByConstantMassOdeSystem<DIM>(this->mDensity,
                                                                          eval_point_index,
                                                                          mpSourceModel);

//todo: set state var to mGrowthValuesAtVertices(i), once this has been interped correctly

        
                            Point<DIM> position = element_iter->vertex(i);
                            mpSourceModel->AddEvaluationPoint(eval_point_index,
                                                              position,
                                                              vertex_index);
                            eval_point_index++;                
                        }
                    }
                }        
        
                element_iter++;    
            }

        this->mWriteOutput = true;
        this->OutputResults(counter);
        WriteBasicOutput(counter);
        counter++;

        }
            
              
        //////////////////////////////////////////////////////
        // update the growth values
        //////////////////////////////////////////////////////
        TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            if(mGrowthOdeSystems[vertex_index]!=NULL)
            {
                mGrowthValuesAtVertices(vertex_index) = mGrowthOdeSystems[vertex_index]->rGetStateVariables()[0];
            }
            else
            {
                assert(fabs(mGrowthValuesAtVertices(vertex_index)-1)<1e-6);
            }
            
            vertex_iter.Next();
        }
        
    
        ////////////////////////////////////////////////////////
        // solve the (quasi-static) finite elasticity problem
        ////////////////////////////////////////////////////////
        this->mWriteOutput = false;
        this->Solve();
        
        
        ////////////////////////////////////////////////////////
        // update the new position at the evaluation points 
        // in the source model
        ////////////////////////////////////////////////////////
        mpSourceModel->UpdateEvaluationPointsNewPosition(this);

        ////////////////////////////////////////////////////////
        // output results
        ////////////////////////////////////////////////////////
        this->mWriteOutput = true;
        this->OutputResults(counter);
        WriteBasicOutput(counter);
        counter++;
        time += mOdeDt;
    }
};
#endif // FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_


