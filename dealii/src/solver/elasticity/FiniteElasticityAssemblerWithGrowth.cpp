#ifndef FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_
#define FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_

#include "FiniteElasticityAssemblerWithGrowth.hpp"
#include "TriangulationVertexIterator.hpp"

#include <dofs/dof_tools.h>

#define COVERAGE_IGNORE

template<unsigned DIM>
FiniteElasticityAssemblerWithGrowth<DIM>::FiniteElasticityAssemblerWithGrowth(Triangulation<DIM>* pMesh,
                                                                              AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                                                              Vector<double> bodyForce,
                                                                              double density,
                                                                              std::string outputDirectory,
                                                                              AbstractGrowingTumourSourceModel<DIM>* pSourceModel,
                                                                              unsigned degreeOfBasesForPosition,
                                                                              unsigned degreeOfBasesForPressure)  :
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
    mNoRefinement = false;   
    ///////////////////////////////////////////////////////////
    // initialise growth variables
    ///////////////////////////////////////////////////////////
    mGrowthValuesAtVertices.reinit(this->mpMesh->n_vertices());
    mGrowthOdeSystems.resize(this->mpMesh->n_vertices());
    
    for (unsigned i=0; i<this->mpMesh->n_vertices(); i++)
    {
        mGrowthOdeSystems[i] = NULL;
        mGrowthValuesAtVertices(i) = 1.0;
    }
    
    /////////////////////////////////////////////////////////////
    // find growing region and create odes for each node in the
    // region
    /////////////////////////////////////////////////////////////
    bool found_growing_region = false;
    
    typename Triangulation<DIM>::active_cell_iterator element_iter = this->mpMesh->begin_active();
    
    double total_element_size = 0;
    
    // loop over all the elements in the mesh..
    while (element_iter!=this->mpMesh->end())
    {
        total_element_size += element_iter->measure();
        
        unsigned region = element_iter->material_id();
        // check the element is set as GROWING_REGION or NON_GROWING_REGION
        if ( (region!=GROWING_REGION) && (region!=NON_GROWING_REGION))
        {
            std::string err_message("Found element in mesh with is does not have it's material ");
            err_message += "id set to either GROWING_REGION or NON_GROWING_REGION";
            EXCEPTION(err_message);
        }
        
        // if the element is a growing region element (eg tumour)
        if (region == GROWING_REGION)
        {
            found_growing_region = true;
            
            // loop over all vertices..
            for (unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
            {
                unsigned vertex_index = element_iter->vertex_index(i);
                // create a growth ode system for the vertex, assuming one has not
                // been created already, and an evaluation point in the source model
                if (mGrowthOdeSystems[vertex_index]==NULL)
                {
                    mGrowthOdeSystems[vertex_index]
                    = new GrowthByConstantMassOdeSystem<DIM>(this->mDensity,
                                                             vertex_index,
                                                             mpSourceModel);
                                                             
                    Point<DIM> position = element_iter->vertex(i);
                    mpSourceModel->AddEvaluationPoint(vertex_index,
                                                      position);
                }
            }
        }
        
        element_iter++;
    }
    
    mAverageElementVolume = total_element_size/this->mpMesh->n_active_cells();
    //std::cout << "\nAverage element volume = " << mAverageElementVolume << "\n";
    
    // check there was at least one growing element...
    if (!found_growing_region)
    {
        EXCEPTION("No elements in the mesh was labelled as growing");
    }
}

template<unsigned DIM>
FiniteElasticityAssemblerWithGrowth<DIM>::~FiniteElasticityAssemblerWithGrowth()
{
    for (unsigned i=0; i<this->mpMesh->n_vertices(); i++)
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


//////////////////////////////////////////////////////////////////////////////////////////
// AssembleOnElement
//////////////////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
        Vector<double>&       elementRhs,
        FullMatrix<double>&   elementMatrix,
        bool                  assembleResidual,
        bool                  assembleJacobian)
{
    static QGauss<DIM>   quadrature_formula(3);
    static QGauss<DIM-1> face_quadrature_formula(3);
    
    const unsigned n_q_points    = quadrature_formula.n_quadrature_points;
    const unsigned n_face_q_points = face_quadrature_formula.n_quadrature_points;
    

    FE_Q<DIM> linear_fe(1);
    FEValues<DIM> linear_fe_values(linear_fe, quadrature_formula,
                                   UpdateFlags(update_values));
    unsigned linear_dofs_per_element = linear_fe.dofs_per_cell;
    linear_fe_values.reinit(elementIter);

    
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
    
    static std::vector< Vector<double> >               local_solution_values(n_q_points);
    static std::vector< std::vector< Tensor<1,DIM> > > local_solution_gradients(n_q_points);
    
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
    
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law = GetMaterialLawForElement(elementIter); 
    
    double element_volume = 0;

    for (unsigned q_point=0; q_point<n_q_points; q_point++)
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
        
        double growth_term_g = 0;
        for (unsigned i=0; i<linear_dofs_per_element; i++)
        {
            // a little iffy - don't know for certain linear_fe_values.shape_value(i,q)
            // corresponds to the basis function for node i;
            unsigned node_index = elementIter->vertex_index(i);
            growth_term_g += mGrowthValuesAtVertices(node_index)*linear_fe_values.shape_value(i,q_point);
        }

                               
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                full_F[i][j] = identity[i][j] + grad_u_p[i][j];
                F[i][j] = full_F[i][j]/growth_term_g;
            }
        }
        

        element_volume +=  determinant(full_F) * fe_values.JxW(q_point);
        ////////////////////////////////////////////////////////
        // no more changes after this, until userflag setting //
        ////////////////////////////////////////////////////////
        
        
        C = transpose(F) * F;
        inv_C = invert(C);
        inv_F = invert(F);
        
        double detF = determinant(F);        
       
        p_material_law->ComputeStressAndStressDerivative(C,inv_C,p,T,this->dTdE,assembleJacobian);
        
        for (unsigned i=0; i<dofs_per_element; i++)
        {
            const unsigned component_i = this->mFeSystem.system_to_component_index(i).first;
            
            if (assembleJacobian)
            {
                for (unsigned j=0; j<dofs_per_element; j++)
                {
                    const unsigned component_j = this->mFeSystem.system_to_component_index(j).first;
                    
                    if ((component_i<this->PRESSURE_COMPONENT_INDEX) &&(component_j<this->PRESSURE_COMPONENT_INDEX) )
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
                    else if ((component_i<this->PRESSURE_COMPONENT_INDEX) &&(component_j==this->PRESSURE_COMPONENT_INDEX) )
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
                    else if ((component_i==this->PRESSURE_COMPONENT_INDEX) &&(component_j<this->PRESSURE_COMPONENT_INDEX) )
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
                if (component_i<this->PRESSURE_COMPONENT_INDEX)
                {
                    elementRhs(i) += - this->mDensity * this->mBodyForce(component_i)
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
                        const unsigned component_i = this->mFeSystem.system_to_component_index(i).first;
                        
                        if (component_i < this->PRESSURE_COMPONENT_INDEX)
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
    

    if (element_volume > pow(2.0,DIM)*mAverageElementVolume)
    {
        elementIter->set_refine_flag();
    }
    else
    {
        elementIter->clear_refine_flag();
    }

    if (element_volume < (1.0/pow(2.0,DIM))*mAverageElementVolume)
    {
        double mean_g = 0;
        for(unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
        {
            mean_g += mGrowthValuesAtVertices(elementIter->vertex_index(i));
        }
        mean_g /= GeometryInfo<DIM>::vertices_per_cell;

        // don't coarsen if g is greater than 1
        if(mean_g < 1)
        {
            elementIter->set_coarsen_flag();
        }
    }
    else
    {
        elementIter->clear_coarsen_flag();
    }

    first = false;
}





template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::WriteStresses(unsigned counter)
{
    // only write output if the flag mWriteOutput has been set
    if (!this->mWriteOutput)
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

    FE_Q<DIM> linear_fe(1);
    FEValues<DIM> linear_fe_values(linear_fe, quadrature_formula,
                                   UpdateFlags(update_values));
    unsigned linear_dofs_per_element = linear_fe.dofs_per_cell;


    FEValues<DIM> fe_values(this->mFeSystem, quadrature_formula,
                            UpdateFlags(update_values    |
                                        update_gradients |
                                        update_q_points  |     // needed for interpolating u and u' on the quad point
                                        update_JxW_values));
                                        
    //const unsigned dofs_per_element = this->mFeSystem.dofs_per_cell;

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
        linear_fe_values.reinit(element_iter);

        fe_values.reinit(element_iter); // compute fe values for this element
        fe_values.get_function_values(this->mCurrentSolution, local_solution_values);
        fe_values.get_function_grads(this->mCurrentSolution, local_solution_gradients);

        AbstractIncompressibleMaterialLaw<DIM>* p_material_law;
        if (!this->mHeterogeneous)
        {
            p_material_law = this->mMaterialLaws[0];
        }
        else
        {
            unsigned index = this->GetMaterialLawIndexFromMaterialId(element_iter->material_id());
            p_material_law = this->mMaterialLaws[index];
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
            
            double p = local_solution_values[q_point](this->PRESSURE_COMPONENT_INDEX);
            
            
            static Tensor<2,DIM> F;
            static Tensor<2,DIM> C;
            static Tensor<2,DIM> inv_C;
            static Tensor<2,DIM> inv_F;
            static SymmetricTensor<2,DIM> T;
            
            double growth_term_g = 0;
            for (unsigned i=0; i<linear_dofs_per_element; i++)
            {
                // a little iffy - don't know for certain linear_fe_values.shape_value(i,q)
                // corresponds to the basis function for node i;
                unsigned node_index = element_iter->vertex_index(i);
                growth_term_g += mGrowthValuesAtVertices(node_index)*linear_fe_values.shape_value(i,q_point);
            }
            
            for (unsigned i=0; i<DIM; i++)
            {
                for (unsigned j=0; j<DIM; j++)
                {
                    //full_F[i][j] = identity[i][j] + grad_u_p[i][j];
                    F[i][j] = (identity[i][j] + grad_u_p[i][j])/growth_term_g;
                }
            }

            
            C = transpose(F) * F;
            inv_C = invert(C);
            inv_F = invert(F);
            
            p_material_law->ComputeStressAndStressDerivative(C,inv_C,p,T,this->dTdE,false);
            
            stress_output << elem_number++ << " " << T[0][0] << " " << T[1][0] << " " << T[1][1] << "\n";
        }
        
        element_iter++;
    }
    stress_output.close();
}







template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::SetTimes(double Tstart, double Tend, double odeDt)
{
    mTstart = Tstart;
    mTend   = Tend;
    mOdeDt  = odeDt;
    
    if (mTstart >= mTend)
    {
        EXCEPTION("Starting time has to less than ending time");
    }
    if (mOdeDt <= 0)
    {
        EXCEPTION("Time step has to be greater than zero");
    }
    
    assert(mOdeDt <= mTend - mTstart + 1e-10);
    
    mTimesSet = true;
}



template<unsigned DIM>
bool FiniteElasticityAssemblerWithGrowth<DIM>::RefineOvergrownElements()
{
    if(mNoRefinement)
    {
        return false;
    }
    typename Triangulation<DIM>::active_cell_iterator element_iter = this->mpMesh->begin_active();
    
    // determine if there any elements to refine
    bool elements_to_refine = false;
    unsigned ref_counter=0;
    unsigned coa_counter=0;
    unsigned total = 0;

    std::cout << "RC:\n\n";
    while (element_iter!=this->mpMesh->end())
    {
        if (element_iter->refine_flag_set() || element_iter->coarsen_flag_set())
        {
            elements_to_refine = true;
        }

        if (element_iter->refine_flag_set())
        {
            ref_counter++;
//            Point<2> posn = element_iter->center();
//            std::cout << posn[0] << " " << posn[1] << " 1\n";
        }

        if (element_iter->coarsen_flag_set())
        {
            coa_counter++;
//            Point<2> posn = element_iter->center();
//            std::cout << posn[0] << " " << posn[1] << " 2\n";
        }
        
        total++;

        element_iter++;
    } 
    
    std::cout << "r,c,t = " << ref_counter << ", " << coa_counter << " " << total << "\n";
    
    static int counter_again = 0;
       
    if (elements_to_refine)
    {
        std::cout << "\n\n************\nRefining\n************\n\n" << std::flush;
        
        WriteGrowthValuesAtVertices(counter_again+200);
        
        this->AddVectorForInterpolation(&mGrowthValuesAtVertices);
        this->RefineCoarsen();        

        WriteGrowthValuesAtVertices(counter_again+300);

        
        ///////////////////////////////////////////////////////////////
        // recalculate the boundary values
        ///////////////////////////////////////////////////////////////
        this->mBoundaryValues.clear();
        std::vector<bool> component_mask(DIM+1);
        for (unsigned i=0; i<DIM; i++)
        {
            component_mask[i] = true;
        }
        component_mask[DIM] = false;
        VectorTools::interpolate_boundary_values(this->mDofHandler,
                                                 FIXED_BOUNDARY,
                                                 ZeroFunction<DIM>(DIM+1),  // note the "+1" here! - number of components
                                                 this->mBoundaryValues,
                                                 component_mask);
                                                 

        ////////////////////////////////////////////////////////////////
        // a check that all elements still have a correct material id
        ////////////////////////////////////////////////////////////////
        element_iter = this->mpMesh->begin_active();
        while (element_iter!=this->mpMesh->end())
        {
            unsigned region = element_iter->material_id();
            if ( (region!=GROWING_REGION) && (region!=NON_GROWING_REGION))
            {
                element_iter->set_material_id(GROWING_REGION);
                assert(DIM==2);// this is just to point out that the above line might not be correct
            }
            
            element_iter++;
        }
        
        return true;
    }
    else
    {
        return false;
    }
}



template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::Run()
{
    if (!mTimesSet)
    {
        EXCEPTION("Start time, end time, dt have not been set. Call SetTimes() before Solve()");
    }
    
    this->WriteOutput(0);
    this->WriteStresses(0);

    unsigned counter=1;
    
    double time = mTstart;
    while (time < mTend)
    {
        this->mWriteOutput = true;
        this->WriteOutput(counter,false);
        
        // check everything is still fine
        assert(this->mpMesh->n_vertices()==mGrowthOdeSystems.size());
        assert(mGrowthValuesAtVertices.size()==mGrowthOdeSystems.size());
        
        
        std::cout << "=======================\n";
        std::cout << "  Time = " << time << "\n";
        std::cout << "=======================\n";
        
        //////////////////////////////////////////////////////
        // Run the source model up to the next timestep
        //////////////////////////////////////////////////////
        mpSourceModel->Run(time, time+mOdeDt, this);
        
        //////////////////////////////////////////////////////
        // integrate the odes
        //////////////////////////////////////////////////////
        for (unsigned i=0; i<mGrowthOdeSystems.size(); i++)
        {
            if (mGrowthOdeSystems[i]!=NULL)
            {
                mOdeSolver.SolveAndUpdateStateVariable(mGrowthOdeSystems[i],
                                                       time,
                                                       time+mOdeDt,
                                                       mOdeDt/10);
            }
        }
   
        
        //////////////////////////////////////////////////////
        // update the growth values
        //////////////////////////////////////////////////////
        TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();

//            Point<2> centre;
//            Point<2>& position = vertex_iter.GetVertex();
//            Point<2> diff = position - centre;
    
//            double distance_to_centre = std::sqrt(diff.square());
//            double source_value = 5*(distance_to_centre - 0.7);
            
//            double source_value = 0;
//            if((position[0]>20) && (position[0]<30))
//            {
//                source_value = 3;
//            }
 
//            double source_value = 1; //10*exp(-0.5*(position[0]-25)*(position[0]-25));
//            mGrowthValuesAtVertices(vertex_index) += this->mOdeDt*(1.0/2.0)*source_value*mGrowthValuesAtVertices(vertex_index);


            if (mGrowthOdeSystems[vertex_index]!=NULL)
            {
                mGrowthValuesAtVertices(vertex_index) = mGrowthOdeSystems[vertex_index]->rGetStateVariables()[0];
            }
            else
            {
                assert(fabs(mGrowthValuesAtVertices(vertex_index)-1)<1e-6);
            }
            
            
            vertex_iter.Next();
        }

        unsigned num_vertices_before = this->mpMesh->n_vertices();
        
// temporary - just to compute volumes.. - make volume calc/flag setting safe..
        this->AssembleSystem(true,false);

        bool refined = RefineOvergrownElements();
        
        if (refined)
        {
            unsigned num_vertices_after = this->mpMesh->n_vertices();
            
            mGrowthOdeSystems.resize(num_vertices_after);
            for (unsigned i=num_vertices_before; i<num_vertices_after; i++)
            {
                mGrowthOdeSystems[i] = NULL;
            }
            
            typename Triangulation<DIM>::active_cell_iterator element_iter 
               = this->mpMesh->begin_active();
                        
            // loop over all the elements in the mesh..
            while (element_iter!=this->mpMesh->end())
            {
                unsigned region = element_iter->material_id();
                if (region == GROWING_REGION)
                {
                    // loop over all vertices..
                    for (unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
                    {
                        unsigned vertex_index = element_iter->vertex_index(i);
                        // create a growth ode system for the vertex, assuming one has not
                        // been created already, and an evaluation point in the source model
                        if (mGrowthOdeSystems[vertex_index]==NULL)
                        {
                            mGrowthOdeSystems[vertex_index]
                                = new GrowthByConstantMassOdeSystem<DIM>(this->mDensity,
                                                                         vertex_index,
                                                                         mpSourceModel);
                            
                            mGrowthOdeSystems[vertex_index]->rGetStateVariables()[0] 
                                = mGrowthValuesAtVertices(vertex_index);

                            Point<DIM> position = element_iter->vertex(i);
                            mpSourceModel->AddEvaluationPoint(vertex_index,
                                                              position);
                        }
                    }
                }
                
                element_iter++;
            }           
            
            // temporary, for debugging
            this->mWriteOutput = true;            
            this->WriteOutput(counter+100);
            this->WriteStresses(counter+100);
            WriteGrowthValuesAtVertices(counter+100); 
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
        mpSourceModel->UpdateEvaluationPointsNewPosition(this->mDeformedPosition);
        
        ////////////////////////////////////////////////////////
        // output results
        ////////////////////////////////////////////////////////
        this->mWriteOutput = true;
        this->WriteOutput(counter);
        this->WriteStresses(counter);
        WriteGrowthValuesAtVertices(counter);
        
        counter++;
        time += mOdeDt;
    }
}


template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::WriteGrowthValuesAtVertices(unsigned counter)
{
    // only write output if the flag mWriteOutput has been set
    if (!this->mWriteOutput)
    {
        return;
    }
    
    std::stringstream ss;
    ss << this->mOutputDirectoryFullPath << "/finiteelas_solution_" << counter << ".growth";
    std::string growth_vals_filename = ss.str();
    std::ofstream growth_vals_output(growth_vals_filename.c_str());

    TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
    while (!vertex_iter.ReachedEnd())
    {
        Point<DIM> posn = vertex_iter.GetVertex();
        unsigned index = vertex_iter.GetVertexGlobalIndex();
        
        growth_vals_output << index << " "; 
        for (unsigned i=0; i<DIM; i++)
        {
            growth_vals_output << posn[i] << " ";
        }
        growth_vals_output << mGrowthValuesAtVertices(index)<<"\n";
        vertex_iter.Next();
    } 
    
    growth_vals_output.close();
}


template<unsigned DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::SetNoRefinement()
{
    mNoRefinement = true;
}

#undef COVERAGE_IGNORE

#endif // FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_

