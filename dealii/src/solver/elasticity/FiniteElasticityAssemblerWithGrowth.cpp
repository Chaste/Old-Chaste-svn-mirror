#ifndef FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_
#define FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_

#include "FiniteElasticityAssemblerWithGrowth.hpp"
#include "TriangulationVertexIterator.hpp"

#include <dofs/dof_tools.h>


// to be removed later
#include "ConstantTumourSourceModel.hpp"



template<int DIM>
FiniteElasticityAssemblerWithGrowth<DIM>::FiniteElasticityAssemblerWithGrowth(Triangulation<DIM>* pMesh,
                                                                              AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                                                              Vector<double> bodyForce,
                                                                              double density,
                                                                              std::string outputDirectory,
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
    
    mpSourceModel = new ConstantTumourSourceModel<DIM>(0.1);

    /////////////////////////////////////////////////////////////
    // find growing region and create odes for each node in the
    // region 
    /////////////////////////////////////////////////////////////
    bool found_growing_region = false;
    unsigned eval_point_index = 0;
    
    Triangulation<2>::active_cell_iterator element_iter = this->mpMesh->begin_active();

    // loop over all the elements in the mesh..
    while(element_iter!=this->mpMesh->end())
    {
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
            for (unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
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
//            const Point<2> vector_to_centre = (element_iter->vertex(vertex) - centre);
//            const double distance_from_centre = std::sqrt(vector_to_centre.square());
//          
//            if (distance_from_centre < 0.2)
//            {
//                element_iter->set_refine_flag();
//                std::cout << "refining element ";
//                break;
//            }
        element_iter++;    
    }
      //  mesh.execute_coarsening_and_refinement();

    // check there was at least one growing element...
    if(!found_growing_region)
    {
        EXCEPTION("No elements in the mesh was labelled as growing");
    }
}

template<int DIM>
FiniteElasticityAssemblerWithGrowth<DIM>::~FiniteElasticityAssemblerWithGrowth()
{       
    for(unsigned i=0; i<this->mpMesh->n_vertices(); i++)
    {
        delete mGrowthOdeSystems[i];
    }
}

template<int DIM>
bool FiniteElasticityAssemblerWithGrowth<DIM>::IsGrowingNode(unsigned vertexIndex)
{
    assert(vertexIndex < mGrowthOdeSystems.size());
    return (mGrowthOdeSystems[vertexIndex]!=NULL);
}


//////////////////////////////////////////////////////////////////////////////////////////
// AssembleOnElement
//////////////////////////////////////////////////////////////////////////////////////////
template<int DIM>
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
  
        double growth_term_g = 0.25*(   mGrowthValuesAtVertices(n0)
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
        ////////////////////////////////////////////////////////
        //              no more changes after this            //
        ////////////////////////////////////////////////////////
        
        
        C = transpose(F) * F;
        inv_C = invert(C);
        inv_F = invert(F);

        double detF = determinant(F);
        
        static SymmetricTensor<2,DIM> T2;
        this->mpMaterialLaw->ComputeStressAndStressDerivative(C,inv_C,p,T,this->dTdE,assembleJacobian);

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

    first = false;
}



template<int DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::RefineOvergrownElements()
{
}
    


template<int DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::SetTimes(double Tstart, double Tend, double odeDt)
{
    mTstart = Tstart;
    mTend   = Tend;
    mOdeDt     = odeDt;
    
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


template<int DIM>
void FiniteElasticityAssemblerWithGrowth<DIM>::Run()
{
    if(!mTimesSet)
    {
        EXCEPTION("Start time, end time, dt have not been set. Call SetTimes() before Solve()");
    }
    

    this->OutputResults(0);


    unsigned counter=1;
    double time = mTstart;
    while(time < mTend)
    {
        std::cout << "===========================\n";
        std::cout << "Time = " << time << "\n"; 
        std::cout << "===========================\n";
        
        //////////////////////////////////////////////////////
        // Run the source model up to the next timestep
        //////////////////////////////////////////////////////
        mpSourceModel->Run(time, time+mOdeDt);


        //////////////////////////////////////////////////////
        // integrate the odes
        //////////////////////////////////////////////////////
        for(unsigned i=0; i<this->mpMesh->n_vertices(); i++)
        {
            if(mGrowthOdeSystems[i]!=NULL)
            {
                mOdeSolver.SolveAndUpdateStateVariable(mGrowthOdeSystems[i],
                                                       time,
                                                       time+mOdeDt,
                                                       mOdeDt);
            }
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
            
            vertex_iter.Next();
        }


//        for(unsigned i=0; i<this->mpMesh->n_vertices(); i++)
//        {
//            std::cout << mGrowthValuesAtVertices(i) << " ";
//        }
//        std::cout << "\n";
    
    
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
        counter++;
        
        time += mOdeDt;
    }
};
#endif // FINITEELASTICITYASSEMBLERWITHGROWTH_CPP_


