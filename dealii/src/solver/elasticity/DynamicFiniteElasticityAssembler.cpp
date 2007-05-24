#ifndef DYNAMICFINITEELASTICITYASSEMBLER_CPP_
#define DYNAMICFINITEELASTICITYASSEMBLER_CPP_

#include "DynamicFiniteElasticityAssembler.hpp"

#include <dofs/dof_tools.h>

template<unsigned DIM>
DynamicFiniteElasticityAssembler<DIM>::DynamicFiniteElasticityAssembler(Triangulation<DIM>* pMesh,
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
    
    mSolutionAtLastTimestep.reinit( this->mCurrentSolution.size() );
    mSolutionAtLastTimestep = 0;
}

template<unsigned DIM>
DynamicFiniteElasticityAssembler<DIM>::~DynamicFiniteElasticityAssembler()
{}


//////////////////////////////////////////////////////////////////////////////////////////
// AssembleOnElement
//////////////////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>
void DynamicFiniteElasticityAssembler<DIM>::AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
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
    
    static std::vector< Vector<double> >                  local_solution_values(n_q_points);
    static std::vector< std::vector< Tensor<1,DIM> > >    local_solution_gradients(n_q_points);
    
    static std::vector< Vector<double> >                  local_solution_values_last_timestep(n_q_points);
    
    static Tensor<2,DIM> identity;         // how do you do this properly??
    
    static bool first = true;
    
    
    if (first)
    {
        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            local_solution_values[q_point].reinit(DIM+1);
            local_solution_values_last_timestep[q_point].reinit(DIM+1);
            
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
    fe_values.get_function_values(mSolutionAtLastTimestep, local_solution_values_last_timestep);
    fe_values.get_function_grads(this->mCurrentSolution, local_solution_gradients);
    
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law = this->GetMaterialLawForElement(elementIter);

    
    for (unsigned q_point=0; q_point<n_q_points; q_point++)
    {
        const std::vector< Tensor<1,DIM> >& grad_u_p = local_solution_gradients[q_point];
        
        double p = local_solution_values[q_point](this->PRESSURE_COMPONENT_INDEX);
        
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
        
        static SymmetricTensor<2,DIM> T2;
        p_material_law->ComputeStressAndStressDerivative(C,inv_C,p,T,this->dTdE,assembleJacobian);
        
        
        for (unsigned i=0; i<dofs_per_element; i++)
        {
            const unsigned component_i = this->mFeSystem.system_to_component_index(i).first;
            
            if (assembleJacobian)
            {
                for (unsigned j=0; j<dofs_per_element; j++)
                {
                    const unsigned component_j = this->mFeSystem.system_to_component_index(j).first;
                    
                    if ( (component_i<this->PRESSURE_COMPONENT_INDEX) && (component_j<this->PRESSURE_COMPONENT_INDEX) )
                    {
                        // time derivative part
                        elementMatrix(i,j) +=    this->mDensity
                                                 * mDtInverse
                                                 * fe_values.shape_value(i,q_point)
                                                 * fe_values.shape_value(j,q_point)
                                                 * identity[component_i][component_j]
                                                 * fe_values.JxW(q_point);
                                                 
                        // stress part
                        for (unsigned M=0; M<DIM; M++)
                        {
                            for (unsigned N=0; N<DIM; N++)
                            {
                                // dFdU part
                                elementMatrix(i,j) +=   T[M][N]
                                                        * fe_values.shape_grad(j,q_point)[M]
                                                        * fe_values.shape_grad(i,q_point)[N]
                                                        * identity[component_i][component_j]
                                                        * fe_values.JxW(q_point);
                                                        
                                // dTdE part
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
                    else if ((component_i<this->PRESSURE_COMPONENT_INDEX) && (component_j==this->PRESSURE_COMPONENT_INDEX) )
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
                    // time derivative part
                    elementRhs(i) +=   this->mDensity
                                       * (
                                           local_solution_values[q_point](component_i)
                                           - local_solution_values_last_timestep[q_point](component_i)
                                       )
                                       * mDtInverse
                                       * fe_values.shape_value(i,q_point)
                                       * fe_values.JxW(q_point);
                                       
                                       
                                       
                    // body force part
                    elementRhs(i) += - this->mDensity * this->mBodyForce(component_i)
                                     * fe_values.shape_value(i,q_point)
                                     * fe_values.JxW(q_point);
                                     
                    // stress part
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
    
    first = false;
}


template<unsigned DIM>
void DynamicFiniteElasticityAssembler<DIM>::SetTimes(double Tstart, double Tend, double dt)
{
    mTstart = Tstart;
    mTend   = Tend;
    mDt     = dt;
    mDtInverse = 1.0/dt;
    
    if (mTstart >= mTend)
    {
        EXCEPTION("Starting time has to less than ending time");
    }
    if (mDt <= 0)
    {
        EXCEPTION("Time step has to be greater than zero");
    }
    
    assert(mDt <= mTend - mTstart + 1e-10);
    
    mTimesSet = true;
}


template<unsigned DIM>
void DynamicFiniteElasticityAssembler<DIM>::Solve()
{
    if (!mTimesSet)
    {
        EXCEPTION("Start time, end time, dt have not been set. Call SetTimes() before Solve()");
    }
    
    
    double time = mTstart;
    
    this->OutputResultsGMV(0);
    unsigned time_counter=1;
    
    mSolutionAtLastTimestep = this->mCurrentSolution;
    
    
    // compute residual
    this->AssembleSystem(true, false);
    double norm_resid = this->CalculateResidualNorm();
    std::cout << "\nNorm of residual is " << norm_resid << "\n";
    
    // use the larger of the tolerances formed from the absolute or
    // relative possibilities
    double tol = NEWTON_ABS_TOL;
    if ( tol < NEWTON_REL_TOL*norm_resid )
    {
        tol = NEWTON_REL_TOL*norm_resid;
    }
    std::cout << "Solving with tolerance " << tol << "\n";
    
    
    
    while (time < mTend)
    {
        std::cout << "\n===================\n"
                  <<   "Time = " << time
                  << "\n===================\n";
        
        // compute residual
        this->AssembleSystem(true, false);
        double norm_resid = this->CalculateResidualNorm();
        std::cout << "\nNorm of residual is " << norm_resid << "\n";
        
        unsigned newton_counter = 1;
        
        while (norm_resid > tol)
        {
            std::cout <<  "\n-------------------\n"
                      <<   "Newton iteration " << newton_counter
                      << ":\n-------------------\n";
            
            this->TakeNewtonStep();
            
            this->AssembleSystem(true, false);
            norm_resid = this->CalculateResidualNorm();
            
            std::cout << "Norm of residual is " << norm_resid << "\n";
            
            newton_counter++;
            if (newton_counter==20)
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
        
        time += mDt;
        
        mSolutionAtLastTimestep = this->mCurrentSolution;
        
        DofVertexIterator<DIM> vertex_iter(this->mpMesh, &this->mDofHandler);
        while (!vertex_iter.End())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<DIM> old_posn = vertex_iter.GetVertex();
            
            for (unsigned i=0; i<DIM; i++)
            {
                this->mDeformedPosition[i](vertex_index) =   old_posn(i)
                                                           + this->mCurrentSolution(vertex_iter.GetDof(i));
            }
            
            vertex_iter.Next();
        }
        
        this->OutputResultsGMV(time_counter);
        time_counter++;
    }
}




#endif // DYNAMICFINITEELASTICITYASSEMBLER_CPP_


