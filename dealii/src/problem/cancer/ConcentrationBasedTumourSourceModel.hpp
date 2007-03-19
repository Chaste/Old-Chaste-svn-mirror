#ifndef CONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_
#define CONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_

#include "AbstractGrowingTumourSourceModel.hpp"
#include "TriangulationVertexIterator.hpp"

// todo: lots. move laplaces eqn class to somewhere sensible (refactor assemblers 
// this, test update mesh, test copy mesh, test source terms.


template<unsigned DIM>
class LaplacesEquation
{
protected:
    Triangulation<DIM>*   mpMesh;
    FE_Q<DIM>             mFe;            
    DoFHandler<DIM>       mDofHandler;

    ConstraintMatrix      mHangingNodeConstraints;

    SparsityPattern       mSparsityPattern;

    SparseMatrix<double>  mSystemMatrix;
    Vector<double>        mRhsVector;
    Vector<double>        mSolution;
    
    Vector<double>        mSolutionAtVertices;

    virtual void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter, 
                                   Vector<double>&                                 elementRhs,
                                   FullMatrix<double>&                             elementMatrix)
    {
        static QGauss<DIM>   quadrature_formula(3);
        const unsigned n_q_points    = quadrature_formula.n_quadrature_points;
        
        
        // would want this to be static too (slight speed up), but causes errors
        // in debug mode (upon destruction of the class, in 2d, or something)
        FEValues<DIM> fe_values(mFe, quadrature_formula, 
                                UpdateFlags(update_values    |
                                            update_gradients |
                                            update_JxW_values));

        
        const unsigned dofs_per_element = mFe.dofs_per_cell;
    
        std::vector<unsigned> local_dof_indices(dofs_per_element);

        elementMatrix = 0;
        elementRhs = 0;
          
        elementIter->get_dof_indices(local_dof_indices);
    
        fe_values.reinit(elementIter); // compute fe values for this element
    
        double source_term = 1;

// !!!!!!!!!!!!!!!! GROWING_REGION == 99
        if(elementIter->material_id()==99)
        {
            source_term = 1;
        }
    
        for(unsigned q_point=0; q_point<n_q_points; q_point++)
        {       
            for(unsigned i=0; i<dofs_per_element; i++)
            {
                for(unsigned j=0; j<dofs_per_element; j++)
                {
                    elementMatrix(i,j) +=    fe_values.shape_grad(i,q_point) 
                                           * fe_values.shape_grad(j,q_point) 
                                           * fe_values.JxW(q_point);
                }
                
                
                elementRhs(i) +=    fe_values.shape_value(i,q_point) 
                                  * source_term
                                  * fe_values.JxW (q_point);
            }
        }
    }


    void AssembleSystem()
    {
        const unsigned       dofs_per_element = mFe.dofs_per_cell;
      
        FullMatrix<double>   element_matrix(dofs_per_element, dofs_per_element);
        Vector<double>       element_rhs(dofs_per_element);
    
        // the dofs associated with the nodes of an element
        std::vector<unsigned> local_dof_indices(dofs_per_element);
    
        typename DoFHandler<DIM>::active_cell_iterator  element_iter = mDofHandler.begin_active();
        
        mRhsVector = 0;
        mSystemMatrix = 0;
 
        while(element_iter!=mDofHandler.end())  
        {
            // zero the small matrix and vector
            element_matrix = 0;
            element_rhs = 0;
          
            element_iter->get_dof_indices(local_dof_indices);
    
            AssembleOnElement(element_iter,
                              element_rhs, 
                              element_matrix);                    
    
            for(unsigned i=0; i<dofs_per_element; i++)
            {
                for(unsigned j=0; j<dofs_per_element; j++)
                {
                    mSystemMatrix.add(local_dof_indices[i],
                                      local_dof_indices[j], 
                                      element_matrix(i,j));
                }

                mRhsVector(local_dof_indices[i]) += element_rhs(i);
            }
    
            element_iter++;
        }

        mHangingNodeConstraints.condense(mSystemMatrix);
        mHangingNodeConstraints.condense(mRhsVector);
        
        ApplyDirichletBoundaryConditions();
    }
    

    void ApplyDirichletBoundaryConditions()
    {
        std::map<unsigned int,double> boundary_values;
        VectorTools::interpolate_boundary_values(mDofHandler,
                                                 0,
                                                 ZeroFunction<2>(),
                                                 boundary_values);
    
        MatrixTools::apply_boundary_values(boundary_values,
                                           mSystemMatrix,
                                           mSolution,
                                           mRhsVector);
    }

public :
    LaplacesEquation(Triangulation<DIM>* pMesh) 
      : mFe(1),
        mDofHandler(*pMesh)
    {
        assert(pMesh!=NULL);
        mpMesh = pMesh;

        typename Triangulation<DIM>::cell_iterator element_iter = mpMesh->begin_active();
    
        while(element_iter != mpMesh->end())
        {
            for(unsigned face_index=0; face_index<GeometryInfo<DIM>::faces_per_cell; face_index++)
            {
                if(element_iter->face(face_index)->at_boundary()) 
                {
                    element_iter->face(face_index)->set_boundary_indicator(0);
                }
            }
            element_iter++;
        }

        // distribute dofs
        mDofHandler.distribute_dofs(mFe);

        // HANGIGN NODES, SEE DEALII TUTORIAL 2
        // clear the constrait matrix
        mHangingNodeConstraints.clear();
        // form constraints 
        DoFTools::make_hanging_node_constraints(mDofHandler, mHangingNodeConstraints);
        // some postprocessing
        mHangingNodeConstraints.close();
    
        // form sparsity pattern
        mSparsityPattern.reinit(mDofHandler.n_dofs(), 
                                mDofHandler.n_dofs(),
                                mDofHandler.max_couplings_between_dofs());
    
        DoFTools::make_sparsity_pattern(mDofHandler, mSparsityPattern);
    
        // see dealii tutorial 2
        mHangingNodeConstraints.condense(mSparsityPattern);
    
    
        mSparsityPattern.compress();
        
        // initialise vectors and matrices
        mSystemMatrix.reinit(mSparsityPattern);
        mSolution.reinit(mDofHandler.n_dofs());
        mRhsVector.reinit(mDofHandler.n_dofs());
    }

    void Solve()
    {
        AssembleSystem();
        
        SolverControl solver_control(1000, 1e-12);
        PrimitiveVectorMemory<> vector_memory;
        SolverCG<>              cg(solver_control, vector_memory);
  
        cg.solve(mSystemMatrix, mSolution, mRhsVector, PreconditionIdentity());
        
        
        mSolutionAtVertices.reinit(mpMesh->n_vertices());
        
        for(unsigned i=0;i<mSolutionAtVertices.size(); i++)
        {
            mSolutionAtVertices(i) = -1e200;
        }
        
       
        DofVertexIterator<DIM> vertex_iter(mpMesh, &mDofHandler);
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();

            mSolutionAtVertices(vertex_index) = mSolution(vertex_iter.GetDof(0));

            Point<DIM> posn = vertex_iter.GetVertex();
            std::cout << vertex_index << " " << posn(0) << " " << posn(1) << " " << mSolutionAtVertices(vertex_index) << "\n" << std::flush;

            vertex_iter.Next();
        }
    }

    Vector<double> GetSolutionAtVertices()
    {
        return mSolutionAtVertices;
    }
};










/** 
 * A simple tumour source model for which s = constant at every evaluation 
 * point, where the constant is taken in in the constructor
 */
template<unsigned DIM>
class ConcentrationBasedTumourSourceModel : public AbstractGrowingTumourSourceModel<DIM>
{
private :
    friend class TestConcentrationBasedTumourSourceModel;
    
    Triangulation<DIM> mMesh;
    Vector<double> mConcentrations;

    void UpdateMesh(FiniteElasticityAssembler<DIM>* pFiniteElasticityAssembler)
    {
        
        Vector<double>& solution = pFiniteElasticityAssembler->rGetSolutionVector();
        DoFHandler<DIM>& dof_handler = pFiniteElasticityAssembler->rGetDofHandler();

        
        TriangulationVertexIterator<DIM> vertex_iter(&mMesh);        
        DofVertexIterator<DIM> fe_vertex_iter(pFiniteElasticityAssembler->GetMesh(), &dof_handler);
        
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            unsigned fe_vertex_index = fe_vertex_iter.GetVertexGlobalIndex();
            
            // check both meshes are in sync
            assert(vertex_index==fe_vertex_index);
            
            Point<DIM>& old_posn = vertex_iter.GetVertex();
            for(unsigned i=0; i<DIM; i++)
            {
                vertex_iter.GetVertex()[i] = old_posn(i)+solution(fe_vertex_iter.GetDof(i));
            } 
            vertex_iter.Next();
            fe_vertex_iter.Next();
        }
    }

public :
    ConcentrationBasedTumourSourceModel(Triangulation<DIM>& initialMesh)
    {
        mMesh.copy_triangulation(initialMesh);
    }

    void Run(double tStart, double tEnd, FiniteElasticityAssembler<DIM>* pFiniteElasticityAssembler)
    {   
        Triangulation<DIM>* p_fe_mesh = pFiniteElasticityAssembler->GetMesh();
        if(p_fe_mesh->n_vertices()==mMesh.n_vertices())
        {
            UpdateMesh(pFiniteElasticityAssembler);
        }
        else
        {
            mMesh.clear(); 
            mMesh.copy_triangulation(*p_fe_mesh);
            UpdateMesh(pFiniteElasticityAssembler);
        }
        
        LaplacesEquation<DIM> laplace(&mMesh);

        
        laplace.Solve();
        // proper solution! not solution vector.. not reference too..
        mConcentrations = laplace.GetSolutionAtVertices();

        static double death_threshold;
        static double birth_threshold;
        static bool first = true;
        
        if(first==true)
        {
            first = false;
            double max_conc = -1e10;
            double min_conc = 1e10;
            for(unsigned i=0; i<mConcentrations.size(); i++)
            {
                if(mConcentrations(i) > max_conc)
                {
                    max_conc = mConcentrations(i);
                }
                if(mConcentrations(i) < min_conc)
                {
                    min_conc = mConcentrations(i);
                }
            }
            
            death_threshold = (2.0/3.0)*min_conc + (1.0/3.0)*max_conc;
            birth_threshold = (1.0/3.0)*min_conc + (2.0/3.0)*max_conc;
            
            std::cout << "max, min = " << max_conc << " " << min_conc << "\n";
            std::cout << "d, b     = " << death_threshold << " " << birth_threshold << "\n";
        }
            

        std::cout << "\n\nSource values:\n";
        typename std::map<unsigned,EvaluationPointInfo<DIM> >::iterator iter 
           = this->mEvaluationPoints.begin();
        while(iter!=this->mEvaluationPoints.end())
        {
            unsigned mesh_index = iter->second.MeshGlobalIndex;
            double concentration = mConcentrations(mesh_index);
            
            if(concentration < death_threshold)
            {
                iter->second.SourceValue = -0.5;
            }
            else if(concentration < birth_threshold)
            {
                iter->second.SourceValue =  0.0;
            }
            else
            {
                iter->second.SourceValue =  0.5;
            }
                
            std::cout << concentration << " " << iter->second.MeshGlobalIndex << " " << iter->second.OldPosition[0] << " " << iter->second.OldPosition[1] << ": " << iter->second.SourceValue << "\n";    
                
            iter++;
        }
    }
};


#endif /*CONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_*/
