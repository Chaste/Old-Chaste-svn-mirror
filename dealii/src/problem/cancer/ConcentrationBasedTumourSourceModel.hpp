/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_
#define CONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_

#include "AbstractGrowingTumourSourceModel.hpp"
#include "TriangulationVertexIterator.hpp"
#include "AbstractDealiiAssembler.hpp"

// todo: lots. move laplaces eqn class to somewhere sensible (refactor assemblers
// this, test update mesh, test copy mesh, test source terms.

#define COVERAGE_IGNORE

template<unsigned DIM>
class LaplacesEquation : public AbstractDealiiAssembler<DIM>
{
protected:
    FE_Q<DIM>        mFe;    
    Vector<double>   mSolutionAtVertices;
    
    void DistributeDofs()
    {
        this->mDofHandler.distribute_dofs(mFe);
    }
    
    virtual void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                                   Vector<double>&                                 elementRhs,
                                   FullMatrix<double>&                             elementMatrix,
                                   bool                                            assembleVector,
                                   bool                                            assembleMatrix)
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
        
        Point<1> conc;
        
        double source_term = 1;
        
//// !!!!!!!!!!!!!!!! GROWING_REGION == 99
//        if (elementIter->material_id()==99)
//        {
//            source_term = 1;
//        }
//        

        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            for (unsigned i=0; i<dofs_per_element; i++)
            {
                if(assembleMatrix)
                {
                    for (unsigned j=0; j<dofs_per_element; j++)
                    {
                        elementMatrix(i,j) +=   fe_values.shape_grad(i,q_point)
                                              * fe_values.shape_grad(j,q_point)
                                              * fe_values.JxW(q_point);
                    }
                }
                if(assembleVector)
                {
                    elementRhs(i) +=    fe_values.shape_value(i,q_point)
                                      * source_term
                                      * fe_values.JxW (q_point);
                }
            }
        }
    }

    void ApplyDirichletBoundaryConditions()
    {
        std::map<unsigned int,double> boundary_values;
        VectorTools::interpolate_boundary_values(this->mDofHandler,
                                                 0,
                                                 ZeroFunction<2>(),
                                                 boundary_values);
                                                 
        MatrixTools::apply_boundary_values(boundary_values,
                                           this->mSystemMatrix,
                                           this->mCurrentSolution,
                                           this->mRhsVector);
    }
    
public :
    LaplacesEquation(Triangulation<DIM>* pMesh)
            :   AbstractDealiiAssembler<DIM>(pMesh),
                mFe(1)
    {  
        // distribute dofs
        this->mDofHandler.distribute_dofs(mFe);
        this->InitialiseMatricesVectorsAndConstraints();
        
        this->mDofsPerElement = mFe.dofs_per_cell;
        
        assert(pMesh!=NULL);

        typename Triangulation<DIM>::cell_iterator element_iter = this->mpMesh->begin_active();
        
        while (element_iter != this->mpMesh->end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<DIM>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    element_iter->face(face_index)->set_boundary_indicator(0);
                }
            }
            element_iter++;
        }
    }
    
    void Solve()
    {
        this->AssembleSystem(true, true);
        
        SolverControl solver_control(1000, 1e-12);
        PrimitiveVectorMemory<> vector_memory;
        SolverCG<>              cg(solver_control, vector_memory);
        
        cg.solve(this->mSystemMatrix, this->mCurrentSolution, this->mRhsVector, PreconditionIdentity());
        
        
        mSolutionAtVertices.reinit(this->mpMesh->n_vertices());
        for (unsigned i=0;i<mSolutionAtVertices.size(); i++)
        {
            mSolutionAtVertices(i) = -1e200;
        }
        
        
        DofVertexIterator<DIM> vertex_iter(this->mpMesh, &this->mDofHandler);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            
            mSolutionAtVertices(vertex_index) = mCurrentSolution(vertex_iter.GetDof(0));
            
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
    
    Triangulation<DIM> mDeformedMesh;
    Vector<double> mConcentrations;
    

    void UpdateMesh(FiniteElasticityAssembler<DIM>* pFiniteElasticityAssembler)
    {
        std::vector<Vector<double> >& deformed_position = pFiniteElasticityAssembler->rGetDeformedPosition();
        
        assert(deformed_position[0].size()==mDeformedMesh.n_vertices());

        TriangulationVertexIterator<DIM> vertex_iter(&mDeformedMesh);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            for (unsigned i=0; i<DIM; i++)
            {
                vertex_iter.GetVertex()[i] = deformed_position[i](vertex_index);
            }
            vertex_iter.Next();
        }
    }

    
public :
    ConcentrationBasedTumourSourceModel(Triangulation<DIM>& mesh)
    {
        mDeformedMesh.copy_triangulation(mesh);
    }
    
    void Run(double tStart, double tEnd, FiniteElasticityAssembler<DIM>* pFiniteElasticityAssembler)
    {
        Triangulation<DIM>* p_fe_mesh = pFiniteElasticityAssembler->GetMesh();
        if (p_fe_mesh->n_vertices()==mDeformedMesh.n_vertices())
        {
            UpdateMesh(pFiniteElasticityAssembler);
        }
        else
        {
            mDeformedMesh.clear();
            mDeformedMesh.copy_triangulation(*p_fe_mesh);
            UpdateMesh(pFiniteElasticityAssembler);
        }
        
        LaplacesEquation<DIM> laplace(&mDeformedMesh);
        
        laplace.Solve();

        // proper solution, not solution vector.. not reference too..
        mConcentrations = laplace.GetSolutionAtVertices();
        
        static double death_threshold;
        static double birth_threshold;
        static bool first = true;
        
        if (first==true)
        {
            first = false;
            double max_conc = -1e10;
            double min_conc = 1e10;
            for (unsigned i=0; i<mConcentrations.size(); i++)
            {
                if (mConcentrations(i) > max_conc)
                {
                    max_conc = mConcentrations(i);
                }
                if (mConcentrations(i) < min_conc)
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
        while (iter!=this->mEvaluationPoints.end())
        {
            unsigned mesh_index = iter->first;

            double concentration = mConcentrations(mesh_index);
            
            if (concentration < death_threshold)
            {
                iter->second.SourceValue = -0.5;
            }
            else if (concentration < birth_threshold)
            {
                iter->second.SourceValue =  0.0;
            }
            else
            {
                iter->second.SourceValue =  0.5;
            }
            
            std::cout << concentration << " " << iter->first << " " << iter->second.OldPosition[0] << " " << iter->second.OldPosition[1] << ": " << iter->second.SourceValue << "\n";
            
            iter++;
        }
    }
};

#undef COVERAGE_IGNORE


#endif /*CONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_*/
