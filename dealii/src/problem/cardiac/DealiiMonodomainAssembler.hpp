#ifndef DEALIIMONODOMAINASSEMBLER_HPP_
#define DEALIIMONODOMAINASSEMBLER_HPP_

#include "AbstractDealiiAssembler.hpp"
#include "AbstractCardiacCell.hpp"
#include "SimpleDataWriter.hpp"

/**
 *  Dealii Monodomain Assembler
 * 
 *  TEMPORARY
 * 
 *  this class is a stripped down basic dealii version of MonodomainPde and MonoDg0Assembler,
 *  and is only being used as a stop gap (for coupling mechanics to monodomain) until dealii can
 *  be connected to chaste properly (ticket:442). It is not properly tested or parallel etc
 * 
 *  It will be used by CardioElectroMechanicsProblem until the [Mono/Bi]Dg0Assembler can be 
 *  plugged in instead.
 */
template<unsigned DIM>
class DealiiMonodomainAssembler : public AbstractDealiiAssembler<DIM>
{
private:
    FE_Q<DIM> mFe;

    std::vector< AbstractCardiacCell* >& mrCells;
    
    std::vector<double> mIionicCache;
    std::vector<double> mIntracellularStimulusCache;
    
    std::vector<double> mVoltage;
    
    double mDt;

    void ApplyDirichletBoundaryConditions()
    {
    }
    
    virtual void DistributeDofs()
    {
        this->mDofHandler.distribute_dofs(mFe);
    }
    
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                           Vector<double>&        elementRhs,
                           FullMatrix<double>&    elementMatrix,
                           bool                   assembleVector,
                           bool                   assembleMatrix)
    {
        double surface_area_to_volume_ratio = 1400; // 1/cm 
        double capacitance = 1.0; // uF/cm^2
        double conductivity = 1.75;
        double dudt_coeff_term =  capacitance*surface_area_to_volume_ratio;  
        
        static QGauss<DIM>   quadrature_formula(3);
        const unsigned n_q_points = quadrature_formula.n_quadrature_points;
        

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
        

        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            // interpolate source term
            double source_term=0;
            for (unsigned int j=0; j<dofs_per_element; ++j)
            {
               unsigned cell_index = local_dof_indices[j];
               double non_linear_source_term_at_node= - surface_area_to_volume_ratio*mIionicCache[cell_index]
                                                      - mIntracellularStimulusCache[cell_index];
               source_term += fe_values.shape_value(j, q_point) * non_linear_source_term_at_node;
            }
               
            double voltage=0;
            for (unsigned int j=0; j<dofs_per_element; ++j)
            {
                voltage += this->mCurrentSolution(local_dof_indices[j])*fe_values.shape_value(j,q_point); 
            }


            for (unsigned i=0; i<dofs_per_element; i++)
            {
                if(assembleMatrix)
                {
                    for (unsigned j=0; j<dofs_per_element; j++)
                    {
                        elementMatrix(i,j) +=   conductivity
                                              * fe_values.shape_grad (i, q_point) 
                                              * fe_values.shape_grad (j, q_point) 
                                              * fe_values.JxW (q_point);

                        elementMatrix(i,j) +=   fe_values.shape_value(i, q_point) 
                                              * fe_values.shape_value(j, q_point) 
                                              * fe_values.JxW (q_point) 
                                              * dudt_coeff_term / mDt;
                    }
                }
                if(assembleVector)
                {
                    elementRhs(i) +=   fe_values.shape_value(i, q_point) 
                                     * fe_values.JxW (q_point) 
                                     * (voltage * dudt_coeff_term / mDt + source_term);
                }
            }
        }
    }


public:
    DealiiMonodomainAssembler(Triangulation<DIM>* pMesh, std::vector<AbstractCardiacCell*>& rCells)
        : AbstractDealiiAssembler<DIM>(pMesh),
          mrCells(rCells),
          mFe(1)
    {
        assert(rCells.size() == pMesh->n_vertices());

        double initial_voltage = -83;
            
        mIionicCache.resize(rCells.size(), 0.0);
        mIntracellularStimulusCache.resize(rCells.size(), 0.0);
        mVoltage.resize(rCells.size(), initial_voltage);
        
        this->mDofsPerElement = mFe.dofs_per_cell;
        DistributeDofs();
        this->InitialiseMatricesVectorsAndConstraints();
        
        for(unsigned i=0; i<this->mCurrentSolution.size(); i++)
        {
            this->mCurrentSolution(i) = initial_voltage;
        }
    
    }
    
    void Solve(double startTime, double endTime, unsigned numTimeSteps)
    {
        mDt = (endTime-startTime)/numTimeSteps;
 
        SimpleDataWriter writer("DealiiMonodomain", "voltage_0.dat", mVoltage);
        
        std::vector<double> x(this->mpMesh->n_vertices());
        std::vector<double> y(this->mpMesh->n_vertices());
        TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
        while(!vertex_iter.ReachedEnd())
        {
            x[vertex_iter.GetVertexGlobalIndex()] = vertex_iter.GetVertex()[0]; 
            y[vertex_iter.GetVertexGlobalIndex()] = vertex_iter.GetVertex()[1]; 
            
            vertex_iter.Next();
        }
        SimpleDataWriter writer2("DealiiMonodomain", "mesh.dat", x, y, false);
        
   
        bool assemble_matrix = true; 
        for(unsigned time_counter = 0; time_counter < numTimeSteps; time_counter++)
        {
            double time = startTime + time_counter*mDt;
        
            // solve cells
            for(unsigned cell_index=0; cell_index<mrCells.size(); cell_index++)
            {
                // set voltage
                mrCells[cell_index]->SetVoltage( mVoltage[cell_index] );
                
                // solve
                mrCells[cell_index]->ComputeExceptVoltage(time, time+mDt);
                
                // update caches
                mIionicCache[cell_index] = mrCells[cell_index]->GetIIonic();
                mIntracellularStimulusCache[cell_index] = mrCells[cell_index]->GetIntracellularStimulus(time+mDt);
            }
            
            // assemble system
            this->AssembleSystem(true, assemble_matrix);
            assemble_matrix = false;
            
            // solve for voltage
            SolverControl solver_control(1000, 1e-12, false, false);
            PrimitiveVectorMemory<> vector_memory;
            SolverCG<> cg(solver_control, vector_memory);
        
            cg.solve(this->mSystemMatrix, this->mCurrentSolution, this->mRhsVector, PreconditionIdentity());
            
            // set up voltage store (mCurrentSolution is indexed by dof index not vertex index)
            DofVertexIterator<DIM> vertex_iter(this->mpMesh, &this->mDofHandler);
            while (!vertex_iter.ReachedEnd())
            {
                unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
                mVoltage[vertex_index] = this->mCurrentSolution(vertex_iter.GetDof(0));
                vertex_iter.Next();
            }
            
            std::stringstream file;
            file << "voltage_" << time_counter+1 << ".dat";
            SimpleDataWriter writer("DealiiMonodomain", file.str(), mVoltage, false);
        }
    }
    
    std::vector<double>& GetVoltage()
    {
        return mVoltage;
    }
};

#endif /*DEALIIMONODOMAINASSEMBLER_HPP_*/
