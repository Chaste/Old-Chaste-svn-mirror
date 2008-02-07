#ifndef _TISSUESIMULATIONWITHNUTRIENTS_CPP_
#define _TISSUESIMULATIONWITHNUTRIENTS_CPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "TissueSimulationWithNutrients.hpp"
#include "SimpleLinearEllipticAssembler.hpp"

//////////////////////////////////////////////////////////////////////////////
//                          Setup/AfterSolve methods                        //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetupSolve()
{
    if (this->mrTissue.Begin() != this->mrTissue.End())
    {
        SetupWriteNutrient();
        double current_time = SimulationTime::Instance()->GetDimensionalisedTime();            
        WriteNutrient(current_time);                      
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetupWriteNutrient()
{ 
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/vis_results/",false);
    if (output_file_handler.IsMaster())
    {
        mpNutrientResultsFile = output_file_handler.OpenOutputFile("results.viznutrient");
        *this->mpSetupFile << "Nutrient \n" ;
        if (mWriteAverageRadialNutrientResults)
        {
            mpAverageRadialNutrientResultsFile = output_file_handler.OpenOutputFile("radial_dist.dat");                     
        }
    }        
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::UseCoarseNutrientMesh(double coarseGrainScaleFactor)
{
    assert(mpAveragedSinksPde);
    CreateCoarseNutrientMesh(coarseGrainScaleFactor);
}    

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::CreateCoarseNutrientMesh(double coarseGrainScaleFactor)
{
    // Set up coarse nutrient mesh
//    \todo: we could instead use the disk with 984 elements etc.
//    \todo: translate the mesh so it's centre is over the centre of the 
//           cells - could then also automatically calculate the scale 
//           factor from the initial dimensions of the cells and the end time
    TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
    mpCoarseNutrientMesh = new ConformingTetrahedralMesh<2,2>;
    mpCoarseNutrientMesh->ConstructFromMeshReader(mesh_reader);
    
    // Scale mesh
    mpCoarseNutrientMesh->Scale(coarseGrainScaleFactor, coarseGrainScaleFactor);
}    

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::AfterSolve()
{
    if (this->mrTissue.Begin() != this->mrTissue.End() // if there are any cells
    && PetscTools::AmMaster())
    {
        mpNutrientResultsFile->close();
        
        if (mWriteAverageRadialNutrientResults)
        {
            WriteAverageRadialNutrientDistribution(SimulationTime::Instance()->GetDimensionalisedTime(), mNumRadialIntervals);
            mpAverageRadialNutrientResultsFile->close();
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//                             PostSolve methods                            // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SolveNutrientPde()
{
    if(mpCoarseNutrientMesh!=NULL)
    {
        SolveNutrientPdeUsingCoarseMesh();
        return ;
    }
    assert(mpAveragedSinksPde == NULL);
    assert(mpPde);
    ConformingTetrahedralMesh<DIM,DIM>& r_mesh = this->mrTissue.rGetMesh();
    CellwiseData<DIM>::Instance()->ReallocateMemory();
    std::set<unsigned> ghost_node_indices = this->mrTissue.GetGhostNodeIndices();
    
    // We shouldn't have any ghost nodes in a TissueSimulationWithNutrients
    assert(ghost_node_indices.size()==0);
            
    // Set up boundary conditions
    BoundaryConditionsContainer<DIM,DIM,1> bcc;
    ConstBoundaryCondition<DIM>* p_boundary_condition = new ConstBoundaryCondition<DIM>(1.0);
    for (typename ConformingTetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = r_mesh.GetBoundaryNodeIteratorBegin();
         node_iter != r_mesh.GetBoundaryNodeIteratorEnd();
         ++node_iter)
    {
        bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
    }
    
    // Set up assembler - note this is a purpose-made elliptic assembler
    // that interpolates the source terms from node onto gauss points,
    // as for a nutrients simulation the source will only be known at the
    // cells (nodes), not the gauss points
    TissueSimulationWithNutrientsAssembler<DIM> assembler(&r_mesh,mpPde,&bcc);
    
    PetscInt size_of_soln_previous_step = 0;
    
    if(mNutrientSolution)
    {
        VecGetSize(mNutrientSolution, &size_of_soln_previous_step);
    }
    
    if (size_of_soln_previous_step == (int)r_mesh.GetNumNodes())
    {
        // We make an initial guess which gets copied by the Solve method of
        // SimpleLinearSolver, so we need to delete it too.
        Vec initial_guess;
        VecDuplicate(mNutrientSolution, &initial_guess);
        VecCopy(mNutrientSolution, initial_guess);
        
        // Use current solution as the initial guess
        VecDestroy(mNutrientSolution);    // Solve method makes its own mNutrientSolution
        mNutrientSolution = assembler.Solve(initial_guess);
        VecDestroy(initial_guess);
    }
    else
    {
        if (mNutrientSolution)
        {
            VecDestroy(mNutrientSolution);
        }
        mNutrientSolution = assembler.Solve();
    }            

    ReplicatableVector result_repl(mNutrientSolution);

    // Update cellwise data
    for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
    {
        double oxygen_conc = result_repl[i];
        CellwiseData<DIM>::Instance()->SetValue(oxygen_conc, r_mesh.GetNode(i));
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SolveNutrientPdeUsingCoarseMesh()
{
    assert(mpPde==NULL);
    assert(mpAveragedSinksPde);
    ConformingTetrahedralMesh<DIM,DIM>& r_mesh = *mpCoarseNutrientMesh;
    CellwiseData<DIM>::Instance()->ReallocateMemory();
    // We shouldn't have any ghost nodes in a TissueSimulationWithNutrients
    std::set<unsigned> ghost_node_indices = this->mrTissue.GetGhostNodeIndices();
    assert(ghost_node_indices.size()==0);
    
    // Loop over cells and calculate centre of distribution
    c_vector<double, DIM> centre;
    for(typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        centre += cell_iter.rGetLocation();
    }
    centre /= this->mrTissue.GetNumRealCells();
    
    // Find max radius
    double max_radius =0.0;
    for(typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        double radius = norm_2(centre - cell_iter.rGetLocation() ) ;
        if (radius > max_radius)
        {
            max_radius = radius;
        }
    }    
    // Set up boundary conditions
    BoundaryConditionsContainer<DIM,DIM,1> bcc;
    ConstBoundaryCondition<DIM>* p_boundary_condition = new ConstBoundaryCondition<DIM>(1.0);
    for (unsigned i=0; i< r_mesh.GetNumNodes(); i++)
    {
        double distance_from_centre = norm_2(r_mesh.GetNode(i)->rGetLocation() - centre);
        if( distance_from_centre > max_radius )
        {
            bcc.AddDirichletBoundaryCondition(r_mesh.GetNode(i), p_boundary_condition,0,false);
        }
    }
    
//    for (typename ConformingTetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = r_mesh.GetBoundaryNodeIteratorBegin();
//         node_iter != r_mesh.GetBoundaryNodeIteratorEnd();
//         ++node_iter)
//    {
//        bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
//    }
    
    PetscInt size_of_soln_previous_step = 0;
    
    if(mNutrientSolution)
    {
        VecGetSize(mNutrientSolution, &size_of_soln_previous_step);
    }
    
    mpAveragedSinksPde->SetupSourceTerms(*mpCoarseNutrientMesh);
    
    SimpleLinearEllipticAssembler<DIM,DIM> assembler(mpCoarseNutrientMesh, mpAveragedSinksPde, &bcc);
    if (size_of_soln_previous_step == (int)r_mesh.GetNumNodes())
    {
        // We make an initial guess which gets copied by the Solve method of
        // SimpleLinearSolver, so we need to delete it too.
        Vec initial_guess;
        VecDuplicate(mNutrientSolution, &initial_guess);
        VecCopy(mNutrientSolution, initial_guess);
        
        // Use current solution as the initial guess
        VecDestroy(mNutrientSolution);    // Solve method makes its own mNutrientSolution
        mNutrientSolution = assembler.Solve(initial_guess);
        VecDestroy(initial_guess);
    }
    else
    {
        if (mNutrientSolution)
        {
            VecDestroy(mNutrientSolution);
        }
        mNutrientSolution = assembler.Solve();
    }            

    // Update cellwise data - since the cells are not nodes on the coarse
    // mesh we have to interpolate from the nodes of the coarse mesh onto
    // the cell locations
    ReplicatableVector nutrient_repl(mNutrientSolution);

    for(typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = cell_iter.rGetLocation();
        unsigned elem_index = mpCoarseNutrientMesh->GetContainingElementIndex(r_position_of_cell);

        Element<DIM,DIM>* p_element = mpCoarseNutrientMesh->GetElement(elem_index);
    
        c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(r_position_of_cell);
         
        double interpolated_nutrient = 0.0;
        for(unsigned i=0; i<DIM+1/*num_nodes*/; i++)
        {         
            double nodal_value = nutrient_repl[ p_element->GetNodeGlobalIndex(i) ];
            interpolated_nutrient += nodal_value*weights(i);
        }
    
        CellwiseData<DIM>::Instance()->SetValue( interpolated_nutrient, cell_iter.GetNode() );
    }
}


template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::PostSolve()
{
    SolveNutrientPde();
    
    // Save results to file
    SimulationTime* p_time = SimulationTime::Instance();
    
    double time_next_step = p_time->GetDimensionalisedTime() + p_time->GetTimeStep();
    
    if ((p_time->GetTimeStepsElapsed()+1)%this->mSamplingTimestepMultiple == 0)
    {
        WriteNutrient(time_next_step);
    }
        
#define COVERAGE_IGNORE
    // Note: The number of timesteps per day is equal to 2880=24*120
    if ( mWriteDailyAverageRadialNutrientResults &&
         (p_time->GetTimeStepsElapsed()+1)%2880==0 )
    {
        WriteAverageRadialNutrientDistribution(time_next_step, mNumRadialIntervals);
    }
#undef COVERAGE_IGNORE
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::WriteNutrient(double time)
{
    if (PetscTools::AmMaster())
    {
        // Since there are no ghost nodes, the number of nodes must equal the number of real cells 
        assert(this->mrTissue.rGetMesh().GetNumNodes()==this->mrTissue.GetNumRealCells());
        
        (*mpNutrientResultsFile) << time << "\t";
        
        unsigned global_index; 
        double x;
        double y;
        double nutrient;

        for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
             cell_iter != this->mrTissue.End();
             ++cell_iter)
        {
            global_index = cell_iter.GetNode()->GetIndex();
            x = cell_iter.rGetLocation()[0];
            y = cell_iter.rGetLocation()[1];                    
            nutrient = CellwiseData<DIM>::Instance()->GetValue(&(*cell_iter));                    
            
            (*mpNutrientResultsFile) << global_index << " " << x << " " << y << " " << nutrient << " ";
        }            
        (*mpNutrientResultsFile) << "\n";
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetWriteAverageRadialNutrientResults(unsigned numRadialIntervals, bool writeDailyResults)
{
    mWriteAverageRadialNutrientResults = true;
    mNumRadialIntervals = numRadialIntervals;
    mWriteDailyAverageRadialNutrientResults = writeDailyResults;
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::WriteAverageRadialNutrientDistribution(double time, unsigned numRadialIntervals)
{
    (*mpAverageRadialNutrientResultsFile) << time << " "; 
    
    // Get reference to the mesh and its size
    ConformingTetrahedralMesh<DIM,DIM>& r_mesh = this->mrTissue.rGetMesh();
    unsigned num_nodes = r_mesh.GetNumNodes();
    
    // Calculate the centre of the tissue
    c_vector<double,DIM> centre = zero_vector<double>(DIM);
    for (unsigned i=0; i< num_nodes; i++)
    {
        centre += r_mesh.GetNode(i)->rGetLocation();    
    }
    centre /= (double) num_nodes;
    
    // Calculate the distance between each node and the centre 
    // of the tissue, as well as the maximum of these
    std::map<double, TissueCell*> distance_cell_map;
    
    double max_distance_from_centre = 0.0;
            
    for (unsigned i=0; i<this->mrTissue.GetNumRealCells(); i++)
    {
        double distance = norm_2(r_mesh.GetNode(i)->rGetLocation()-centre);            
        distance_cell_map[distance] = &(this->mrTissue.rGetCellAtNodeIndex(i));
        
        if (distance > max_distance_from_centre)
        {
            max_distance_from_centre = distance;
        }
    }

    // Create vector of radius intervals                                  
    std::vector<double> radius_intervals;
    for (unsigned i=0; i<numRadialIntervals; i++)
    {
        double upper_radius = max_distance_from_centre*((double) i+1)/((double) numRadialIntervals);
        radius_intervals.push_back(upper_radius);       
    }
    
    // Calculate nutrient concentration in each radial interval
    double lower_radius = 0.0;
    for (unsigned i=0; i<numRadialIntervals; i++)
    {   
        unsigned counter = 0;
        double average_conc = 0.0;
        
        for (std::map<double, TissueCell*>::iterator iter=distance_cell_map.begin();
             iter != distance_cell_map.end();
             ++iter)
        {
            if ((*iter).first > lower_radius && (*iter).first <= radius_intervals[i])
            {
                average_conc += CellwiseData<DIM>::Instance()->GetValue((*iter).second); 
                counter++;
            }
        }
        if (counter > 0)
        {
            average_conc /= (double) counter;
        }
        
        // Write results to file
        (*mpAverageRadialNutrientResultsFile) << radius_intervals[i] << " " << average_conc << " ";
        lower_radius = radius_intervals[i];
    }
    (*mpAverageRadialNutrientResultsFile) << "\n";
}

//////////////////////////////////////////////////////////////////////////////
//                           Save/Load methods                              // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::Save()
{
    CommonSave(this);
}

template<unsigned DIM>
TissueSimulationWithNutrients<DIM>* TissueSimulationWithNutrients<DIM>::Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
{
    std::string archive_filename =
        TissueSimulation<DIM>::GetArchivePathname(rArchiveDirectory, rTimeStamp);

    // Create an input archive
    std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
    boost::archive::text_iarchive input_arch(ifs);

    TissueSimulation<DIM>::CommonLoad(input_arch);

    TissueSimulationWithNutrients<DIM>* p_sim; 
    input_arch >> p_sim;
     
    if (p_sim->rGetTissue().rGetMesh().GetNumNodes()!=p_sim->rGetTissue().rGetCells().size()) 
    { 
        #define COVERAGE_IGNORE 
        std::stringstream string_stream; 
        string_stream << "Error in Load(), number of nodes (" << p_sim->rGetTissue().rGetMesh().GetNumNodes() 
                      << ") is not equal to the number of cells (" << p_sim->rGetTissue().rGetCells().size()  
                      << ")"; 
        EXCEPTION(string_stream.str()); 
        #undef COVERAGE_IGNORE 
    } 
      
    return p_sim;         
}       


#endif //_TISSUESIMULATIONWITHNUTRIENTS_CPP_
