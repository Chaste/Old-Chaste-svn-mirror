/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#include "TissueSimulationWithPdes.hpp"
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "SimpleDataWriter.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "TissueSimulationWithPdesAssembler.hpp"
#include "CellwiseData.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"
#include "TrianglesMeshReader.hpp"

template<unsigned DIM>
TissueSimulationWithPdes<DIM>::TissueSimulationWithPdes(AbstractTissue<DIM>& rTissue,
					                                    std::vector<AbstractForce<DIM>*> forceCollection,
					                                    AbstractLinearEllipticPde<DIM,DIM>* pPde,
							                            double boundaryValue,
							                            bool isNeumannBoundaryCondition,
					                                    bool deleteTissueAndForceCollection,
					                                    bool initialiseCells)
    : TissueSimulation<DIM>(rTissue,
                            forceCollection,
                            deleteTissueAndForceCollection,
                            initialiseCells),
      mCurrentPdeSolution(NULL),
      mpPde(pPde),
      mWriteAverageRadialPdeSolution(false),
      mWriteDailyAverageRadialPdeSolution(false),
      mNumRadialIntervals(0), // 'unset' value
      mpCoarsePdeMesh(NULL),
      mBoundaryValue(boundaryValue),
      mIsNeumannBoundaryCondition(isNeumannBoundaryCondition)
{
    // We must be using a mesh-based tissue
    assert(dynamic_cast<MeshBasedTissue<DIM>*>(&(this->mrTissue)) != NULL);

    // We must not have any ghost nodes
    assert(dynamic_cast<MeshBasedTissueWithGhostNodes<DIM>*>(&(this->mrTissue)) == NULL);
}

template<unsigned DIM>
TissueSimulationWithPdes<DIM>::~TissueSimulationWithPdes()
{
    if (mCurrentPdeSolution)
    {
        VecDestroy(mCurrentPdeSolution);
    }
    if (mpCoarsePdeMesh)
    {
        delete mpCoarsePdeMesh;
    }
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::SetPde(AbstractLinearEllipticPde<DIM,DIM>* pPde)
{
    mpPde = pPde;
}

template<unsigned DIM>
Vec TissueSimulationWithPdes<DIM>::GetCurrentPdeSolution()
{
    return mCurrentPdeSolution;
}

//////////////////////////////////////////////////////////////////////////////
//                          Setup/AfterSolve methods                        //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::WriteVisualizerSetupFile()
{
    for (unsigned i=0; i<this->mForceCollection.size(); i++)
    {
        if (dynamic_cast<AbstractTwoBodyInteractionForce<DIM>*>(this->mForceCollection[i]))
        {
            double cutoff = (static_cast<AbstractTwoBodyInteractionForce<DIM>*>(this->mForceCollection[i]))->GetCutoffPoint();
            *(this->mpVizSetupFile) << "Cutoff\t" << cutoff << "\n";
        }
    }
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::SetupSolve()
{
    if (mpCoarsePdeMesh != NULL)
    {
        InitialiseCoarsePdeMesh();
    }

    // We must initially have at least one cell in the tissue simulation
    assert(this->mrTissue.GetNumRealCells() != 0);

    SetupWritePdeSolution();
    double current_time = SimulationTime::Instance()->GetTime();
    WritePdeSolution(current_time);
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::SetupWritePdeSolution()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/", false);
    if (PetscTools::AmMaster())
    {
        mpVizPdeSolutionResultsFile = output_file_handler.OpenOutputFile("results.vizpdesolution");
        *this->mpVizSetupFile << "PDE \n";
        if (mWriteAverageRadialPdeSolution)
        {
            mpAverageRadialPdeSolutionResultsFile = output_file_handler.OpenOutputFile("radial_dist.dat");
        }
    }
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::UseCoarsePdeMesh(double coarseGrainScaleFactor)
{
    assert(dynamic_cast<AveragedSinksPde<DIM>*>(mpPde) != NULL);
    CreateCoarsePdeMesh(coarseGrainScaleFactor);
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::CreateCoarsePdeMesh(double coarseGrainScaleFactor)
{
    EXCEPTION("This method is only implemented in 2D");
}

/**
 * The CreateCoarsePdeMesh method is currently only implemented in 2D, hence there
 * are two definitions to this method (one templated and one not).
 *
 * @param coarseGrainScaleFactor the ratio of the width of the coarse PDE mesh to the initial width of the tissue
 */
template<>
void TissueSimulationWithPdes<2>::CreateCoarsePdeMesh(double coarseGrainScaleFactor)
{
    // Create coarse PDE mesh (can use a larger mesh if required, e.g. disk_984_elements)
    TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
    mpCoarsePdeMesh = new TetrahedralMesh<2,2>;
    mpCoarsePdeMesh->ConstructFromMeshReader(mesh_reader);

    // Find centre of tissue
    c_vector<double,2> centre_of_tissue = zero_vector<double>(2);
    for (AbstractTissue<2>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        centre_of_tissue += this->mrTissue.GetLocationOfCellCentre(*cell_iter);
    }
    centre_of_tissue /= this->mrTissue.GetNumRealCells();

    // Find max radius of tissue
    double max_tissue_radius = 0.0;
    for (AbstractTissue<2>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        double radius = norm_2(centre_of_tissue - this->mrTissue.GetLocationOfCellCentre(*cell_iter));
        if (radius > max_tissue_radius)
        {
            max_tissue_radius = radius;
        }
    }

    // Find centre of coarse PDE mesh
    c_vector<double,2> centre_of_coarse_mesh = zero_vector<double>(2);
    for (unsigned i=0; i<mpCoarsePdeMesh->GetNumNodes(); i++)
    {
        centre_of_coarse_mesh += mpCoarsePdeMesh->GetNode(i)->rGetLocation();
    }
    centre_of_coarse_mesh /= mpCoarsePdeMesh->GetNumNodes();

    // Find max radius of coarse PDE mesh
    double max_mesh_radius = 0.0;
    for (unsigned i=0; i<mpCoarsePdeMesh->GetNumNodes(); i++)
    {
        double radius = norm_2(centre_of_coarse_mesh - mpCoarsePdeMesh->GetNode(i)->rGetLocation());
        if (radius > max_mesh_radius)
        {
            max_mesh_radius = radius;
        }
    }

    // Translate centre of coarse PDE mesh to the origin
    mpCoarsePdeMesh->Translate(-centre_of_coarse_mesh[0], -centre_of_coarse_mesh[1]);

    // Scale coarse PDE mesh
    double scale_factor = (max_tissue_radius/max_mesh_radius)*coarseGrainScaleFactor;
    mpCoarsePdeMesh->Scale(scale_factor, scale_factor);

    // Translate centre of coarse PDE mesh to centre of the tissue
    mpCoarsePdeMesh->Translate(centre_of_tissue[0], centre_of_tissue[1]);
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::InitialiseCoarsePdeMesh()
{
    mCellPdeElementMap.clear();

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        // Find the element of mpCoarsePdeMesh that contains this cell
        const ChastePoint<DIM>& r_position_of_cell = this->mrTissue.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = mpCoarsePdeMesh->GetContainingElementIndex(r_position_of_cell);
        mCellPdeElementMap[&(*cell_iter)] = elem_index;
    }
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::AfterSolve()
{
    if (this->mrTissue.GetNumRealCells() != 0 && PetscTools::AmMaster())
    {
        mpVizPdeSolutionResultsFile->close();

        if (mWriteAverageRadialPdeSolution)
        {
            WriteAverageRadialPdeSolution(SimulationTime::Instance()->GetTime(), mNumRadialIntervals);
            mpAverageRadialPdeSolutionResultsFile->close();
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//                             PostSolve methods                            //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::SolvePde()
{
    if (mpCoarsePdeMesh != NULL)
    {
        SolvePdeUsingCoarseMesh();
        return;
    }

    assert(mpPde);
    assert(dynamic_cast<AveragedSinksPde<DIM>*>(mpPde) == NULL);

    // Note: If not using a coarse PDE mesh, we MUST be using a MeshBasedTissue
    // Make sure the mesh is in a nice state
    this->mrTissue.Update();

    TetrahedralMesh<DIM,DIM>& r_mesh = static_cast<MeshBasedTissue<DIM>*>(&(this->mrTissue))->rGetMesh();
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    // Set up boundary conditions
    BoundaryConditionsContainer<DIM,DIM,1> bcc;
    ConstBoundaryCondition<DIM>* p_bc = new ConstBoundaryCondition<DIM>(mBoundaryValue);
    if (mIsNeumannBoundaryCondition) ///\todo test this! (#1465)
    {
    	for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = r_mesh.GetBoundaryElementIteratorBegin();
    	     elem_iter != r_mesh.GetBoundaryElementIteratorEnd();
	         ++elem_iter)
	    {
	        bcc.AddNeumannBoundaryCondition(*elem_iter, p_bc);
	    }
    }
    else
    {
    	for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = r_mesh.GetBoundaryNodeIteratorBegin();
	         node_iter != r_mesh.GetBoundaryNodeIteratorEnd();
	         ++node_iter)
	    {
            bcc.AddDirichletBoundaryCondition(*node_iter, p_bc);
         }
    }

    /*
     * Set up assembler. This is a purpose-made elliptic assembler which must
     * interpolate contributions to source terms from nodes onto Gauss points,
     * because the PDE solution is only stored at the cells (nodes).
     */
    TissueSimulationWithPdesAssembler<DIM> assembler(&r_mesh, mpPde, &bcc);

    PetscInt size_of_soln_previous_step = 0;

    if (mCurrentPdeSolution)
    {
        VecGetSize(mCurrentPdeSolution, &size_of_soln_previous_step);
    }
    if (size_of_soln_previous_step == (int)r_mesh.GetNumNodes())
    {
        // We make an initial guess which gets copied by the Solve method of
        // SimpleLinearSolver, so we need to delete it too.
        Vec initial_guess;
        VecDuplicate(mCurrentPdeSolution, &initial_guess);
        VecCopy(mCurrentPdeSolution, initial_guess);

        // Use current solution as the initial guess
        VecDestroy(mCurrentPdeSolution);    // Solve method makes its own mCurrentPdeSolution
        mCurrentPdeSolution = assembler.Solve(initial_guess);
        VecDestroy(initial_guess);
    }
    else
    {
        if (mCurrentPdeSolution)
        {
            assert(size_of_soln_previous_step != 0);
            VecDestroy(mCurrentPdeSolution);
        }
        mCurrentPdeSolution = assembler.Solve();
    }

    ReplicatableVector solution_repl(mCurrentPdeSolution);

    // Update cellwise data
    for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
    {
        double solution = solution_repl[i];
        unsigned index = r_mesh.GetNode(i)->GetIndex();
        CellwiseData<DIM>::Instance()->SetValue(solution, index);
    }
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::SolvePdeUsingCoarseMesh()
{
    assert(dynamic_cast<AveragedSinksPde<DIM>*>(mpPde) != NULL);

    TetrahedralMesh<DIM,DIM>& r_mesh = *mpCoarsePdeMesh;
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    // Loop over cells and calculate centre of distribution
    c_vector<double, DIM> centre = zero_vector<double>(DIM);
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        centre += this->mrTissue.GetLocationOfCellCentre(*cell_iter);
    }
    centre /= this->mrTissue.GetNumRealCells();

    // Find max radius
    double max_radius = 0.0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        double radius = norm_2(centre - this->mrTissue.GetLocationOfCellCentre(*cell_iter));
        if (radius > max_radius)
        {
            max_radius = radius;
        }
    }

    // Set up boundary conditions
    BoundaryConditionsContainer<DIM,DIM,1> bcc;
    ConstBoundaryCondition<DIM>* p_bc = new ConstBoundaryCondition<DIM>(mBoundaryValue);
    if (mIsNeumannBoundaryCondition) ///\todo test this! (#1465)
    {
    	EXCEPTION("Neumann BCs not yet implemented when using a coarse PDE mesh");
    }
    else
    {
	    // Get the set of coarse element indices that contain tissue cells
	    std::set<unsigned> coarse_element_indices_in_map;
	    for (typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
	        cell_iter != this->mrTissue.End();
	        ++cell_iter)
	    {
	        coarse_element_indices_in_map.insert(mCellPdeElementMap[&(*cell_iter)]);
	    }
	
	    // Find the node indices that associated with elements whose
	    // indices are NOT in the set coarse_element_indices_in_map
	    std::set<unsigned> coarse_mesh_boundary_node_indices;
	
	    for (unsigned i=0; i<r_mesh.GetNumElements(); i++)
	    {
	        // If the element index is NOT in the set...
	        if (coarse_element_indices_in_map.find(i) == coarse_element_indices_in_map.end())
	        {
	            // ... then get the element...
	            Element<DIM,DIM>* p_element = r_mesh.GetElement(i);
	
	            // ... and add its associated nodes to coarse_mesh_boundary_node_indices
	            for (unsigned local_index=0; local_index<DIM+1; local_index++)
	            {
	                unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
	                coarse_mesh_boundary_node_indices.insert(node_index);
	            }
	        }
	    }
	
	    // Apply boundary condition to the nodes in the set coarse_mesh_boundary_node_indices
	    for (std::set<unsigned>::iterator iter = coarse_mesh_boundary_node_indices.begin();
	         iter != coarse_mesh_boundary_node_indices.end();
	         ++iter)
	    {
	        bcc.AddDirichletBoundaryCondition(r_mesh.GetNode(*iter), p_bc, 0, false);
	    }
    }

    PetscInt size_of_soln_previous_step = 0;

    if (mCurrentPdeSolution)
    {
        VecGetSize(mCurrentPdeSolution, &size_of_soln_previous_step);
    }

    static_cast<AveragedSinksPde<DIM>*>(mpPde)->SetupSourceTerms(*mpCoarsePdeMesh);

    SimpleLinearEllipticAssembler<DIM,DIM> assembler(mpCoarsePdeMesh, mpPde, &bcc);

    if (size_of_soln_previous_step == (int)r_mesh.GetNumNodes())
    {
        // We make an initial guess which gets copied by the Solve method of
        // SimpleLinearSolver, so we need to delete it too.
        Vec initial_guess;
        VecDuplicate(mCurrentPdeSolution, &initial_guess);
        VecCopy(mCurrentPdeSolution, initial_guess);

        // Use current solution as the initial guess
        VecDestroy(mCurrentPdeSolution);    // Solve method makes its own mCurrentPdeSolution
        mCurrentPdeSolution = assembler.Solve(initial_guess);
        VecDestroy(initial_guess);
    }
    else
    {
        assert(mCurrentPdeSolution == NULL);
        /*
         * Eventually we will enable the coarse PDE mesh to change size, for example
         * in the case of a spheroid that grows a lot (see #630). In this case we should
         * uncomment the following code.
         *
        if (mCurrentPdeSolution)
        {
            assert(0);
            VecDestroy(mCurrentPdeSolution);
        }
        *
        */
        mCurrentPdeSolution = assembler.Solve();
    }

    /*
     * Update cellwise data - since the cells are not nodes on the coarse
     * mesh, we have to interpolate from the nodes of the coarse mesh onto
     * the cell locations.
     */
    ReplicatableVector solution_repl(mCurrentPdeSolution);

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        // Find coarse mesh element containing cell
        unsigned elem_index = FindCoarseElementContainingCell(*cell_iter);

        Element<DIM,DIM>* p_element = mpCoarsePdeMesh->GetElement(elem_index);

        const ChastePoint<DIM>& r_position_of_cell = this->mrTissue.GetLocationOfCellCentre(*cell_iter);

        c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(r_position_of_cell);

        double interpolated_solution = 0.0;
        for (unsigned i=0; i<DIM+1/*num_nodes*/; i++)
        {
            double nodal_value = solution_repl[ p_element->GetNodeGlobalIndex(i) ];
            interpolated_solution += nodal_value*weights(i);
        }

        unsigned index = this->mrTissue.GetLocationIndexUsingCell(*cell_iter);
        CellwiseData<DIM>::Instance()->SetValue(interpolated_solution, index);
    }
}

template<unsigned DIM>
unsigned TissueSimulationWithPdes<DIM>::FindCoarseElementContainingCell(TissueCell& rCell)
{
    // Get containing element at last timestep from mCellPdeElementMap
    unsigned old_element_index = mCellPdeElementMap[&rCell];

    // Create a std::set of guesses for the current containing element
    std::set<unsigned> test_elements;
    test_elements.insert(old_element_index);

    Element<DIM,DIM>* p_element = mpCoarsePdeMesh->GetElement(old_element_index);

    for (unsigned local_index=0; local_index<DIM+1; local_index++)
    {
        std::set<unsigned> element_indices = p_element->GetNode(local_index)->rGetContainingElementIndices();

        for (std::set<unsigned>::iterator iter = element_indices.begin();
             iter != element_indices.end();
             ++iter)
        {
            test_elements.insert(*iter);
        }
    }

    // Find new element, using the previous one as a guess
    const ChastePoint<DIM>& r_cell_position = this->mrTissue.GetLocationOfCellCentre(rCell);
    unsigned new_element_index = mpCoarsePdeMesh->GetContainingElementIndex(r_cell_position, false, test_elements);

    // Update mCellPdeElementMap
    mCellPdeElementMap[&rCell] = new_element_index;

    return new_element_index;
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::PostSolve()
{
    SolvePde();

    // Save results to file
    SimulationTime* p_time = SimulationTime::Instance();

    double time_next_step = p_time->GetTime() + p_time->GetTimeStep();

    if ((p_time->GetTimeStepsElapsed()+1)%this->mSamplingTimestepMultiple == 0)
    {
        WritePdeSolution(time_next_step);
    }

#define COVERAGE_IGNORE
    if (mWriteDailyAverageRadialPdeSolution)
    {
        ///\todo Worry about round-off errors
        unsigned num_timesteps_per_day = (unsigned) (DBL_EPSILON + 24/SimulationTime::Instance()->GetTimeStep());

        if ((p_time->GetTimeStepsElapsed()+1) % num_timesteps_per_day == 0)
        {
            WriteAverageRadialPdeSolution(time_next_step, mNumRadialIntervals);
        }
    }
#undef COVERAGE_IGNORE

}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::WritePdeSolution(double time)
{
    if (PetscTools::AmMaster())
    {
        // Since there are no ghost nodes, the number of nodes must equal the number of real cells
        assert(this->mrTissue.GetNumNodes() == this->mrTissue.GetNumRealCells());

        (*mpVizPdeSolutionResultsFile) << time << "\t";

        for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
             cell_iter != this->mrTissue.End();
             ++cell_iter)
        {
            unsigned global_index = this->mrTissue.GetLocationIndexUsingCell(*cell_iter);
            (*mpVizPdeSolutionResultsFile) << global_index << " ";

            const c_vector<double,DIM>& position = this->mrTissue.GetLocationOfCellCentre(*cell_iter);
            for (unsigned i=0; i<DIM; i++)
            {
                (*mpVizPdeSolutionResultsFile) << position[i] << " ";
            }

            double solution = CellwiseData<DIM>::Instance()->GetValue(*cell_iter);
            (*mpVizPdeSolutionResultsFile) << solution << " ";
        }
        (*mpVizPdeSolutionResultsFile) << "\n";
    }
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::SetWriteAverageRadialPdeSolution(unsigned numRadialIntervals, bool writeDailyResults)
{
    mWriteAverageRadialPdeSolution = true;
    mNumRadialIntervals = numRadialIntervals;
    mWriteDailyAverageRadialPdeSolution = writeDailyResults;
}

template<unsigned DIM>
void TissueSimulationWithPdes<DIM>::WriteAverageRadialPdeSolution(double time, unsigned numRadialIntervals)
{
    (*mpAverageRadialPdeSolutionResultsFile) << time << " ";

    // Calculate the centre of the tissue
    c_vector<double,DIM> centre = zero_vector<double>(DIM);
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
         cell_iter != this->mrTissue.End();
         ++cell_iter)
    {
       centre += (this->mrTissue.GetLocationOfCellCentre(*cell_iter));
    }
    centre /= ((double) this->mrTissue.GetNumNodes());

    // Calculate the distance between each node and the centre of the tissue, as well as the maximum of these
    std::map<double, TissueCell*> distance_cell_map;
    double max_distance_from_centre = 0.0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
         cell_iter != this->mrTissue.End();
         ++cell_iter)
    {
        double distance = norm_2(this->mrTissue.GetLocationOfCellCentre(*cell_iter) - centre);
        distance_cell_map[distance] = &(*cell_iter);

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

    // Calculate PDE solution in each radial interval
    double lower_radius = 0.0;
    for (unsigned i=0; i<numRadialIntervals; i++)
    {
        unsigned counter = 0;
        double average_solution = 0.0;

        for (std::map<double, TissueCell*>::iterator iter = distance_cell_map.begin();
             iter != distance_cell_map.end();
             ++iter)
        {
            if (iter->first > lower_radius && iter->first <= radius_intervals[i])
            {
                average_solution += CellwiseData<DIM>::Instance()->GetValue(*(iter->second));
                counter++;
            }
        }
        if (counter > 0)
        {
            average_solution /= (double) counter;
        }

        // Write results to file
        (*mpAverageRadialPdeSolutionResultsFile) << radius_intervals[i] << " " << average_solution << " ";
        lower_radius = radius_intervals[i];
    }
    (*mpAverageRadialPdeSolutionResultsFile) << "\n";
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class TissueSimulationWithPdes<1>;
template class TissueSimulationWithPdes<2>;
template class TissueSimulationWithPdes<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TissueSimulationWithPdes)
