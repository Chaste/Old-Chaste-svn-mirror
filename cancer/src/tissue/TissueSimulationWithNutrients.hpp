/*

Copyright (C) University of Oxford, 2008

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
#ifndef TISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TISSUESIMULATIONWITHNUTRIENTS_HPP_

#include "TissueSimulation.hpp"
#include "SimpleDataWriter.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "TissueSimulationWithNutrientsAssembler.hpp"
#include "CellwiseData.hpp"
#include "PetscTools.hpp"
#include "AveragedSinksPde.hpp"

/// \todo This class needs documenting (see #736)
template<unsigned DIM>
class TissueSimulationWithNutrients : public TissueSimulation<DIM>
{
    // Allow tests to access private members, in order to test computation of private functions
    friend class TestTissueSimulationWithNutrients;

private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<TissueSimulation<DIM> >(*this);
        archive & mWriteAverageRadialNutrientResults;
        archive & mWriteDailyAverageRadialNutrientResults;
        archive & mNumRadialIntervals;
        archive & mCellNutrientElementMap;
    }

    /**
     *  Current nutrient concentration, for use as an initial guess
     *  when solving the nutrient PDE.
     */
    Vec mNutrientSolution;

    /**
     *  Pointer to the PDE satisfied by the nutrient.
     */
    AbstractLinearEllipticPde<DIM>* mpPde;

    /**
     *  Pointer to the averaged sink PDE satisfied by the nutrient.
     */
    AveragedSinksPde<DIM>* mpAveragedSinksPde;

    /**
     *  File that the nutrient values are written out to.
     */
    out_stream mpNutrientResultsFile;

    /**
     *  File that the average radial nutrient distribution is written out to.
     */
    out_stream mpAverageRadialNutrientResultsFile;

    /**
     *  Whether to write to file the average radial nutrient distribution.
     */
    bool mWriteAverageRadialNutrientResults;

    /**
     *  Whether to write the average radial nutrient distribution DAILY.
     */
    bool mWriteDailyAverageRadialNutrientResults;

    /**
     *  Number of radial 'bins' used to calculate the average
     *  radial nutrient distribution.
     */
    unsigned mNumRadialIntervals;

    /**
     *  Coarse nutrient mesh on which to solve the nutrient PDE.
     */
    TetrahedralMesh<DIM,DIM>* mpCoarseNutrientMesh;

    /**
     * Map between cells and the elements of the coarse nutrient mesh containing them.
     */
    std::map<TissueCell*, unsigned> mCellNutrientElementMap;

    /**
     *  Overridden SetupSolve() method.
     */
    void SetupSolve();

    /**
     *  Set up the nutrient writer.
     */
    void SetupWriteNutrient();

    /**
     *  Write the nutrient distribution to file at a specified time.
     *
     * @param time The time at which to record the nutrient distribution
     */
    void WriteNutrient(double time);

    /**
     *  Write the average radial nutrient distribution to file at a specified time.
     *
     * @param time The time at which to record the average radial nutrient distribution
     * @param numIntervals  The number of radial intervals in which the average nutrient concentration is calculated
     */
    void WriteAverageRadialNutrientDistribution(double time, unsigned numIntervals);

    /**
     *  Solve the nutrient PDE.
     */
    void SolveNutrientPde();

    /**
     *  Solve the nutrient PDE on a coarse mesh.
     */
    void SolveNutrientPdeUsingCoarseMesh();

    /**
     *  Find the index of the coarse mesh element containing rCell.
     */
    unsigned FindElementContainingCell(TissueCell& rCell);

    /**
     *  Overridden PostSolve() method.
     */
    void PostSolve();

    /**
     *  Overridden AfterSolve() method.
     */
    void AfterSolve();

    /**
     *  Create a coarse mesh on which to solve the nutrient PDE.
     */
    void CreateCoarseNutrientMesh(double coarseGrainScaleFactor);

    /**
     *  Initialise the std::map mCellNutrientElementMap.
     */
    void InitialiseCoarseNutrientMesh();

public:

    /**
     * Constructor
     *
     * @param rTissue A tissue facade class (contains a mesh and cells)
     * @param forceCollection The mechanics to use in the simulation
     * @param pPde The PDE for the nutrient concentration(s)
     * @param pAveragedSinksPde The PDE for the nutrient concentration(s)
     * @param deleteTissue whether to delete the tissue on destruction to free up memory
     * @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     *
     */
     TissueSimulationWithNutrients(AbstractTissue<DIM>& rTissue,
                                   std::vector<AbstractForce<DIM>*> forceCollection,
                                   AbstractLinearEllipticPde<DIM>* pPde=NULL,
                                   AveragedSinksPde<DIM>* pAveragedSinksPde=NULL,
                                   bool deleteTissueAndForceCollection=false,
                                   bool initialiseCells=true);
                     
    /**
     * Destructor
     *
     * Free any memory allocated by the constructor.
     * This frees the current nutrient distribution, if it exists.
     */
    ~TissueSimulationWithNutrients();

    /**
     * A small hack until we fully archive this class -
     * needed to set the PDE after loading a simulation
     * from an archive.
     */
    void SetPde(AbstractLinearEllipticPde<DIM>* pPde);

    /**
     * A small hack until we fully archive this class -
     * needed to set the PDE after loading a simulation
     * from an archive.
     */
    void SetAveragedSinksPde(AveragedSinksPde<DIM>* pAveragedSinksPde);

    /**
     *  Get the current nutrient solution
     */
    Vec GetNutrientSolution();

    /**
     * Write the final (and optionally also the daily) average
     * radial nutrient distribution to file.
     *
     * @param numRadialIntervals The number of radial intervals in which the average nutrient concentration is calculated
     * @param writeDailyResults Whether to record the average radial nutrient distribution at the end of each day of the simulation
     */

    void SetWriteAverageRadialNutrientResults(unsigned numRadialIntervals=10,
                                              bool writeDailyResults=false);

    /**
     * Solve the nutrient PDE on a coarse mesh.
     *
     * @param coarseGrainScaleFactor The ratio of the width of the coarse nutrient mesh to the initial width of the tissue
     */
    void UseCoarseNutrientMesh(double coarseGrainScaleFactor=10.0);

};
                   
template<unsigned DIM>
TissueSimulationWithNutrients<DIM>::TissueSimulationWithNutrients(AbstractTissue<DIM>& rTissue,
                                   std::vector<AbstractForce<DIM>*> forceCollection,
                                   AbstractLinearEllipticPde<DIM>* pPde,
                                   AveragedSinksPde<DIM>* pAveragedSinksPde,
                                   bool deleteTissueAndForceCollection,
                                   bool initialiseCells)
    : TissueSimulation<DIM>(rTissue, 
                            forceCollection, 
                            deleteTissueAndForceCollection, 
                            initialiseCells),
      mNutrientSolution(NULL),
      mpPde(pPde),
      mpAveragedSinksPde(pAveragedSinksPde),
      mWriteAverageRadialNutrientResults(false),
      mWriteDailyAverageRadialNutrientResults(false),
      mNumRadialIntervals(0), // 'unset' value
      mpCoarseNutrientMesh(NULL)
{
}

template<unsigned DIM>
TissueSimulationWithNutrients<DIM>::~TissueSimulationWithNutrients()
{
    if (mNutrientSolution)
    {
        VecDestroy(mNutrientSolution);
    }
    if (mpCoarseNutrientMesh)
    {
        delete mpCoarseNutrientMesh;
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetPde(AbstractLinearEllipticPde<DIM>* pPde)
{
    mpPde = pPde;
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetAveragedSinksPde(AveragedSinksPde<DIM>* pAveragedSinksPde)
{
    mpAveragedSinksPde = pAveragedSinksPde;
}

template<unsigned DIM>
Vec TissueSimulationWithNutrients<DIM>::GetNutrientSolution()
{
    return mNutrientSolution;
}

//////////////////////////////////////////////////////////////////////////////
//                          Setup/AfterSolve methods                        //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetupSolve()
{
    if (mpCoarseNutrientMesh!=NULL)
    {
        InitialiseCoarseNutrientMesh();
    }
    if (this->mrTissue.Begin() != this->mrTissue.End())
    {
        SetupWriteNutrient();
        double current_time = SimulationTime::Instance()->GetTime();
        WriteNutrient(current_time);
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetupWriteNutrient()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/",false);
    if (output_file_handler.IsMaster())
    {
        mpNutrientResultsFile = output_file_handler.OpenOutputFile("results.viznutrient");
        *this->mpSetupFile << "Nutrient \n";
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
    // Create coarse nutrient mesh (can use a larger mesh if required, e.g. disk_984_elements)
    TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
    mpCoarseNutrientMesh = new TetrahedralMesh<2,2>;
    mpCoarseNutrientMesh->ConstructFromMeshReader(mesh_reader);

    // Find centre of tissue
    c_vector<double,2> centre_of_tissue = zero_vector<double>(2);
    for (typename MeshBasedTissue<2>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        centre_of_tissue += cell_iter.rGetLocation();
    }
    centre_of_tissue /= this->mrTissue.GetNumRealCells();

    // Find max radius of tissue
    double max_tissue_radius = 0.0;
    for (typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        double radius = norm_2(centre_of_tissue - cell_iter.rGetLocation() );
        if (radius > max_tissue_radius)
        {
            max_tissue_radius = radius;
        }
    }

    // Find centre of coarse nutrient mesh
    c_vector<double,2> centre_of_nutrient_mesh = zero_vector<double>(2);

    for (unsigned i=0; i<mpCoarseNutrientMesh->GetNumNodes(); i++)
    {
        centre_of_nutrient_mesh += mpCoarseNutrientMesh->GetNode(i)->rGetLocation();
    }
    centre_of_nutrient_mesh /= mpCoarseNutrientMesh->GetNumNodes();

    // Find max radius of coarse nutrient mesh
    double max_mesh_radius = 0.0;
    for (unsigned i=0; i<mpCoarseNutrientMesh->GetNumNodes(); i++)
    {
        double radius = norm_2(centre_of_nutrient_mesh - mpCoarseNutrientMesh->GetNode(i)->rGetLocation());
        if (radius > max_mesh_radius)
        {
            max_mesh_radius = radius;
        }
    }

    // Translate centre of coarse nutrient mesh to the origin
    mpCoarseNutrientMesh->Translate(-centre_of_nutrient_mesh[0], -centre_of_nutrient_mesh[1]);

    // Scale nutrient mesh
    double scale_factor = (max_tissue_radius/max_mesh_radius)*coarseGrainScaleFactor;
    mpCoarseNutrientMesh->Scale(scale_factor, scale_factor);

    // Translate centre of coarse nutrient mesh to centre of the tissue
    mpCoarseNutrientMesh->Translate(centre_of_tissue[0], centre_of_tissue[1]);
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::InitialiseCoarseNutrientMesh()
{
    mCellNutrientElementMap.clear();

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        // Find the element of mpCoarseNutrientMesh that contains this cell
        const ChastePoint<DIM>& r_position_of_cell = cell_iter.rGetLocation();
        unsigned elem_index = mpCoarseNutrientMesh->GetContainingElementIndex(r_position_of_cell);
        mCellNutrientElementMap[&(*cell_iter)] = elem_index;
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::AfterSolve()
{
    TissueSimulation<DIM>::AfterSolve();

    if (this->mrTissue.Begin() != this->mrTissue.End() // if there are any cells
    && PetscTools::AmMaster())
    {
        mpNutrientResultsFile->close();

        if (mWriteAverageRadialNutrientResults)
        {
            WriteAverageRadialNutrientDistribution(SimulationTime::Instance()->GetTime(), mNumRadialIntervals);
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
    if (mpCoarseNutrientMesh!=NULL)
    {
        SolveNutrientPdeUsingCoarseMesh();
        return;
    }

    assert(mpAveragedSinksPde == NULL);
    assert(mpPde);

    // Note: If not using a coarse nutrient mesh, we MUST be using a MeshBasedTissue

    TetrahedralMesh<DIM,DIM>& r_mesh = static_cast<MeshBasedTissue<DIM>*>(&(this->mrTissue))->rGetMesh();
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    // We shouldn't have any ghost nodes in a TissueSimulationWithNutrients
    assert(this->mrTissue.HasGhostNodes()==false);

    // Set up boundary conditions
    BoundaryConditionsContainer<DIM,DIM,1> bcc;
    ConstBoundaryCondition<DIM>* p_boundary_condition = new ConstBoundaryCondition<DIM>(1.0);
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = r_mesh.GetBoundaryNodeIteratorBegin();
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

    if (mNutrientSolution)
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
            assert(size_of_soln_previous_step != 0);
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

    TetrahedralMesh<DIM,DIM>& r_mesh = *mpCoarseNutrientMesh;
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    // We shouldn't have any ghost nodes in a TissueSimulationWithNutrients
    assert(this->mrTissue.HasGhostNodes()==false);

    // Loop over cells and calculate centre of distribution
    c_vector<double, DIM> centre = zero_vector<double>(DIM);
    for (typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        centre += cell_iter.rGetLocation();
    }
    centre /= this->mrTissue.GetNumRealCells();

    // Find max radius
    double max_radius = 0.0;
    for (typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        double radius = norm_2(centre - cell_iter.rGetLocation());
        if (radius > max_radius)
        {
            max_radius = radius;
        }
    }

    // Set up boundary conditions
    BoundaryConditionsContainer<DIM,DIM,1> bcc;
    ConstBoundaryCondition<DIM>* p_boundary_condition = new ConstBoundaryCondition<DIM>(1.0);

    // Get the set of coarse element indices that contain tissue cells
    std::set<unsigned> coarse_element_indices_in_map;
    for (typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        coarse_element_indices_in_map.insert(mCellNutrientElementMap[&(*cell_iter)]);
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
                unsigned node_index = p_element->GetNode(local_index)->GetIndex();
                coarse_mesh_boundary_node_indices.insert(node_index);
            }
        }
    }

    // Apply boundary condition to the nodes in the set coarse_mesh_boundary_node_indices
    for (std::set<unsigned>::iterator iter = coarse_mesh_boundary_node_indices.begin();
         iter != coarse_mesh_boundary_node_indices.end();
         ++iter)
    {
        bcc.AddDirichletBoundaryCondition(r_mesh.GetNode(*iter), p_boundary_condition, 0, false);
    }

    PetscInt size_of_soln_previous_step = 0;

    if (mNutrientSolution)
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
        assert(mNutrientSolution == NULL);
        /**
         * Eventually we will enable the coarse nutrient mesh to change size, for example
         * in the case of a spheroid that grows a lot (see #630). In this case we should
         * uncomment the following code.
         *
        if (mNutrientSolution)
        {
            assert(0);
            VecDestroy(mNutrientSolution);
        }
        *
        */
        mNutrientSolution = assembler.Solve();
    }

    // Update cellwise data - since the cells are not nodes on the coarse
    // mesh, we have to interpolate from the nodes of the coarse mesh onto
    // the cell locations
    ReplicatableVector nutrient_repl(mNutrientSolution);

    for (typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        // Find coarse mesh element containing cell
        unsigned elem_index = FindElementContainingCell(*cell_iter);

        Element<DIM,DIM>* p_element = mpCoarseNutrientMesh->GetElement(elem_index);

        const ChastePoint<DIM>& r_position_of_cell = cell_iter.rGetLocation();

        c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(r_position_of_cell);

        double interpolated_nutrient = 0.0;
        for (unsigned i=0; i<DIM+1/*num_nodes*/; i++)
        {
            double nodal_value = nutrient_repl[ p_element->GetNodeGlobalIndex(i) ];
            interpolated_nutrient += nodal_value*weights(i);
        }

        CellwiseData<DIM>::Instance()->SetValue(interpolated_nutrient, cell_iter.GetNode());
    }
}

template<unsigned DIM>
unsigned TissueSimulationWithNutrients<DIM>::FindElementContainingCell(TissueCell& rCell)
{
    // Get containing element at last timestep from mCellNutrientElementMap
    unsigned old_element_index = mCellNutrientElementMap[&rCell];

    // Create a std::set of guesses for the current containing element
    std::set<unsigned> test_elements;
    test_elements.insert(old_element_index);

    Element<DIM,DIM>* p_element = mpCoarseNutrientMesh->GetElement(old_element_index);

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
    const ChastePoint<DIM>& r_cell_position = this->mrTissue.GetLocationOfCell(rCell);
    unsigned new_element_index = mpCoarseNutrientMesh->GetContainingElementIndex(r_cell_position, false, test_elements);

    // Update mCellNutrientElementMap
    mCellNutrientElementMap[&rCell] = new_element_index;

    return new_element_index;
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::PostSolve()
{
    SolveNutrientPde();

    // Save results to file
    SimulationTime* p_time = SimulationTime::Instance();

    double time_next_step = p_time->GetTime() + p_time->GetTimeStep();

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
        assert(this->mrTissue.GetNumNodes()==this->mrTissue.GetNumRealCells());

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
    TetrahedralMesh<DIM,DIM>& r_mesh = static_cast<MeshBasedTissue<DIM>*>(&(this->mrTissue))->rGetMesh();
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
        distance_cell_map[distance] = &(this->mrTissue.rGetCellUsingLocationIndex(i));

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




#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TissueSimulationWithNutrients)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulationWithNutrients.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TissueSimulationWithNutrients<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractForce<DIM>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TissueSimulationWithNutrients<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<DIM>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialize instance
    ::new(t)TissueSimulationWithNutrients<DIM>(*p_tissue, force_collection, NULL, NULL, true, false);
}
}
} // namespace ...


#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
