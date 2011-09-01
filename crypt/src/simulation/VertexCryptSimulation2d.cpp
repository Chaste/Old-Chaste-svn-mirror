/*

Copyright (C) University of Oxford, 2005-2011

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

#include "VertexCryptSimulation2d.hpp"
#include "WntConcentration.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"

VertexCryptSimulation2d::VertexCryptSimulation2d(AbstractCellPopulation<2>& rCellPopulation,
                                                 bool deleteCellPopulationAndForceCollection,
                                                 bool initialiseCells)
    : OffLatticeSimulation<2>(rCellPopulation,
                             deleteCellPopulationAndForceCollection,
                             initialiseCells),
      mWriteBetaCatenin(false)
{
    /*
     * To check if beta-catenin results will be written to file, we test if the first
     * cell has a cell-cycle model that is a subclass of AbstractVanLeeuwen2009WntSwatCellCycleModel.
     * In doing so, we assume that all cells in the simulation have the same cell cycle
     * model.
     */
    if (dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(mrCellPopulation.Begin()->GetCellCycleModel()))
    {
        mWriteBetaCatenin = true;
    }

    if (!mDeleteCellPopulationAndForcesAndBCsInDestructor)
    {
		// Pass a CryptSimulationBoundaryCondition object into mBoundaryConditions
		CryptSimulationBoundaryCondition<2>* p_boundary_condition = new CryptSimulationBoundaryCondition<2>(&rCellPopulation);
		AddCellPopulationBoundaryCondition(p_boundary_condition);
    }
}

VertexCryptSimulation2d::~VertexCryptSimulation2d()
{
	// Delete the CryptSimulationBoundaryCondition object from mBoundaryConditions
	for (std::vector<AbstractCellPopulationBoundaryCondition<2>*>::iterator it=mBoundaryConditions.begin();
         it != mBoundaryConditions.end();
         ++it)
    {
		delete *it;
    }
    // Now clear the container in case another piece of code tries to delete it
    mBoundaryConditions.clear();
}

c_vector<double, 2> VertexCryptSimulation2d::CalculateCellDivisionVector(CellPtr pParentCell)
{
    c_vector<double, 2> axis_of_division = zero_vector<double>(2);

    // We don't need to prescribe how 'stem' cells divide if Wnt is present
    bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<2>::Destroy();
        if (pParentCell->GetCellCycleModel()->GetCellProliferativeType() == STEM)
        {
            axis_of_division(0) = 1.0;
            axis_of_division(1) = 0.0;
        }
    }
    return axis_of_division;
}

void VertexCryptSimulation2d::WriteVisualizerSetupFile()
{
    *mpVizSetupFile << "MeshWidth\t" << mrCellPopulation.GetWidth(0) << "\n";
}

void VertexCryptSimulation2d::SetupWriteBetaCatenin()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);
    mVizBetaCateninResultsFile = output_file_handler.OpenOutputFile("results.vizbetacatenin");
    *mpVizSetupFile << "BetaCatenin\n";
}

void VertexCryptSimulation2d::WriteBetaCatenin(double time)
{
    *mVizBetaCateninResultsFile <<  time << "\t";

    unsigned global_index;
    double x;
    double y;
    double b_cat_membrane;
    double b_cat_cytoplasm;
    double b_cat_nuclear;

    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        global_index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        x = mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[0];
        y = mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

        // We should only be calling this code block if mWriteBetaCatenin has been set to true in the constructor
        assert(mWriteBetaCatenin);

        AbstractVanLeeuwen2009WntSwatCellCycleModel* p_model = dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(cell_iter->GetCellCycleModel());
        b_cat_membrane = p_model->GetMembraneBoundBetaCateninLevel();
        b_cat_cytoplasm = p_model->GetCytoplasmicBetaCateninLevel();
        b_cat_nuclear = p_model->GetNuclearBetaCateninLevel();

        *mVizBetaCateninResultsFile << global_index << " " << x << " " << y << " " << b_cat_membrane << " " << b_cat_cytoplasm << " " << b_cat_nuclear << " ";
    }

    *mVizBetaCateninResultsFile << "\n";
}

void VertexCryptSimulation2d::SetupSolve()
{
    /*
     * If there are any cells in the simulation, and mWriteBetaCatenin has been set to true in the constructor,
     * then set up the beta-catenin results file and write the initial conditions to file.
     */
    bool any_cells_present = (mrCellPopulation.Begin() != mrCellPopulation.End());
    if (any_cells_present && mWriteBetaCatenin)
    {
        SetupWriteBetaCatenin();
        double current_time = SimulationTime::Instance()->GetTime();
        WriteBetaCatenin(current_time);
    }
}

void VertexCryptSimulation2d::PostSolve()
{
    SimulationTime* p_time = SimulationTime::Instance();

    if ((p_time->GetTimeStepsElapsed()+1)%mSamplingTimestepMultiple == 0)
    {
        /*
         * If there are any cells in the simulation, and mWriteBetaCatenin has been set
         * to true in the constructor, then set up the beta-catenin results file and
         * write the initial conditions to file.
         */
        bool any_cells_present = (mrCellPopulation.Begin() != mrCellPopulation.End());
        if (any_cells_present && mWriteBetaCatenin)
        {
            double time_next_step = p_time->GetTime() + p_time->GetTimeStep();
            WriteBetaCatenin(time_next_step);
        }
    }
}

void VertexCryptSimulation2d::AfterSolve()
{
    /*
     * If there are any cells in the simulation, and mWriteBetaCatenin has been set
     * to true in the constructor, then close the beta-catenin results file.
     */
    bool any_cells_present = (mrCellPopulation.Begin() != mrCellPopulation.End());
    if (any_cells_present && mWriteBetaCatenin)
    {
        mVizBetaCateninResultsFile->close();
    }
}

void VertexCryptSimulation2d::UseJiggledBottomCells()
{
    // The CryptSimulationBoundaryCondition object is the first element of mBoundaryConditions
    static_cast<CryptSimulationBoundaryCondition<2>*>(mBoundaryConditions[0])->SetUseJiggledBottomCells(true);
}

void VertexCryptSimulation2d::SetBottomCellAncestors()
{
    unsigned index = 0;
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        if (mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[1] < 1.0)
        {
            cell_iter->SetAncestor(index++);
        }
    }
}

void VertexCryptSimulation2d::OutputSimulationParameters(out_stream& rParamsFile)
{
    double width = mrCellPopulation.GetWidth(0);
    bool use_jiggled_bottom_cells = static_cast<CryptSimulationBoundaryCondition<2>*>(mBoundaryConditions[0])->GetUseJiggledBottomCells();

    *rParamsFile << "\t\t<CryptCircumference>" << width << "</CryptCircumference>\n";
    *rParamsFile << "\t\t<UseJiggledBottomCells>" << use_jiggled_bottom_cells << "</UseJiggledBottomCells>\n";

    // Call method on direct parent class
    OffLatticeSimulation<2>::OutputSimulationParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VertexCryptSimulation2d)
