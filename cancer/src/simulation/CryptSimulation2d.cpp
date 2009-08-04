/*

Copyright (C) University of Oxford, 2005-2009

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

#include "CryptSimulation2d.hpp"
#include "WntConcentration.hpp"
#include "IngeWntSwatCellCycleModel.hpp"


CryptSimulation2d::CryptSimulation2d(AbstractTissue<2>& rTissue,
                  std::vector<AbstractForce<2>*> forceCollection,
                  bool deleteTissueAndForceCollection,
                  bool initialiseCells)
    : TissueSimulation<2>(rTissue,
                          forceCollection,
                          deleteTissueAndForceCollection,
                          initialiseCells),
      mUseJiggledBottomCells(false)
{
    mpStaticCastTissue = static_cast<MeshBasedTissueWithGhostNodes<2>*>(&mrTissue);
}


c_vector<double, 2> CryptSimulation2d::CalculateCellDivisionVector(TissueCell* pParentCell)
{
    // Location of parent and daughter cells
    c_vector<double, 2> parent_coords = mpStaticCastTissue->GetLocationOfCellCentre(pParentCell);
    c_vector<double, 2> daughter_coords;

    // Get separation parameter
    double separation = TissueConfig::Instance()->GetDivisionSeparation();

    // Make a random direction vector of the required length
    c_vector<double, 2> random_vector;

    /*
     * Pick a random direction and move the parent cell backwards by 0.5*separation
     * in that direction and return the position of the daughter cell 0.5*separation
     * forwards in that direction.
     */

    double random_angle = RandomNumberGenerator::Instance()->ranf();
    random_angle *= 2.0*M_PI;

    random_vector(0) = 0.5*separation*cos(random_angle);
    random_vector(1) = 0.5*separation*sin(random_angle);

    c_vector<double, 2> proposed_new_parent_coords = parent_coords - random_vector;
    c_vector<double, 2> proposed_new_daughter_coords = parent_coords + random_vector;

    if (   (proposed_new_parent_coords(1) >= 0.0)
        && (proposed_new_daughter_coords(1) >= 0.0))
    {
        // We are not too close to the bottom of the tissue, so move parent
        parent_coords = proposed_new_parent_coords;
        daughter_coords = proposed_new_daughter_coords;
    }
    else
    {
        proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
        while (proposed_new_daughter_coords(1) < 0.0)
        {
            random_angle = RandomNumberGenerator::Instance()->ranf();
            random_angle *= 2.0*M_PI;

            random_vector(0) = separation*cos(random_angle);
            random_vector(1) = separation*sin(random_angle);
            proposed_new_daughter_coords = parent_coords + random_vector;
        }
        daughter_coords = proposed_new_daughter_coords;
    }

    assert(daughter_coords(1) >= 0.0); // to make sure dividing cells stay in the tissue
    assert(parent_coords(1) >= 0.0);   // to make sure dividing cells stay in the tissue

    // Set the parent to use this location
    ChastePoint<2> parent_coords_point(parent_coords);

    unsigned node_index = mpStaticCastTissue->GetLocationIndexUsingCell(pParentCell);
    mrTissue.SetNode(node_index, parent_coords_point);

    return daughter_coords;
}


void CryptSimulation2d::WriteVisualizerSetupFile()
{
    *mpSetupFile << "MeshWidth\t" << mpStaticCastTissue->rGetMesh().GetWidth(0u) << "\n";// get furthest distance between nodes in the x-direction
}


void CryptSimulation2d::SetupWriteBetaCatenin()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);
    mBetaCatResultsFile = output_file_handler.OpenOutputFile("results.vizbetacatenin");
    *mpSetupFile << "BetaCatenin\n";
}


void CryptSimulation2d::WriteBetaCatenin(double time)
{
    *mBetaCatResultsFile <<  time << "\t";

    unsigned global_index;
    double x;
    double y;
    double b_cat_membrane;
    double b_cat_cytoplasm;
    double b_cat_nuclear;

    for (AbstractTissue<2>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        global_index = mpStaticCastTissue->GetLocationIndexUsingCell(&(*cell_iter));
        x = mpStaticCastTissue->GetLocationOfCellCentre(&(*cell_iter))[0];
        y = mpStaticCastTissue->GetLocationOfCellCentre(&(*cell_iter))[1];

        // If writing beta-catenin, the model has to be an IngeWntSwatCellCycleModel
        IngeWntSwatCellCycleModel *p_model = static_cast<IngeWntSwatCellCycleModel*>(cell_iter->GetCellCycleModel());

        b_cat_membrane = p_model->GetMembraneBoundBetaCateninLevel();
        b_cat_cytoplasm = p_model->GetCytoplasmicBetaCateninLevel();
        b_cat_nuclear = p_model->GetNuclearBetaCateninLevel();

        *mBetaCatResultsFile << global_index << " " << x << " " << y << " " << b_cat_membrane << " " << b_cat_cytoplasm << " " << b_cat_nuclear << " ";
    }

    *mBetaCatResultsFile << "\n";
}


void CryptSimulation2d::SetupSolve()
{
    if (   (mrTissue.Begin() != mrTissue.End()) // there are any cells
        && (dynamic_cast<IngeWntSwatCellCycleModel*>(mrTissue.Begin()->GetCellCycleModel())) ) // assume all the cells are the same
    {
        SetupWriteBetaCatenin();
        double current_time = SimulationTime::Instance()->GetTime();
        WriteBetaCatenin(current_time);
    }
}


void CryptSimulation2d::PostSolve()
{
    SimulationTime *p_time = SimulationTime::Instance();

    if ((p_time->GetTimeStepsElapsed()+1)%mSamplingTimestepMultiple==0)
    {
        if (   (mrTissue.Begin() != mrTissue.End()) // there are any cells
            && (dynamic_cast<IngeWntSwatCellCycleModel*>(mrTissue.Begin()->GetCellCycleModel())) ) // assume all the cells are the same
        {
            double time_next_step = p_time->GetTime() + p_time->GetTimeStep();
            WriteBetaCatenin(time_next_step);
        }
    }
}


void CryptSimulation2d::AfterSolve()
{
    if (   (mrTissue.Begin() != mrTissue.End()) // there are any cells
        && (dynamic_cast<IngeWntSwatCellCycleModel*>(mrTissue.Begin()->GetCellCycleModel())) ) // assume all the cells are the same
    {
        mBetaCatResultsFile->close();
    }
}


void CryptSimulation2d::UseJiggledBottomCells()
{
    mUseJiggledBottomCells = true;
}


void CryptSimulation2d::ApplyTissueBoundaryConditions(const std::vector< c_vector<double, 2> >& rOldLocations)
{
    bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<2>::Destroy();
    }

    // Iterate over all nodes associated with real cells to update their positions
    // according to any tissue boundary conditions
    for (AbstractTissue<2>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = mpStaticCastTissue->GetLocationIndexUsingCell(&(*cell_iter));

        // Get pointer to this node
        Node<2> *p_node = mpStaticCastTissue->GetNodeCorrespondingToCell(&(*cell_iter));

        if (!is_wnt_included)
        {
            /**
             * If WntConcentration is not set up then stem cells must be pinned,
             * so we reset the location of each stem cell.
             */
            if (cell_iter->GetCellType()==STEM)
            {
                // Get old node location
                c_vector<double, 2> old_node_location = rOldLocations[node_index];

                // Return node to old location
                p_node->rGetModifiableLocation()[0] = old_node_location[0];
                p_node->rGetModifiableLocation()[1] = old_node_location[1];
            }
        }

        // Any cell that has moved below the bottom of the crypt must be moved back up
        if (p_node->rGetLocation()[1] < 0.0)
        {
            p_node->rGetModifiableLocation()[1] = 0.0;
            if (mUseJiggledBottomCells)
            {
               /*
                * Here we give the cell a push upwards so that it doesn't
                * get stuck on the bottom of the crypt (as per #422).
                *
                * Note that all stem cells may get moved to the same height, so
                * we use a random perturbation to help ensure we are not simply
                * faced with the same problem at a different height!
                */
                p_node->rGetModifiableLocation()[1] = 0.05*mpRandomGenerator->ranf();
            }
        }
        assert(p_node->rGetLocation()[1] >= 0.0);
    }
}
