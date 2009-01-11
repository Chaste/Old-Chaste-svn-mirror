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


#include "CryptSimulation2d.hpp"
#include "WntConcentration.hpp"

c_vector<double, 2> CryptSimulation2d::CalculateDividingCellCentreLocations(AbstractTissue<2>::Iterator parentCell)
{
    double separation = CancerParameters::Instance()->GetDivisionSeparation();
    c_vector<double, 2> parent_coords = parentCell.rGetLocation();
    c_vector<double, 2> daughter_coords;

    // Pick a random direction and move the parent cell backwards by 0.5*sep in that
    // direction and return the position of the daughter cell (0.5*sep forwards in the
    // random vector direction

    // Make a random direction vector of the required length
    c_vector<double, 2> random_vector;

    double random_angle = RandomNumberGenerator::Instance()->ranf();
    random_angle *= 2.0*M_PI;

    random_vector(0) = 0.5*separation*cos(random_angle);
    random_vector(1) = 0.5*separation*sin(random_angle);

    c_vector<double, 2> proposed_new_parent_coords = parent_coords-random_vector;
    c_vector<double, 2> proposed_new_daughter_coords = parent_coords+random_vector;

    if (   (proposed_new_parent_coords(1) >= 0.0)
        && (proposed_new_daughter_coords(1) >= 0.0))
    {
        // We are not too close to the bottom of the tissue
        // move parent
        parent_coords = proposed_new_parent_coords;
        daughter_coords = proposed_new_daughter_coords;
    }
    else
    {
        proposed_new_daughter_coords = parent_coords+2.0*random_vector;
        while (proposed_new_daughter_coords(1) < 0.0)
        {
            random_angle = RandomNumberGenerator::Instance()->ranf();
            random_angle *= 2.0*M_PI;

            random_vector(0) = separation*cos(random_angle);
            random_vector(1) = separation*sin(random_angle);
            proposed_new_daughter_coords = parent_coords+random_vector;
        }
        daughter_coords = proposed_new_daughter_coords;
    }

    assert(daughter_coords(1)>=0.0); // to make sure dividing cells stay in the tissue
    assert(parent_coords(1)>=0.0);   // to make sure dividing cells stay in the tissue

    // Set the parent to use this location
    ChastePoint<2> parent_coords_point(parent_coords);
    mrTissue.MoveCell(parentCell, parent_coords_point);
    return daughter_coords;
}


void CryptSimulation2d::WriteVisualizerSetupFile()
{
    *mpSetupFile << "MeshWidth\t" << mpStaticCastTissue->rGetMesh().GetWidth(0u) << "\n";// get furthest distance between nodes in the x-direction
}


void CryptSimulation2d::SetupWriteBetaCatenin()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/",false);
    mBetaCatResultsFile = output_file_handler.OpenOutputFile("results.vizbetacatenin");
    *mpSetupFile << "BetaCatenin\n";
}


void CryptSimulation2d::WriteBetaCatenin(double time)
{
    *mBetaCatResultsFile <<  time << "\t";

    double global_index;
    double x;
    double y;
    double b_cat_membrane;
    double b_cat_cytoplasm;
    double b_cat_nuclear;

    for (AbstractTissue<2>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        global_index = (double) cell_iter.GetNode()->GetIndex();
        x = cell_iter.rGetLocation()[0];
        y = cell_iter.rGetLocation()[1];

        // If writing beta-catenin, the model has be be IngeWntSwatCellCycleModel
        IngeWntSwatCellCycleModel* p_model = dynamic_cast<IngeWntSwatCellCycleModel*>(cell_iter->GetCellCycleModel());

        b_cat_membrane = p_model->GetMembraneBoundBetaCateninLevel();
        b_cat_cytoplasm = p_model->GetCytoplasmicBetaCateninLevel();
        b_cat_nuclear = p_model->GetNuclearBetaCateninLevel();

        *mBetaCatResultsFile << global_index << " " << x << " " << y << " " << b_cat_membrane << " " << b_cat_cytoplasm << " " << b_cat_nuclear << " ";
    }

    *mBetaCatResultsFile << "\n";
}


void CryptSimulation2d::SetupSolve()
{
    if (   ( mrTissue.Begin() != mrTissue.End() )  // there are any cells
        && ( mrTissue.Begin()->GetCellCycleModel()->UsesBetaCat()) ) // assume all the cells are the same
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
        if (   ( mrTissue.Begin() != mrTissue.End() )  // there are any cells
            && ( mrTissue.Begin()->GetCellCycleModel()->UsesBetaCat()) ) // assume all the cells are the same
        {
            double time_next_step = p_time->GetTime() + p_time->GetTimeStep();
            WriteBetaCatenin(time_next_step);
        }
    }
}


void CryptSimulation2d::AfterSolve()
{
    if (   ( mrTissue.Begin() != mrTissue.End() )  // there are any cells
        && ( mrTissue.Begin()->GetCellCycleModel()->UsesBetaCat()) ) // assume all the cells are the same
    {
        mBetaCatResultsFile->close();
    }

    TissueSimulation<2>::AfterSolve();
}


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


void CryptSimulation2d::UseJiggledBottomCells()
{
    mUseJiggledBottomCells = true;
}

void CryptSimulation2d::ApplyTissueBoundaryConditions(TissueCell& rCell, ChastePoint<2>& rPoint)
{
    bool is_wnt_included = WntConcentration::Instance()->IsWntSetUp();

    if (!is_wnt_included) 
    {
        WntConcentration::Destroy();
        
        /**
         * If WntConcentration is not set up then stem cells must be pinned,
         * so we reset the x-coordinate of each stem cell.
         */
        if (rCell.GetCellType()==STEM)
        {
            unsigned index = rCell.GetLocationIndex();
            rPoint.rGetLocation()[0] = mrTissue.GetNode(index)->rGetLocation()[0];
            rPoint.rGetLocation()[1] = mrTissue.GetNode(index)->rGetLocation()[1];
        }
    }
    
    // Any cell that has moved below the bottom of the crypt must be moved back up
    if (rPoint.rGetLocation()[1] < 0.0)
    {
        rPoint.rGetLocation()[1] = 0.0;
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
            rPoint.rGetLocation()[1] = 0.05*mpRandomGenerator->ranf();
        }
    }
    assert(rPoint[1]>=0.0);    
}
