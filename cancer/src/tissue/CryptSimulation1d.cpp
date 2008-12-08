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

#include "CryptSimulation1d.hpp"

/** Constructor
 *  @param rMesh
 *  @param cells is defaulted to the empty vector, in which case SetIncludeRandomBirth()
 *  should be called for any birth to happen.
 *  @param pGen is a RandomNumberGenerator class.  If it's not given then a new one is
 *  constructed and random numbers are reseeded with srandom(0).
 */
CryptSimulation1d::CryptSimulation1d(MutableMesh<1,1> &rMesh,
                                     std::vector<TissueCell> cells)
        : mrMesh(rMesh),
        mCells(cells)
{
    mpSimulationTime = SimulationTime::Instance();
    mpGen = RandomNumberGenerator::Instance();
    mpParams = CancerParameters::Instance();
    mpParams->SetSpringStiffness(30.0);
    mDt = 1.0/(120.0); // Timestep of 30 seconds (as per Meineke)
    mEndTime = 120.0; // hours

    mIncludeVariableRestLength = false;
    mOutputDirectory = "";

    mpSimulationTime = SimulationTime::Instance();
    if (!mpSimulationTime->IsStartTimeSetUp())
    {
        EXCEPTION("Start time not set in simulation time singleton object");
    }
}

/**
 * Free any memory allocated by the constructor
 */
CryptSimulation1d::~CryptSimulation1d()
{
    SimulationTime::Destroy();
}

/**
 * Change the time step used for mechanics
 */
void CryptSimulation1d::SetDt(double dt)
{
    assert(dt>0);
    mDt=dt;
}

/**
 * Set the end time of the next simualtion.Solve() call.
 */
void CryptSimulation1d::SetEndTime(double endTime)
{
    assert(endTime>0);
    mEndTime=endTime;
}

/**
 * The subdirectory of CHASTETESTOUTPUT where the results will be stored.
 */
void CryptSimulation1d::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
}

/**
 *  Call this before Solve() to simulate cell growth after cell division.
 *  (will eventually become SetIncludeCellBirth() and then become the default)
 */
void CryptSimulation1d::SetIncludeVariableRestLength()
{
    mIncludeVariableRestLength = true;
}

/**
 * The maximum number of cells that will be in this simulation.
 */
void CryptSimulation1d::SetMaxCells(unsigned maxCells)
{
    mMaxCells = maxCells;
}

/**
 *  Get the cells vector
 */
std::vector<TissueCell> CryptSimulation1d::GetCells()
{
    assert(mCells.size()>0);
    return mCells;
}

/**
 * Main Solve method.
 *
 * Once CryptSimulation1d object has been set up, call this to run simulation
 */
void CryptSimulation1d::Solve()
{
    if (mOutputDirectory=="")
    {
        EXCEPTION("OutputDirectory not set");
    }

    // Creating Column Data Writer Handler
    ColumnDataWriter tabulated_writer(mOutputDirectory, "tabulated_results");
    unsigned time_var_id = tabulated_writer.DefineUnlimitedDimension("Time","hours");

    std::vector<unsigned> type_var_ids;
    std::vector<unsigned> position_var_ids;

    type_var_ids.resize(mMaxCells);
    position_var_ids.resize(mMaxCells);

    // Set up columns
    for (unsigned cell=0; cell<mMaxCells; cell++)
    {
        std::stringstream cell_type_var_name;
        std::stringstream cell_position_var_name;
        cell_type_var_name << "cell_type_" << cell;
        cell_position_var_name << "cell_position_" << cell;
        type_var_ids[cell]=tabulated_writer.DefineVariable(cell_type_var_name.str(),"dimensionless");
        position_var_ids[cell]=tabulated_writer.DefineVariable(cell_position_var_name.str(),"cell_lengths");
    }
    tabulated_writer.EndDefineMode();

    unsigned num_time_steps = (unsigned)(mEndTime/mDt+0.5);
    mpSimulationTime->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);

    double time_since_last_birth = 15.0; // 15 hours - only used in non-random birth

    unsigned num_births = 0;
    unsigned num_deaths = 0;

    std::vector<double> new_point_position(mrMesh.GetNumAllNodes());

    // Creating Simple File Handler
    OutputFileHandler output_file_handler(mOutputDirectory, false);
    out_stream p_results_file = output_file_handler.OpenOutputFile("results");
    while ( mpSimulationTime->GetTimeStepsElapsed() < num_time_steps)
    {
        // Cell birth
        if (!mCells.empty())
        {
            for (unsigned i=0; i<mCells.size(); i++)
            {
                if (mrMesh.GetNode(i)->IsDeleted()) continue; // Skip deleted cells
                // Check for this cell dividing
                if (mCells[i].ReadyToDivide())
                {
                    // Create new cell
                    TissueCell new_cell = mCells[i].Divide();

                    // Add new node to mesh
                    Node<1> *p_our_node = mrMesh.GetNode(i);

                    // Note: May need to check which side element is put esp. at the ends
                    Element<1,1> *p_element = mrMesh.GetElement(*(p_our_node->ContainingElementsBegin()));

                    unsigned new_node_index = AddNodeToElement(p_element,mpSimulationTime->GetDimensionalisedTime());

                    // Update cells
                    new_cell.SetLocationIndex(new_node_index);
                    if (new_node_index == mCells.size())
                    {
                        mCells.push_back(new_cell);
                    }
                    else
                    {
                        mCells[new_node_index] = new_cell;
                    }
                    num_births++;
                }
            }
        }

        // Calculate node velocities
        std::vector<double> drdt(mrMesh.GetNumAllNodes());
        if (mIncludeVariableRestLength && !mCells.empty())
        {
            for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
            {
                Element<1,1>* element = mrMesh.GetElement(elem_index);
                if (!element->IsDeleted())
                {
                    c_vector<double, 2> drdt_contributions;
                    double distance_between_nodes = fabs(element->GetNodeLocation(1,0) - element->GetNodeLocation(0,0));
                    double unit_vector_backward = -1;
                    double unit_vector_forward = 1;
                    double age0 = mCells[element->GetNode(0)->GetIndex()].GetAge();
                    double age1 = mCells[element->GetNode(1)->GetIndex()].GetAge();
                    double rest_length = 1.0;

                    if (age0<1.0 && age1<1.0 && fabs(age0-age1)<1e-6)
                    {
                        /* Spring Rest Length Increases to normal rest length from 0.9 to normal rest length, 1.0, over 1 hour
                         * This doesnt happen at present as when the full line is included the tests fail
                         *
                         * This is wrong but due to the model being set up in 1D, when a new cell with a weaker spring is
                         * put in next to other stressed cells, the weaker spring will be compressed too much and lead to
                         * cells being pushed through other ones.  Leading to an exception being thrown in line 319 ish.
                         */
                        rest_length=(0.9+0.1*age0);

                        assert(rest_length<=1.0);
                    }
                    drdt_contributions(0) = ( mpParams->GetSpringStiffness() / mpParams->GetDampingConstantNormal() ) *(  unit_vector_forward  * (distance_between_nodes - rest_length) );
                    drdt_contributions(1) = ( mpParams->GetSpringStiffness() / mpParams->GetDampingConstantNormal() ) *(  unit_vector_backward * (distance_between_nodes - rest_length) );
                    drdt[ element->GetNode(0)->GetIndex() ] += drdt_contributions(0);
                    drdt[ element->GetNode(1)->GetIndex() ] += drdt_contributions(1);
                }
            }
        }
        else
        {
            for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
            {
                Element<1,1>* element = mrMesh.GetElement(elem_index);
                if (!element->IsDeleted())
                {
                    c_vector<double, 2> drdt_contributions;
                    double distance_between_nodes = fabs(element->GetNodeLocation(1,0) - element->GetNodeLocation(0,0));
                    double unit_vector_backward = -1;
                    double unit_vector_forward = 1;

                    drdt_contributions(0) =( mpParams->GetSpringStiffness() / mpParams->GetDampingConstantNormal() ) *(  unit_vector_forward  * (distance_between_nodes - 1.0) );
                    drdt_contributions(1) =( mpParams->GetSpringStiffness() / mpParams->GetDampingConstantNormal() ) *(  unit_vector_backward * (distance_between_nodes - 1.0) );

                    drdt[ element->GetNode(0)->GetIndex() ] += drdt_contributions(0);
                    drdt[ element->GetNode(1)->GetIndex() ] += drdt_contributions(1);
                }
            }
        }

        //Update node positions
        for (unsigned index = 1; index<mrMesh.GetNumAllNodes(); index++)
        {
            // assume stem cells are fixed, or if there are no cells, fix node 0
            if (!mrMesh.GetNode(index)->IsDeleted())
            {
                c_vector<double, 1> node_loc = mrMesh.GetNode(index)->rGetLocation();
                ChastePoint<1> new_point;
                new_point.rGetLocation()[0] = node_loc[0] + mDt*drdt[index]; // new_point_position[index];
                mrMesh.SetNode(index, new_point, false);
            }
        }

        // Remove nodes that are beyond the crypt
        while (true)
        {
            bool sloughed_node = false;
            MutableMesh<1,1>::BoundaryNodeIterator it = mrMesh.GetBoundaryNodeIteratorEnd();
            while (it != mrMesh.GetBoundaryNodeIteratorBegin())
            {
                it--;
                const Node<1> *p_node = *it;
                if (p_node->rGetLocation()[0] > mpParams->GetCryptLength())
                {
                    // It's fallen off
                    mrMesh.DeleteBoundaryNodeAt(p_node->GetIndex());
                    num_deaths++;
                    sloughed_node = true;
                    break;
                }
            }
            if (!sloughed_node) break;
        }
        // Check nodes havent crossed
        mrMesh.RefreshMesh();

        // Increment simulation time here, so results files look sensible
        mpSimulationTime->IncrementTimeOneStep();

        // Writing Results To Tabulated File First And Then To Space Separated File
        tabulated_writer.PutVariable(time_var_id, mpSimulationTime->GetDimensionalisedTime());
        (*p_results_file) << mpSimulationTime->GetDimensionalisedTime() << "\t";

        unsigned cell=0; // NB this is not the index in mCells, but the index in the mesh!
        for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
        {
            if (!mrMesh.GetNode(index)->IsDeleted())
            {
                if (mCells.size() > 0)
                {
                    CellType type  = mCells[index].GetCellType();
                    if (type == STEM)
                    {
                        tabulated_writer.PutVariable(type_var_ids[cell], 0);
                    }
                    else if (type == TRANSIT)
                    {
                        tabulated_writer.PutVariable(type_var_ids[cell], 1);
                    }
                    else if (type == DIFFERENTIATED)
                    {
                        tabulated_writer.PutVariable(type_var_ids[cell], 2);
                    }
                    else
                    {
                        // should be impossible to get here, until cancer cells
                        // are implemented
                        NEVER_REACHED;
                    }
                }
                else
                {
                    tabulated_writer.PutVariable(type_var_ids[cell], -1);
                }

                const c_vector<double, 1> node_loc = mrMesh.GetNode(index)->rGetLocation();
                tabulated_writer.PutVariable(position_var_ids[cell], node_loc[0]);
                (*p_results_file) << node_loc[0] << " ";

                cell++;
            }
        }
        tabulated_writer.AdvanceAlongUnlimitedDimension();
        (*p_results_file) << "\n";

        time_since_last_birth += mDt;
    }
}

/**
 * This handles the birth of a cell in 1D by splitting an existing element in two and placing
 * a new point in the middle of it.
 */
unsigned CryptSimulation1d::AddNodeToElement(Element<1,1>* pElement, double time)
{

    double displacement;
    double left_position= pElement->GetNodeLocation(0,0);
    if (mIncludeVariableRestLength)
    {
        double age0 = mCells[pElement->GetNode(0)->GetIndex()].GetAge();
        double age1 = mCells[pElement->GetNode(1)->GetIndex()].GetAge();

        if (fabs(age0)<1e-6)
        {
            // Place the new node to 0.1 to the right of the left-hand node
            displacement = 0.1;
        }
        else if (fabs(age1)<1e-6)
        {
            // Place the new node to 0.1 to the left of the right-hand node
            double element_length = fabs(pElement->GetNodeLocation(1,0) - pElement->GetNodeLocation(0,0));
            displacement = element_length - 0.1;
        }
        else
        {
            // This called by Tyson Novak cells which might not be age 0 when the simulation divides them.
            double element_length = fabs(pElement->GetNodeLocation(1,0) - pElement->GetNodeLocation(0,0));
            // Pick a random position in the central 60% of the element
            displacement = 0.2 + (mpGen->ranf())*(element_length-0.4);
        }
    }
    else
    {
        double element_length = fabs(pElement->GetNodeLocation(1,0) - pElement->GetNodeLocation(0,0));
        // Pick a random position in the central 60% of the element
        displacement = 0.2 + (mpGen->ranf())*(element_length-0.4);

    }
    ChastePoint<1> new_point(left_position + displacement);

    return mrMesh.RefineElement(pElement, new_point);
}

