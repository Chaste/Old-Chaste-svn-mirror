#ifndef CRYPTSIMULATION_HPP_
#define CRYPTSIMULATION_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
#include "StochasticCellCycleModel.hpp"
#include "Exception.hpp"
#include "ColumnDataWriter.hpp"

#include "MeinekeCryptCellTypes.hpp"

#include <cmath>
#include <ctime>
#include <iostream>

/**
 * Solve a crypt simulation based on the Meineke paper.
 *
 * The spring lengths are governed by the equations
 * dr/dt = stem_cycle_time*(mu/eta) sum_j r_hat_i,j*(|r_i,j|-s0)
 *       = alpha sum_j r_hat_i,j*(|r_i,j|-s0)
 *
 * where alpha = stem_cycle_time*(mu/eta) = stem_cycle_time*meineke_lambda.
 *       s0    = natural length of the spring

 * Length is scaled by natural length
 * Time is scaled by a stem cell cycle time

 * meineke_lambda = mu (spring constant) / eta (damping) = 0.01 (from Meineke - note
 * that the value we use for Meineke lambda is completely different because we have
 * nondimensionalised)
 */

class CryptSimulation
{
private:
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<1,1> &mrMesh;
    
    bool mIncludeRandomBirth;
    bool mIncludeVariableRestLength;
    
    unsigned mMaxCells;
    
    std::string mOutputDirectory;
    
    std::vector<MeinekeCryptCell> mCells;
    
    CancerParameters *mpParams;
    SimulationTime *mpSimulationTime;
    RandomNumberGenerator *mpGen;
    bool mCreatedRng;
    
public:

    /** Constructor
     *  @param rMesh
     *  @param cells is defaulted to the empty vector, in which case SetIncludeRandomBirth()
     *  should be called for any birth to happen.
     *  @param pGen is a RandomNumberGenerator class.  If it's not given then a new one is 
     *  constructed and random numbers are reseeded with srandom(0).
     */
    CryptSimulation(ConformingTetrahedralMesh<1,1> &rMesh,
                    std::vector<MeinekeCryptCell> cells = std::vector<MeinekeCryptCell>(),
                    RandomNumberGenerator *pGen = NULL)
            : mrMesh(rMesh),
            mCells(cells)
    {
    	mpSimulationTime = SimulationTime::Instance();
    	if (pGen!=NULL)
        {
            mpGen = pGen;
            mCreatedRng = false;
        }
        else
        {
            mpGen = new RandomNumberGenerator;
            mCreatedRng = true;
        }
        mpParams = CancerParameters::Instance();
        mDt = 1.0/(120.0); // ie 30 sec NOTE: hardcoded 120?
        mEndTime = 120.0; //hours
        
        mIncludeRandomBirth = false;
        mIncludeVariableRestLength = false;
        mOutputDirectory = "";
    }
    
    /**
     * Free any memory allocated by the constructor
     */
    ~CryptSimulation()
    {
        if (mCreatedRng)
        {
            delete mpGen;
        }
        SimulationTime::Destroy();
    }
    
    void SetDt(double dt)
    {
        assert(dt>0);
        mDt=dt;
    }
    
    void SetEndTime(double endTime)
    {
        assert(endTime>0);
        mEndTime=endTime;
    }
    
    
    /**
     *  Call this before Solve() if no cells have been specified. Randomly adds a new 
     *  node every 1 time unit, starting 0.1
     */
    void SetIncludeRandomBirth()
    {
        mIncludeRandomBirth = true;
    }
    
    void SetOutputDirectory(std::string outputDirectory)
    {
        mOutputDirectory = outputDirectory;
    }
    
    /**
     *  Call this before Solve() to simulate cell growth after cell division.
     *  (will eventually become SetIncludeCellBirth() and then become the default)
     */
    void SetIncludeVariableRestLength()
    {
        mIncludeVariableRestLength = true;
    }
    
    void SetMaxCells(unsigned maxCells)
    {
        mMaxCells = maxCells;
    }
    
    
    /** 
     *  Get the cells vector
     */
    std::vector<MeinekeCryptCell> GetCells()
    {
        assert(mCells.size()>0);
        return mCells;
    }   
    
    
    /**
     * Main Solve method.
     * 
     * Once CryptSimulation object has been set up, call this to run simulation
     */
    void Solve()
    {
    	if (mOutputDirectory=="")
        {
            EXCEPTION("OutputDirectory not set");
        }
        
        //Creating Column Data Writer Handler
        ColumnDataWriter tabulated_writer(mOutputDirectory, "tabulated_results");
        unsigned time_var_id = tabulated_writer.DefineUnlimitedDimension("Time","hours");
        
        std::vector<unsigned> type_var_ids;
        std::vector<unsigned> position_var_ids;
        
        type_var_ids.resize(mMaxCells);
        position_var_ids.resize(mMaxCells);
        
        // set up columns
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
        mpSimulationTime = SimulationTime::Instance();
        mpSimulationTime->SetStartTime(0.0);
        mpSimulationTime->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
                                     
        //double time = 0.0;
        double time_since_last_birth = 15.0;//15 hours - only used in non-random birth
        
        unsigned num_births = 0;
        unsigned num_deaths = 0;
        
        std::vector<double> new_point_position(mrMesh.GetNumAllNodes());
        
        // Creating Simple File Handler
        OutputFileHandler output_file_handler(mOutputDirectory, false);
        out_stream p_results_file = output_file_handler.OpenOutputFile("results");
        while ( mpSimulationTime->GetTimeStepsElapsed() < num_time_steps)
        {
            //std::cout << "Simulation time = " << mpSimulationTime->GetDimensionalisedTime() << "\n" << std::endl;
            // Cell birth
            if (mIncludeRandomBirth && time_since_last_birth > 5.0)// Ie Birth every 5 hours
            {
                unsigned new_node_index = AddRandomNode(mpSimulationTime->GetDimensionalisedTime());
                time_since_last_birth = 0 ;
                // Create new cell note all are Stem Cells and have generation 0 for random birth
                //RandomNumberGenerator *pGen=new RandomNumberGenerator;
                CryptCellType cell_type=STEM ;
                unsigned generation=0;
                
                MeinekeCryptCell new_cell(cell_type, HEALTHY, generation, new StochasticCellCycleModel(mpGen));

                // Update cells vector
                new_cell.SetNodeIndex(new_node_index);
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
            else if (!mCells.empty())
            {
                for (unsigned i=0; i<mCells.size(); i++)
                {
                    if (mrMesh.GetNode(i)->IsDeleted()) continue; // Skip deleted cells
                    // Check for this cell dividing
                    if (mCells[i].ReadyToDivide())
                    {
                        // Create new cell
                        MeinekeCryptCell new_cell = mCells[i].Divide();
                        
                        // Add new node to mesh
                        Node<1> *p_our_node = mrMesh.GetNode(i);
                        
                        //Note: May need to check which side element is put esp. at the ends
                        Element<1,1> *p_element = mrMesh.GetElement(p_our_node->GetNextContainingElementIndex());
                        
                        unsigned new_node_index = AddNodeToElement(p_element,mpSimulationTime->GetDimensionalisedTime());
                        // Update cells        	 variableID unknown vector
                        new_cell.SetNodeIndex(new_node_index);
                        if (new_node_index == mCells.size())
                        {
                            mCells.push_back(new_cell);
                        }
                        else
                        {
                            mCells[new_node_index] = new_cell;
                        }
                        num_births++;
                        //std::cout<< "num_births=" << num_births <<std::endl<< std::flush;
                    }
                }
            }
            
            // calculate node velocities
            std::vector<double> drdt(mrMesh.GetNumAllNodes());
            if (mIncludeVariableRestLength && !mCells.empty())
            {
                //std::cout<< "elements" << mrMesh.GetNumAllElements() <<std::endl<< std::flush;
                
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
                        drdt_contributions(0) = mpParams->GetMeinekeLambda() *(  unit_vector_forward  * (distance_between_nodes - rest_length) );
                        drdt_contributions(1) = mpParams->GetMeinekeLambda() *(  unit_vector_backward * (distance_between_nodes - rest_length) );
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
                        
                        drdt_contributions(0) = mpParams->GetMeinekeLambda() *(  unit_vector_forward  * (distance_between_nodes - 1.0) );
                        drdt_contributions(1) = mpParams->GetMeinekeLambda() *(  unit_vector_backward * (distance_between_nodes - 1.0) );
                        
                        drdt[ element->GetNode(0)->GetIndex() ] += drdt_contributions(0);
                        drdt[ element->GetNode(1)->GetIndex() ] += drdt_contributions(1);
                    }
                }
            }
            
            // update node positions
            for (unsigned index = 1; index<mrMesh.GetNumAllNodes(); index++)
            {
                // assume stem cells are fixed, or if there are no cells, fix node 0
                //if(
                if (!mrMesh.GetNode(index)->IsDeleted())
                {
                    c_vector<double, 1> node_loc = mrMesh.GetNode(index)->rGetLocation();
                    Point<1> new_point;
                    new_point.rGetLocation()[0] = node_loc[0] + mDt*drdt[index]; // new_point_position[index];
                    mrMesh.SetNode(index, new_point, false);
                }
            }
            
            // Remove nodes that are beyond the crypt
            while (true)
            {
                bool sloughed_node = false;
                ConformingTetrahedralMesh<1,1>::BoundaryNodeIterator it = mrMesh.GetBoundaryNodeIteratorEnd();
                while (it != mrMesh.GetBoundaryNodeIteratorBegin())
                {
                    it--;
                    const Node<1> *p_node = *it;
                    if (p_node->rGetLocation()[0] > mpParams->GetCryptLength())
                    {
                        // It's fallen off
                        mrMesh.DeleteBoundaryNodeAt(p_node->GetIndex());
                        num_deaths++;
                        //std::cout<< "num_deaths=" << num_deaths <<std::endl<< std::flush;
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
                        CryptCellType type  = mCells[index].GetCellType();
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
#define COVERAGE_IGNORE
                            assert(0);
#undef COVERAGE_IGNORE
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
    
    
    
private:
    unsigned AddRandomNode(double time)
    {
    
        //Pick an element
        unsigned random_element_number = mpGen->randMod(mrMesh.GetNumAllElements());
        Element<1,1>* p_random_element = mrMesh.GetElement(random_element_number);
        double element_length = fabs(p_random_element->GetNodeLocation(1,0) - p_random_element->GetNodeLocation(0,0));
        //std::cout << "length " <<element_length << "\n";
        
        // keep picking until find an element which is big enough and not deleted
        while (element_length <0.4 || p_random_element->IsDeleted())
        {
            // the following is ignored in coverage as there is a random
            // chance of it not happening
#define COVERAGE_IGNORE
            random_element_number = mpGen->randMod(mrMesh.GetNumAllElements());
            p_random_element = mrMesh.GetElement(random_element_number);
            element_length = fabs(p_random_element->GetNodeLocation(1,0) - p_random_element->GetNodeLocation(0,0));
#undef COVERAGE_IGNORE
            //std::cout << "..too small, trying: length " <<element_length << "\n";
            //double random_displacement = 0.2+mpGen->randf()*0.6;
        }
        // Reset age of left node to zero
        mCells[p_random_element->GetNode(0)->GetIndex()].SetBirthTime(time);
        return AddNodeToElement(p_random_element,time);
    }
    
    
    unsigned AddNodeToElement(Element<1,1>* pElement, double time)
    {
    
        double displacement;
        double left_position= pElement->GetNodeLocation(0,0);
        if (mIncludeVariableRestLength)
        {
            double age0 = mCells[pElement->GetNode(0)->GetIndex()].GetAge();
            double age1 = mCells[pElement->GetNode(1)->GetIndex()].GetAge();
            
            if (fabs(age0)<1e-6)
            {
                // place the new node to 0.1 to the right of the left-hand node
                displacement = 0.1;
            }
            else if (fabs(age1)<1e-6)
            {
                // place the new node to 0.1 to the left of the right-hand node
                double element_length = fabs(pElement->GetNodeLocation(1,0) - pElement->GetNodeLocation(0,0));
                displacement = element_length - 0.1;
            }
            else
            {
                #define COVERAGE_IGNORE
                EXCEPTION("No cell has divided in this element");
                #undef COVERAGE_IGNORE
            }
        }
        else
        {
            double element_length = fabs(pElement->GetNodeLocation(1,0) - pElement->GetNodeLocation(0,0));
            // pick a random position in the central 60% of the element
            displacement = 0.2 + (mpGen->ranf())*(element_length-0.4);
            
        }
        Point<1> new_point(left_position + displacement);
        
        return mrMesh.RefineElement(pElement, new_point);
    }
    
    
     
};

#endif /*CRYPTSIMULATION_HPP_*/
