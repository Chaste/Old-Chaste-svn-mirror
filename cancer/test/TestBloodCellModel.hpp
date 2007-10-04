#ifndef TESTBLOODCELLMODEL_HPP_
#define BLOODCELLMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"
#include <cmath>
#include <vector>
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "HoneycombMeshGenerator.hpp"



class BloodCellModel : public TissueSimulation<2>
{
private:
    static const double mCellRadius = 0.5;

    double mTop;
    double mBottom;
    
    std::vector<c_vector<double,2> > CalculateVelocitiesOfEachNode()
    {
        std::vector<c_vector<double,2> > drdt = TissueSimulation<2>::CalculateVelocitiesOfEachNode();
        
        double damping = CancerParameters::Instance()->GetDampingConstantNormal();
        double stiffness = CancerParameters::Instance()->GetSpringStiffness();
        double body_force = 10;
        double coeff_of_friction = 2;
        
        assert(body_force>0 && coeff_of_friction>0);

        double average_flow_rate = 0;

        for(Tissue<2>::Iterator iter = mrTissue.Begin(); 
            iter != mrTissue.End();
            ++iter)
        {
            // body force
            drdt[iter->GetNodeIndex()](0) += damping*body_force;
            
            // friction
            double y = iter.rGetLocation()[1];
            double dist_to_bottom = y - mBottom;
            double dist_to_top    = mTop - y; 
            
            if( dist_to_bottom < mCellRadius )
            {
                double reaction = stiffness*(mCellRadius - dist_to_bottom);
                
                drdt[iter->GetNodeIndex()](1) += damping * reaction; 
                drdt[iter->GetNodeIndex()](0) -= coeff_of_friction * reaction * damping;
                if(drdt[iter->GetNodeIndex()](0) < 0)
                {
                    drdt[iter->GetNodeIndex()](0) = 0;
                }
            }

            if( dist_to_top < mCellRadius )
            {
                double reaction = stiffness*(mCellRadius - dist_to_top);

                drdt[iter->GetNodeIndex()](1) -= damping * reaction; 
                drdt[iter->GetNodeIndex()](0) -= coeff_of_friction * reaction * damping;
                if(drdt[iter->GetNodeIndex()](0) < 0)
                {
                    drdt[iter->GetNodeIndex()](0) = 0;
                }
            }
            
            average_flow_rate += drdt[iter->GetNodeIndex()](0);
        } 
        
        average_flow_rate /= mrTissue.GetNumRealCells();
        
        std::cout << SimulationTime::Instance()->GetDimensionalisedTime() << ": " << average_flow_rate <<"\n";  
        
        return drdt;
    }    


    void UpdateNodePositions(const std::vector< c_vector<double, 2> >& rDrDt)
    {
        // update ghost positions first because they do not affect the real cells
        mrTissue.UpdateGhostPositions(mDt);
    
        // Iterate over all cells to update their positions.
        for (Tissue<2>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            TissueCell& cell = *cell_iter;
            unsigned index = cell.GetNodeIndex();
            
            ChastePoint<2> new_point(mrTissue.rGetMesh().GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
                
            // for all cells - move up if below the bottom surface
            if (new_point.rGetLocation()[1] < mBottom)
            {
                new_point.rGetLocation()[1] = mBottom;
            }

            if (new_point.rGetLocation()[1] > mTop)
            {
                new_point.rGetLocation()[1] = mTop;
            }
            
            mrTissue.MoveCell(cell_iter, new_point);                    
        }
    }
    
    void PostSolve()
    {
        double stiffness_ratio = 4;
        static bool has_stiffened = false;
        if(SimulationTime::Instance()->GetDimensionalisedTime() > 0.5*mEndTime)
        {
            if(!has_stiffened)
            {
                CancerParameters* p_params = CancerParameters::Instance();
                double stiffness = p_params->GetSpringStiffness();
                p_params->SetSpringStiffness(stiffness_ratio*stiffness);
                has_stiffened = true;
                
                std::cout << "Stiffening from " << stiffness << " to " << stiffness_ratio*stiffness 
                          << " at t=" << SimulationTime::Instance()->GetDimensionalisedTime() << "\n";
            }
        }
    }
    
public:
    BloodCellModel(Tissue<2>& rCrypt, double bottom, double top)
        : TissueSimulation<2>(rCrypt)
    {
        assert(bottom < top);
        mBottom = bottom;
        mTop = top;
         
        for(Tissue<2>::Iterator iter = mrTissue.Begin(); 
            iter != mrTissue.End();
            ++iter)
        {
            double y = iter.rGetLocation()[1];
            assert((y>mBottom) && (y<mTop));
        }

        UseCutoffPoint(1.0);
    }
};



class TestBloodCellModel : public CxxTest::TestSuite
{
private:
    void Run(std::string outputDirectory, double vesselHeight = 10)
    { 
        assert(outputDirectory!="");
        
        int num_cells_depth = 3; 
        int num_cells_width = 10;
        unsigned num_ghosts = 0;
 
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, num_ghosts);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
 
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
    
        std::vector<TissueCell> cells;

        for(unsigned j=0; j<p_mesh->GetNumNodes(); j++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = -10;
            cell.SetNodeIndex(j);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
    
        Tissue<2> crypt(*p_mesh,cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        double vessel_width  = 20;
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->SetCryptLength(vesselHeight);
        p_params->SetCryptWidth(vessel_width);


        RandomNumberGenerator* rng = RandomNumberGenerator::Instance();
        for(Tissue<2>::Iterator iter = crypt.Begin(); 
            iter != crypt.End();
            ++iter)
        {
            // completely random
            double x = rng->ranf() * vessel_width;
            double y = rng->ranf() * vesselHeight;
            
//// random perturbation
//            double x = iter.rGetLocation()[0] + (rng->ranf() * 0.2 - 0.1);
//            double y = iter.rGetLocation()[1] + (rng->ranf() * 0.2 - 0.1);

            ChastePoint<2> new_location(x,y);
            
            crypt.MoveCell(iter, new_location);
        }
    
        crypt.ReMesh();
 
        BloodCellModel blood_cell_model(crypt, 0, vesselHeight);

        blood_cell_model.SetOutputDirectory(outputDirectory);
        blood_cell_model.SetEndTime(5);

        blood_cell_model.Solve();

        SimulationTime::Destroy();
    }
    


public:

    void xTestNormal() throw(Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->Reset();
        
        Run("BloodCellModelNormal");
    }

    void xTestTwiceAsStiff() throw(Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->Reset();
        double stiffness = p_params->GetSpringStiffness();
        p_params->SetSpringStiffness(2*stiffness);
        
        Run("BloodCellModelTwiceAsStiff");
    }
    
    void TestMe() throw(Exception)
    {
        double stiffness_ratio = 1;
        double height = 3;
        
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->Reset();
        double stiffness = p_params->GetSpringStiffness();
        p_params->SetSpringStiffness(stiffness_ratio*stiffness);
        
        std::stringstream ss;
        ss << "BloodCellModel_" << stiffness_ratio << "_" << height;
        
        Run(ss.str(), height);
    }
};


#endif /*TESTBLOODCELLMODEL_HPP_*/
