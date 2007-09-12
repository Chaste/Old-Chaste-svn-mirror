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
        double friction = stiffness;
        
        assert(body_force>0 && friction>0);

        for(Crypt<2>::Iterator iter = mrCrypt.Begin(); 
            iter != mrCrypt.End();
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
                drdt[iter->GetNodeIndex()](0) -= damping*friction;
                drdt[iter->GetNodeIndex()](1) += damping*stiffness*(mCellRadius - dist_to_bottom);
            }

            if( dist_to_top < mCellRadius )
            {
                drdt[iter->GetNodeIndex()](0) -= damping*friction;
                drdt[iter->GetNodeIndex()](1) -= damping*stiffness*(mCellRadius - dist_to_top);
            }
        } 
        
        return drdt;
    }    


    void UpdateNodePositions(const std::vector< c_vector<double, 2> >& rDrDt)
    {
        // update ghost positions first because they do not affect the real cells
        mrCrypt.UpdateGhostPositions(mDt);
    
        // Iterate over all cells to update their positions.
        for (Crypt<2>::Iterator cell_iter = mrCrypt.Begin();
             cell_iter != mrCrypt.End();
             ++cell_iter)
        {
            MeinekeCryptCell& cell = *cell_iter;
            unsigned index = cell.GetNodeIndex();
            
            ChastePoint<2> new_point(mrCrypt.rGetMesh().GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
                
            // for all cells - move up if below the bottom surface
            if (new_point.rGetLocation()[1] < mBottom)
            {
                new_point.rGetLocation()[1] = mBottom;
            }

            if (new_point.rGetLocation()[1] > mTop)
            {
                new_point.rGetLocation()[1] = mTop;
            }
            
            mrCrypt.MoveCell(cell_iter, new_point);                    
        }
    }
    
public:
    BloodCellModel(Crypt<2>& rCrypt, double bottom, double top)
        : TissueSimulation<2>(rCrypt)
    {
        assert(bottom < top);
        mBottom = bottom;
        mTop = top;
         
        for(Crypt<2>::Iterator iter = mrCrypt.Begin(); 
            iter != mrCrypt.End();
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
    void Run(std::string outputDirectory)
    { 
        assert(outputDirectory!="");
        
        int num_cells_depth = 6; 
        int num_cells_width = 20;
        unsigned num_ghosts = 0;
 
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, num_ghosts);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
 
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
    
        std::vector<MeinekeCryptCell> cells;

        for(unsigned j=0; j<p_mesh->GetNumNodes(); j++)
        {
            MeinekeCryptCell cell(DIFFERENTIATED, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = -10;
            cell.SetNodeIndex(j);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
    
        Crypt<2> crypt(*p_mesh,cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        double vessel_height = 10;
        double vessel_width  = 25;
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->SetCryptLength(vessel_height);
        p_params->SetCryptWidth(vessel_width);


        RandomNumberGenerator* rng = RandomNumberGenerator::Instance();
        for(Crypt<2>::Iterator iter = crypt.Begin(); 
            iter != crypt.End();
            ++iter)
        {
            // completely random
            double x = rng->ranf() * vessel_width;
            double y = rng->ranf() * vessel_height;
            
//// random perturbation
//            double x = iter.rGetLocation()[0] + (rng->ranf() * 0.2 - 0.1);
//            double y = iter.rGetLocation()[1] + (rng->ranf() * 0.2 - 0.1);

            ChastePoint<2> new_location(x,y);
            
            crypt.MoveCell(iter, new_location);
        }
    
        crypt.ReMesh();
 
        BloodCellModel blood_cell_model(crypt, 0, vessel_height);

        blood_cell_model.SetOutputDirectory(outputDirectory);
        blood_cell_model.SetEndTime(2);

        blood_cell_model.Solve();

        SimulationTime::Destroy();
    }
    


public:

    void TestNormal() throw(Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->Reset();
        
        Run("BloodCellModelNormal");
    }

    void TestTwiceAsStiff() throw(Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->Reset();
        double stiffness = p_params->GetSpringStiffness();
        p_params->SetSpringStiffness(2*stiffness);
        
        Run("BloodCellModelTwiceAsStiff");
    }

//    void xTestFiveTimesAsStiff() throw(Exception)
//    {
//        CancerParameters* p_params = CancerParameters::Instance();
//        p_params->Reset();
//        double stiffness = p_params->GetSpringStiffness();
//        p_params->SetSpringStiffness(2*stiffness);
//        
//        Run("BloodCellModelFiveTimesAsStiff");
//    }
};


#endif /*TESTBLOODCELLMODEL_HPP_*/
