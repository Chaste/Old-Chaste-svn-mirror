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
    double mTop;
    double mBottom;

    std::vector<c_vector<double,2> > CalculateVelocitiesOfEachNode()
    {
        std::vector<c_vector<double,2> > drdt = TissueSimulation<2>::CalculateVelocitiesOfEachNode();
        
        double damping = CancerParameters::Instance()->GetDampingConstantNormal();
        double body_force = 10;
        double friction = 5;
        
        assert(body_force>0 && friction>0);
        for(Crypt<2>::Iterator iter = mrCrypt.Begin(); 
            iter != mrCrypt.End();
            ++iter)
        {
            drdt[iter->GetNodeIndex()](0) += damping*body_force;
            
            double y = iter.rGetLocation()[1];
            if(y<mBottom + 0.1*(mTop - mBottom))
            {
                drdt[iter->GetNodeIndex()](0) -= damping*friction;
            }

            if(y>mBottom + 0.9*(mTop - mBottom))
            {
                drdt[iter->GetNodeIndex()](0) -= damping*friction;
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
            
            Point<2> new_point(mrCrypt.rGetMesh().GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
                
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
    BloodCellModel(Crypt<2>& rCrypt)
        : TissueSimulation<2>(rCrypt)
    {
        double min_y =  1e200;
        double max_y = -1e200;
         
        for(Crypt<2>::Iterator iter = mrCrypt.Begin(); 
            iter != mrCrypt.End();
            ++iter)
        {
            double y = iter.rGetLocation()[1];
 
            if(y > max_y)
            {
                max_y = y;
            }
            if(y < min_y)
            {
                min_y = y;
            }
            
        }
        
        mBottom = min_y - 0.5;
        mTop    = max_y + 0.5;
    }
};


class TestBloodCellModel : public CxxTest::TestSuite
{
public :
    void testBloodCellModel()
    { 
        int num_cells_depth = 6; 
        int num_cells_width = 20;
 
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u);
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

        BloodCellModel blood_cell_model(crypt);

        blood_cell_model.SetOutputDirectory("BloodCellModel");
        blood_cell_model.SetEndTime(2);

        blood_cell_model.Solve();
    }
};


#endif /*TESTBLOODCELLMODEL_HPP_*/
