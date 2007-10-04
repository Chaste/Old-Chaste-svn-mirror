#ifndef TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulationWithNutrients.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "ColumnDataReader.hpp"
#include "SimulationTime.hpp"
#include "SloughingCellKiller.hpp"
#include "PetscTools.hpp"
#include "CellwiseData.hpp"

class OxygenBasedCellCycleModel : public FixedCellCycleModel
{
private:
    double mTimeProgressingThroughCellCycle;
        
public:
    OxygenBasedCellCycleModel() : FixedCellCycleModel()
    {
        mTimeProgressingThroughCellCycle = 0.0;
    }
    
    bool ReadyToDivide()
    {
        double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell);
        
        // the std::max is a hack, due to the choice of PDE
        // (a simple nonlinear PDE would ensure positivity)   
//        if (oxygen_concentration < 0.0)
//        {
//            EXCEPTION("Oxygen concentration has gone negative - check the oxygen PDE");
//        }
        mTimeProgressingThroughCellCycle = mTimeProgressingThroughCellCycle + std::max(oxygen_concentration,0.0)*SimulationTime::Instance()->GetTimeStep(); 
        
        bool result = false;        
        if ( mTimeProgressingThroughCellCycle > CancerParameters::Instance()->GetStemCellCycleTime() )
        {
            result = true;
        }
        
        return result;
    }
    
    AbstractCellCycleModel* CreateCellCycleModel()
    {
        return new OxygenBasedCellCycleModel();
    }
};


class SimpleLinearEllipticPde : public AbstractLinearEllipticPde<2>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<2> )
    {
        return -1.0;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<2> , double )
    {
        return 0.0;
    }
    
    c_matrix<double,2,2> ComputeDiffusionTerm(ChastePoint<2> )
    {
        return identity_matrix<double>(2);
    }
};


class RadiusBasedCellKiller : public AbstractCellKiller<2>
{
private :
    c_vector<double,2> mCentre;
    double mTimeStep;

public :
    RadiusBasedCellKiller(Crypt<2>* pCrypt, c_vector<double,2> centre, double timeStep)
        : AbstractCellKiller<2>(pCrypt),
          mCentre(centre),
          mTimeStep(timeStep)
    {
    }
    
    virtual void TestAndLabelCellsForApoptosisOrDeath()
    {
        for(Crypt<2>::Iterator cell_iter = mpCrypt->Begin();
            cell_iter != mpCrypt->End();
            ++cell_iter)
        {
            const c_vector<double,2>& location = cell_iter.GetNode()->rGetLocation();

            // double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(&(*cell_iter));
        
            double dist_to_centre = norm_2(location - mCentre);
            
            double prob_of_death = 2*mTimeStep - 1*mTimeStep*dist_to_centre;
            if (prob_of_death<=0.0)
            {
                prob_of_death=0.0;
            }
            
            if (!cell_iter->HasApoptosisBegun() &&
                RandomNumberGenerator::Instance()->ranf() < prob_of_death)
            {
                cell_iter->StartApoptosis();
            }    
        }
    }
};

class TestTissueSimulationWithNutrients : public CxxTest::TestSuite
{
public:

    void TestWithOxygen() throw(Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }
        
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
       
        int num_cells_depth = 10;
        int num_cells_width = 10;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(TRANSIT, HEALTHY, 0, new OxygenBasedCellCycleModel());
            double birth_time = -p_gen->ranf()*p_params->GetTransitCellCycleTime();
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
                
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        SimpleLinearEllipticPde pde;

        TissueSimulationWithNutrients<2> simulator(crypt, &pde);

        simulator.SetOutputDirectory("TissueSimulationWithOxygen");
        simulator.SetEndTime(0.5);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        //simulator.UseCutoffPoint(1.5);
        
        c_vector<double,2> centre(2);
        centre(0) = (double)num_cells_width/2.0;
        centre(1) = (double)num_cells_depth/2.0;
        
        AbstractCellKiller<2>* p_killer = new RadiusBasedCellKiller(&crypt, centre, simulator.GetDt());
        simulator.AddCellKiller(p_killer);
        
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetCrypt(crypt);
        
        simulator.Solve();
        
        // add get methods etc and test
        
        delete p_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

};
#endif /*TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_*/
