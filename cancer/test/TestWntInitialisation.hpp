#ifndef TESTWNTINITIALISATION_HPP_
#define TESTWNTINITIALISATION_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"
#include "WntGradient.hpp"
#include "WntGradientTypes.hpp"
#include "Crypt.cpp"
#include "MeinekeCryptCell.hpp"
#include "WntCellCycleModel.hpp"

class TestWntInitialisation : public CxxTest::TestSuite
{
public:
    void TestWntInitialisationSetup() throw(Exception)
    {
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        WntGradient wnt_gradient(LINEAR);
        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        std::vector<WntCellCycleModel*> models;
        
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel(0.0,wnt_gradient);
            p_model->SetUseWntGradient();
            MeinekeCryptCell cell(STEM, HEALTHY, 0, p_model);
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);

            cells.push_back(cell);
            models.push_back(p_model);
        }
        
        // create the crypt
        Crypt<2> crypt(mesh,cells);
        
        CancerParameters::Instance()->SetCryptLength(1.0);

        //wnt_gradient.SetCrypt(crypt);
        SingletonWntGradient::Instance()->SetType(LINEAR);
        SingletonWntGradient::Instance()->SetCrypt(crypt);
        
        Crypt<2>::Iterator iter = crypt.Begin();
        

        while(iter!=crypt.End())
        {
            const WntCellCycleModel* p_model = (WntCellCycleModel*) iter->GetCellCycleModel();
            std::vector<double> proteins = p_model->GetProteinConcentrations();
        
            if(iter.GetNode()->rGetLocation()[1]==0.0)
            {
                TS_ASSERT_DELTA(proteins[5],4.975124378109454e-03, 1e-3);
                TS_ASSERT_DELTA(proteins[6]+proteins[7],6.002649406788524e-01, 1e-3);
                TS_ASSERT_DELTA(proteins[8],1.00, 1e-3);
            }
            else
            {
                TS_ASSERT_DELTA(iter.GetNode()->rGetLocation()[1],1.0,1e-12);
                TS_ASSERT_DELTA(proteins[5],1.000, 1e-3);
                TS_ASSERT_DELTA(proteins[6]+proteins[7],0.0074, 1e-3);
                TS_ASSERT_DELTA(proteins[8],0.00, 1e-3);
            }

            ++iter;
        }

        SingletonWntGradient::Destroy();
    }
};
#endif /*TESTWNTINITIALISATION_HPP_*/
