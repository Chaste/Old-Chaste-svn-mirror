#ifndef TESTWNTGRADIENT_HPP_
#define TESTWNTGRADIENT_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"
#include "WntGradient.hpp"
#include "Tissue.cpp"
#include "TissueCell.hpp"
#include "WntCellCycleModel.hpp"

class TestWntGradient : public CxxTest::TestSuite
{
public:    
    void TestNoWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        WntGradient* p_wnt_gradient = WntGradient::Instance();
        p_wnt_gradient->SetType(NONE);
        
        TS_ASSERT_EQUALS(p_wnt_gradient->GetType(), NONE);
        
        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        WntGradient::Destroy();
    }
    
    void TestLinearWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        WntGradient* p_wnt_gradient = WntGradient::Instance();
        p_wnt_gradient->SetType(LINEAR);
        
        CancerParameters *params = CancerParameters::Instance();
                
        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 1.0-height/params->GetCryptLength(), 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        
        WntGradient::Destroy();
    }
    
    
    void TestOffsetLinearWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        WntGradient* p_wnt_gradient = WntGradient::Instance();
        p_wnt_gradient->SetType(OFFSET_LINEAR);

        CancerParameters *params = CancerParameters::Instance();
        
        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0 , 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        // under a third of the way up the crypt.
        params->SetCryptLength(22.0);
        height = 7.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 1.0 - height/((1.0/3.0)*params->GetCryptLength()) , 1e-9);
        // more than a third of the way up the crypt.
        height = 10.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        WntGradient::Destroy();
    }
    
    
    void TestRadialWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        WntGradient* p_wnt_gradient = WntGradient::Instance();
        p_wnt_gradient->SetType(RADIAL);

        CancerParameters *params = CancerParameters::Instance();
        
        // Test GetWntLevel(double) method
        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-4);
        
        height = -1e-12;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-4);        
        
        height = 21.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0454 , 1e-4);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        
        params->SetCryptLength(22.0);
        height = 7.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.6818, 1e-4);
        
        
        // Test GetWntLevel(TissueCell*) method
        
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // translate mesh so that its centre is at (0,0)
        mesh.Translate(-0.5,-0.5); 
        
        std::vector<TissueCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_model);
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);        
        }
        
        // create a crypt
        Tissue<2> crypt(mesh,cells);        
        CancerParameters::Instance()->SetCryptLength(1.0);
        p_wnt_gradient->SetTissue(crypt);
        
        Tissue<2>::Iterator cell_iter = crypt.Begin();
        
        double wnt_gradient_at_cell0 = p_wnt_gradient->GetWntLevel(&(*cell_iter));
        
        while(cell_iter!=crypt.End())
        {
            TS_ASSERT_DELTA(p_wnt_gradient->GetWntLevel(&(*cell_iter)), wnt_gradient_at_cell0, 1e-12);
            ++cell_iter;
        }
                
        SimulationTime::Destroy();
        WntGradient::Destroy();
    }
    
    void TestArchiveWntGradient()
    {
        CancerParameters::Instance()->Reset();

        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "wnt_grad.arch";
        
        // Create an ouput archive
        {
            
            WntGradient::Instance()->SetType(LINEAR);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const WntGradient&>(*WntGradient::Instance());
            
            WntGradient::Destroy();
        }
        
        {
            WntGradient* p_wnt = WntGradient::Instance();
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // Restore from the archive
            input_arch >> *p_wnt;
           
            double height = 21.0;
            double wnt_level = p_wnt->GetWntLevel(height);
            
            TS_ASSERT_DELTA(wnt_level, 1.0-height/CancerParameters::Instance()->GetCryptLength(), 1e-9);
            WntGradient::Destroy();
        }
    }
    
    
    void TestSingletonnessOfWntGradient()
    {
        CancerParameters::Instance()->Reset();
        CancerParameters *params = CancerParameters::Instance();
        
        WntGradient* p_wnt_gradient = WntGradient::Instance();
        p_wnt_gradient->SetType(NONE);
        
        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        TS_ASSERT_THROWS_ANYTHING(p_wnt_gradient->SetType(NONE));
        
        WntGradient::Destroy();   
 
        p_wnt_gradient = WntGradient::Instance();
        p_wnt_gradient->SetType(LINEAR);
        
        height = 100;
        wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 1.0-height/params->GetCryptLength(), 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        
        TS_ASSERT_THROWS_ANYTHING(p_wnt_gradient->SetConstantWntValueForTesting(-10));
        
        WntGradient::Destroy();
    }
    
    
    void TestWntInitialisationSetup() throw(Exception)
    {
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
                
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        std::vector<WntCellCycleModel*> models;
        
        std::vector<TissueCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_model);
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);

            cells.push_back(cell);
            models.push_back(p_model);
        }
        
        // create the crypt
        Tissue<2> crypt(mesh,cells);
        
        CancerParameters::Instance()->SetCryptLength(1.0);

        WntGradient::Instance()->SetType(LINEAR);
        WntGradient::Instance()->SetTissue(crypt);
        
        Tissue<2>::Iterator iter = crypt.Begin();
        
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

        WntGradient::Destroy();
    }
};

#endif /*TESTWNTGRADIENT_HPP_*/
