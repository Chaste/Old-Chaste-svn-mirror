#ifndef TESTWNTGRADIENT_HPP_
#define TESTWNTGRADIENT_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"
#include "SingletonWntGradient.hpp"
#include "WntGradientTypes.hpp"
#include "Crypt.cpp"
#include "MeinekeCryptCell.hpp"
#include "WntCellCycleModel.hpp"

class TestWntGradient : public CxxTest::TestSuite
{
public:    
    void TestNoWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        SingletonWntGradient* p_wnt_gradient = SingletonWntGradient::Instance();
        p_wnt_gradient->SetType(NONE);
        
        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        SingletonWntGradient::Destroy();
    }
    
    void TestLinearWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        SingletonWntGradient* p_wnt_gradient = SingletonWntGradient::Instance();
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
        
        SingletonWntGradient::Destroy();
    }
    
    
    void TestOffsetLinearWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        SingletonWntGradient* p_wnt_gradient = SingletonWntGradient::Instance();
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
        
        SingletonWntGradient::Destroy();
    }
//    
//    void TestArchiveWntGradient()
//    {
//        CancerParameters::Instance()->Reset();
//
//        OutputFileHandler handler("archive",false);
//        std::string archive_filename;
//        archive_filename = handler.GetTestOutputDirectory() + "wnt_grad.arch";
//        
//        // Create an ouput archive
//        {
//            WntGradientType this_type = LINEAR;
//            
//            WntGradient* const p_wnt_gradient = new WntGradient(this_type);
//            
//            std::ofstream ofs(archive_filename.c_str());
//            boost::archive::text_oarchive output_arch(ofs);
//            
//            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
//            output_arch << p_wnt_gradient;
//            
//            CancerParameters *inst1 = CancerParameters::Instance();
//            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
//            
//            delete p_wnt_gradient;
//        }
//        
//        {
//            WntGradient* p_wnt;
//            
//            CancerParameters *inst1 = CancerParameters::Instance();
//            
//            inst1->SetSG2MDuration(101.0);
//            
//            // Create an input archive
//            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
//            boost::archive::text_iarchive input_arch(ifs);
//            
//            // restore from the archive
//            input_arch >> *inst1;
//            input_arch >> p_wnt;
//            
//            CancerParameters *inst2 = CancerParameters::Instance();
//            TS_ASSERT_EQUALS(inst1, inst2);
//            
//            // Check
//            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
//            double height = 21.0;
//            double wnt_level = p_wnt->GetWntLevel(height);
//            
//            TS_ASSERT_DELTA(wnt_level, 1.0-height/inst1->GetCryptLength(), 1e-9);
//        }
//    }
//    
//    
    void TestSingletonWntGradient()
    {
        CancerParameters::Instance()->Reset();
        CancerParameters *params = CancerParameters::Instance();
        
        SingletonWntGradient* p_wnt_gradient = SingletonWntGradient::Instance();
        p_wnt_gradient->SetType(NONE);
        
        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        TS_ASSERT_THROWS_ANYTHING(p_wnt_gradient->SetType(NONE));
        
        SingletonWntGradient::Destroy();   
 
        p_wnt_gradient = SingletonWntGradient::Instance();
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
        
        SingletonWntGradient::Destroy();
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
        
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
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

#endif /*TESTWNTGRADIENT_HPP_*/
