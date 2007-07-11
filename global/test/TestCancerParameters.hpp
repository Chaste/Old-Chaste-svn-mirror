#ifndef TESTCANCERPARAMETERS_HPP_
#define TESTCANCERPARAMETERS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"



class TestCancerParameters : public CxxTest::TestSuite
{
private:
    void CheckValuesAreTheDefaultValues()
    {
        CancerParameters *inst = CancerParameters::Instance();
        
        TS_ASSERT_DELTA(inst->GetSG2MDuration(), 10.0 , 1e-12);
        TS_ASSERT_DELTA(inst->GetStemCellCycleTime(), 24.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetTransitCellCycleTime(), 12.0, 1e-12);
        TS_ASSERT_EQUALS(inst->GetMaxTransitGenerations(), 3u);
        TS_ASSERT_DELTA(inst->GetCryptLength(), 22.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetCryptWidth(), 10.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetSpringStiffness(), 15.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetDampingConstantNormal(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetDampingConstantMutant(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetApoptosisTime(), 0.25, 1e-12);
    }

public:

    void TestConstructor()
    {
        CheckValuesAreTheDefaultValues();
    }
        
    void TestReset()
    {
        CancerParameters* inst = CancerParameters::Instance();
        
        inst->SetSG2MDuration(11.0);
        inst->SetStemCellCycleTime(35.0);
        inst->SetTransitCellCycleTime(45.0);
        inst->SetMaxTransitGenerations(666u);
        inst->SetCryptLength(100.0);
        inst->SetSpringStiffness(20.0);
        inst->SetDampingConstantNormal(2.0);
        inst->SetDampingConstantMutant(3.0);
        inst->SetApoptosisTime(0.3);
        
        inst->Reset();

        CheckValuesAreTheDefaultValues();
    }
    

    void TestGettersAndSetters()
    {
        CancerParameters *inst1 = CancerParameters::Instance();
        
        inst1->SetSG2MDuration(11.0);
        inst1->SetStemCellCycleTime(35.0);
        inst1->SetTransitCellCycleTime(45.0);
        inst1->SetMaxTransitGenerations(666u);
        inst1->SetCryptLength(100.0);
        inst1->SetSpringStiffness(20.0);
        inst1->SetDampingConstantNormal(2.0);
        inst1->SetDampingConstantMutant(3.0);
        inst1->SetApoptosisTime(0.3);
        
        CancerParameters *inst2 = CancerParameters::Instance();
        
        TS_ASSERT_DELTA(inst2->GetSG2MDuration(), 11.0 , 1e-12);
        TS_ASSERT_DELTA(inst2->GetStemCellCycleTime(), 35.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetTransitCellCycleTime(), 45.0, 1e-12);
        TS_ASSERT_EQUALS(inst2->GetMaxTransitGenerations(), 666u);
        TS_ASSERT_DELTA(inst2->GetCryptLength(), 100.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetSpringStiffness(), 20.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetDampingConstantNormal(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetDampingConstantMutant(), 3.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetApoptosisTime(), 0.3, 1e-12);
    }
    
    void TestArchiveCancerParameters()
    {
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "cancer_params.arch";
        
        // Create an ouput archive
        {
            CancerParameters *inst1 = CancerParameters::Instance();
            // Mess up the cancer parameters
            inst1->SetSG2MDuration(11.0);
            inst1->SetStemCellCycleTime(35.0);
            inst1->SetTransitCellCycleTime(45.0);
            inst1->SetMaxTransitGenerations(666u);
            inst1->SetCryptLength(100.0);
            inst1->SetSpringStiffness(20.0);
            inst1->SetDampingConstantNormal(2.0);
            inst1->SetDampingConstantMutant(3.0);
            inst1->SetApoptosisTime(0.3);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            // save messed up parameters
            output_arch << static_cast<const CancerParameters&>(*inst1);
            
        }
        
        {
            CancerParameters *inst1 = CancerParameters::Instance();
            // restore to nice parameters
            inst1->SetSG2MDuration(10.0);
            inst1->SetStemCellCycleTime(24.0);
            inst1->SetTransitCellCycleTime(12.0);
            inst1->SetMaxTransitGenerations(3u);
            inst1->SetCryptLength(22.0);
            inst1->SetApoptosisTime(0.25);
            inst1->SetSpringStiffness(30.0);
            inst1->SetDampingConstantNormal(1.0);
            inst1->SetDampingConstantMutant(2.0);

            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore messed up parameters from the archive
            input_arch >> *inst1;
            // check they are messed up.
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(), 11.0 , 1e-12);
            TS_ASSERT_DELTA(inst1->GetStemCellCycleTime(), 35.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetTransitCellCycleTime(), 45.0, 1e-12);
            TS_ASSERT_EQUALS(inst1->GetMaxTransitGenerations(), 666u);
            TS_ASSERT_DELTA(inst1->GetCryptLength(), 100.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetSpringStiffness(), 20.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetDampingConstantNormal(), 2.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetDampingConstantMutant(), 3.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetApoptosisTime(), 0.3, 1e-12);
        }
    }
    
};

#endif /*TESTCANCERPARAMETERS_HPP_*/
