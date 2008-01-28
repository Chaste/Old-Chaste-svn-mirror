#ifndef TESTSPIRALPARAMETERSREADER_HPP_
#define TESTSPIRALPARAMETERSREADER_HPP_

#include <cxxtest/TestSuite.h>
#include <memory>
#include "ChasteParameters.hpp"


using std::auto_ptr;

class TestChasteParametersReader : public CxxTest::TestSuite
{
public:
    void TestRead()
    {
        try
        {
            auto_ptr<chaste_parameters_type> p (ChasteParameters("heart/test/data/ChasteParameters.xml"));
            
            TS_ASSERT_EQUALS(p->SimulationDuration(), 10.0);
            TS_ASSERT_EQUALS(p->SimulationDuration(), 10.0);
            TS_ASSERT_EQUALS(p->SlabWidth(), 2.0);
            TS_ASSERT_EQUALS(p->SlabHeight(), 1.0);
            TS_ASSERT_EQUALS(p->InterNodeSpace(), 0.25);
            TS_ASSERT_EQUALS(p->OutputDirectory(), "ChasteResults");
            TS_ASSERT_EQUALS(p->MeshOutputDirectory(), "Slab");
            TS_ASSERT_EQUALS(p->Domain(), domain_type::Mono);
            TS_ASSERT_EQUALS(p->IonicModel(), ionic_model_type::FaberRudy2000Version3);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }

    // For coverage purposes
    void TestReadConst()
    {
        try
        {
            auto_ptr<const chaste_parameters_type> p (ChasteParameters("heart/test/data/ChasteParameters.xml"));
            
            TS_ASSERT_EQUALS(p->SimulationDuration(), 10.0);
            TS_ASSERT_EQUALS(p->SimulationDuration(), 10.0);
            TS_ASSERT_EQUALS(p->SlabWidth(), 2.0);
            TS_ASSERT_EQUALS(p->SlabHeight(), 1.0);
            TS_ASSERT_EQUALS(p->InterNodeSpace(), 0.25);
            TS_ASSERT_EQUALS(p->OutputDirectory(), "ChasteResults");
            TS_ASSERT_EQUALS(p->MeshOutputDirectory(), "Slab");
            TS_ASSERT_EQUALS(p->Domain(), domain_type::Mono);
            TS_ASSERT_EQUALS(p->IonicModel(), ionic_model_type::FaberRudy2000Version3);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }
};
#endif /*TESTSPIRALPARAMETERSREADER_HPP_*/
