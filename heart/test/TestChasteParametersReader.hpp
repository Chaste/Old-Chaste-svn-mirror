/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTSPIRALPARAMETERSREADER_HPP_
#define TESTSPIRALPARAMETERSREADER_HPP_

#include <cxxtest/TestSuite.h>
#include <memory>
#include "ChasteParameters.hpp"


using std::auto_ptr;

class TestChasteParametersReader : public CxxTest::TestSuite
{
public:
    void TestReadWithSlab()
    {
        try
        {
            auto_ptr<chaste_parameters_type> params (ChasteParameters( "heart/test/data/ChasteParametersFullFormat.xml"));

            simulation_type simulation_params = params->Simulation();

            TS_ASSERT_EQUALS(simulation_params.SimulationDuration().get(), 10.0);

            TS_ASSERT(simulation_params.Mesh().present());
            TS_ASSERT(simulation_params.Mesh().get().LoadMesh() == NULL);
            TS_ASSERT(simulation_params.Mesh().get().Slab() != NULL);

            TS_ASSERT_EQUALS(simulation_params.Mesh().get().Slab()->SlabX(), 4.0);
            TS_ASSERT_EQUALS(simulation_params.Mesh().get().Slab()->SlabY(), 0.1);
            TS_ASSERT_EQUALS(simulation_params.Mesh().get().Slab()->SlabZ(), 2.0);
            TS_ASSERT_EQUALS(simulation_params.Mesh().get().Slab()->InterNodeSpace(), 0.1);

            physiological_type physiological_params = params->Physiological();

            TS_ASSERT(physiological_params.IntracellularConductivities().present());
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().get().longi(), 1.75);
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().get().trans(), 1.75);
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().get().normal(), 1.75);

            TS_ASSERT(physiological_params.ExtracellularConductivities().present());
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().get().longi(), 7.0);
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().get().trans(), 7.0);
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().get().normal(), 7.0);

            TS_ASSERT_EQUALS(simulation_params.OutputDirectory().get(), "ChasteResults");
            TS_ASSERT_EQUALS(simulation_params.Domain().get(), domain_type::Mono);
            TS_ASSERT_EQUALS(simulation_params.IonicModel().get(), ionic_model_type::FaberRudy2000);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }

    void TestReadWithLoadMesh()
    {
        try
        {
            auto_ptr<chaste_parameters_type> params (ChasteParameters("heart/test/data/ChasteParametersLoadMesh.xml"));

            simulation_type simulation_params = params->Simulation();

            TS_ASSERT_EQUALS(simulation_params.SimulationDuration().get(), 10.0);

            TS_ASSERT(simulation_params.Mesh().present());
            TS_ASSERT(simulation_params.Mesh().get().LoadMesh() != NULL);
            TS_ASSERT(simulation_params.Mesh().get().Slab() == NULL);

            TS_ASSERT_EQUALS(simulation_params.Mesh().get().LoadMesh()->name(), "foo");

            physiological_type physiological_params = params->Physiological();

            TS_ASSERT(physiological_params.IntracellularConductivities().present());
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().get().longi(), 1.75);
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().get().trans(), 1.75);
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().get().normal(), 1.75);

            TS_ASSERT(physiological_params.ExtracellularConductivities().present());
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().get().longi(), 7.0);
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().get().trans(), 7.0);
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().get().normal(), 7.0);

            TS_ASSERT_EQUALS(simulation_params.OutputDirectory().get(), "ChasteResults");
            TS_ASSERT_EQUALS(simulation_params.Domain().get(), domain_type::Mono);
            TS_ASSERT_EQUALS(simulation_params.IonicModel().get(), ionic_model_type::FaberRudy2000);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }

    void TestUpdate()
    {
        try
        {
            auto_ptr<chaste_parameters_type> params (ChasteParameters("heart/test/data/ChasteParametersLoadMesh.xml"));
            physiological_type physiological_params = params->Physiological();

            physiological_params.ExtracellularConductivities().get().longi() = 9.0;
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().get().longi(), 9.0);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }

};
#endif /*TESTSPIRALPARAMETERSREADER_HPP_*/
