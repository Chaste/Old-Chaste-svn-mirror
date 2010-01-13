/*

Copyright (C) University of Oxford, 2005-2010

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


#include "CardiacSimulation.hpp"


CardiacSimulation::CardiacSimulation(std::string parameterFileName)
{
    // If we have been passed an XML file then parse the XML file, otherwise throw
    if (parameterFileName == "")
    {
        EXCEPTION("No XML file name given");
    }
    ReadParametersFromFile(parameterFileName);
    Run();
    HeartEventHandler::Headings();
    HeartEventHandler::Report();
}


void CardiacSimulation::ReadParametersFromFile(std::string parameterFileName)
{
    try
    {
        // Try the hardcoded schema location first
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile(parameterFileName);
    }
    catch (Exception& e)
    {
        if (e.CheckShortMessageContains("Missing file parsing configuration") == "")
        {
            // Try using the schema location given in the XML
            HeartConfig::Instance()->Reset();
            HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
            HeartConfig::Instance()->SetParametersFile(parameterFileName);
        }
        else
        {
            throw e;
        }
    }
}


void CardiacSimulation::Run()
{
    switch (HeartConfig::Instance()->GetDomain())
    {
        case cp::domain_type::Mono :
        {
            switch (HeartConfig::Instance()->GetSpaceDimension())
            {
                case 3:
                {
                    CreateAndRun<MonodomainProblem<3>,3>();
                    break;
                }

                case 2:
                {
                    CreateAndRun<MonodomainProblem<2>,2>();
                    break;
                }

                case 1:
                {
                    CreateAndRun<MonodomainProblem<1>,1>();
                    break;
                }
                default :
                    EXCEPTION("Monodomain space dimension not supported: should be 1, 2 or 3");
            }
            break;
        }

        case cp::domain_type::Bi :
        {
            switch (HeartConfig::Instance()->GetSpaceDimension())
            {
                case 3:
                {
                    CreateAndRun<BidomainProblem<3>,3>();
                    break;
                }
                case 2:
                {
                    CreateAndRun<BidomainProblem<2>,2>();
                    break;
                }
                case 1:
                {
                    CreateAndRun<BidomainProblem<1>,1>();
                    break;
                }
                default :
                {
                    EXCEPTION("Bidomain space dimension not supported: should be 1, 2 or 3");
                }
            }
            break;
        }
        default :
        {
            // If the domain is not set correctly then the XML parser will have picked it up before now!
            NEVER_REACHED;
        }
    }
}

