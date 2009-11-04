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



#include "CardiacSimulation.hpp"


 
 
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
    switch(HeartConfig::Instance()->GetDomain())
    {
        case cp::domain_type::Mono :
        {
            switch (HeartConfig::Instance()->GetSpaceDimension())
            {
                case 3:
                {
                    MonodomainProblem<3>* p_mono_problem;
                    
                    if (HeartConfig::Instance()->IsSimulationDefined())
                    {
                        HeartConfigRelatedCellFactory<3> cell_factory;
                        p_mono_problem = new MonodomainProblem<3>(&cell_factory);

                        p_mono_problem->Initialise();

                        p_mono_problem->ConvertOutputToMeshalyzerFormat(true);
                    }
                    else // (HeartConfig::Instance()->IsSimulationResumed())
                    {
                        p_mono_problem = CardiacSimulationArchiver<MonodomainProblem<3> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                    }

                    p_mono_problem->Solve();

                    if (HeartConfig::Instance()->GetSaveSimulation())
                    {
                        std::stringstream directory;
                        directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                        CardiacSimulationArchiver<MonodomainProblem<3> >::Save(*p_mono_problem, directory.str(), false);
                    }

                    delete p_mono_problem;

                    break;
                }

                case 2:
                {
                    MonodomainProblem<2>* p_mono_problem;
                    
                    if (HeartConfig::Instance()->IsSimulationDefined())
                    {
                        HeartConfigRelatedCellFactory<2> cell_factory;
                        p_mono_problem = new MonodomainProblem<2>(&cell_factory);

                        p_mono_problem->Initialise();

                        p_mono_problem->ConvertOutputToMeshalyzerFormat(true);
                    }
                    else // (HeartConfig::Instance()->IsSimulationResumed())
                    {
                        p_mono_problem = CardiacSimulationArchiver<MonodomainProblem<2> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                    }

                    p_mono_problem->Solve();

                    if (HeartConfig::Instance()->GetSaveSimulation())
                    {
                        std::stringstream directory;
                        directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                        CardiacSimulationArchiver<MonodomainProblem<2> >::Save(*p_mono_problem, directory.str(), false);
                    }

                    delete p_mono_problem;

                    break;
                }

                case 1:
                {
                    MonodomainProblem<1>* p_mono_problem;
                    
                    if (HeartConfig::Instance()->IsSimulationDefined())
                    {
                        HeartConfigRelatedCellFactory<1> cell_factory;
                        p_mono_problem = new MonodomainProblem<1>(&cell_factory);

                        p_mono_problem->Initialise();

                        p_mono_problem->ConvertOutputToMeshalyzerFormat(true);
                    }
                    else // (HeartConfig::Instance()->IsSimulationResumed())
                    {
                        p_mono_problem = CardiacSimulationArchiver<MonodomainProblem<1> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                    }

                    p_mono_problem->Solve();

                    if (HeartConfig::Instance()->GetSaveSimulation())
                    {
                        std::stringstream directory;
                        directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                        CardiacSimulationArchiver<MonodomainProblem<1> >::Save(*p_mono_problem, directory.str(), false);
                    }

                    delete p_mono_problem;

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
                    BidomainProblem<3>* p_bi_problem;
                    
                    if (HeartConfig::Instance()->IsSimulationDefined())
                    {
                        HeartConfigRelatedCellFactory<3> cell_factory;
                        p_bi_problem = new BidomainProblem<3>(&cell_factory);

                        p_bi_problem->Initialise();

                        p_bi_problem->ConvertOutputToMeshalyzerFormat(true);
                    }
                    else // (HeartConfig::Instance()->IsSimulationResumed())
                    {
                        p_bi_problem = CardiacSimulationArchiver<BidomainProblem<3> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                    }

                    p_bi_problem->Solve();

                    if (HeartConfig::Instance()->GetSaveSimulation())
                    {
                        std::stringstream directory;
                        directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                        CardiacSimulationArchiver<BidomainProblem<3> >::Save(*p_bi_problem, directory.str(), false);
                    }

                    delete p_bi_problem;

                    break;
                }
                case 2:
                {
                    BidomainProblem<2>* p_bi_problem;
                    
                    if (HeartConfig::Instance()->IsSimulationDefined())
                    {
                        HeartConfigRelatedCellFactory<2> cell_factory;
                        p_bi_problem = new BidomainProblem<2>(&cell_factory);

                        p_bi_problem->Initialise();

                        p_bi_problem->ConvertOutputToMeshalyzerFormat(true);
                    }
                    else // (HeartConfig::Instance()->IsSimulationResumed())
                    {
                        p_bi_problem = CardiacSimulationArchiver<BidomainProblem<2> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                    }

                    p_bi_problem->Solve();

                    if (HeartConfig::Instance()->GetSaveSimulation())
                    {
                        std::stringstream directory;
                        directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                        CardiacSimulationArchiver<BidomainProblem<2> >::Save(*p_bi_problem, directory.str(), false);
                    }

                    delete p_bi_problem;

                    break;
                }
                case 1:
                {
                    BidomainProblem<1>* p_bi_problem;
                    
                    if (HeartConfig::Instance()->IsSimulationDefined())
                    {
                        HeartConfigRelatedCellFactory<1> cell_factory;
                        p_bi_problem = new BidomainProblem<1>(&cell_factory);

                        p_bi_problem->Initialise();

                        p_bi_problem->ConvertOutputToMeshalyzerFormat(true);
                    }
                    else // (HeartConfig::Instance()->IsSimulationResumed())
                    {
                        p_bi_problem = CardiacSimulationArchiver<BidomainProblem<1> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                    }

                    p_bi_problem->Solve();

                    if (HeartConfig::Instance()->GetSaveSimulation())
                    {
                        std::stringstream directory;
                        directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                        CardiacSimulationArchiver<BidomainProblem<1> >::Save(*p_bi_problem, directory.str(), false);
                    }

                    delete p_bi_problem;
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
            //If the domain is not set correctly then the XML parser will have picked it up before now!
            NEVER_REACHED;
        }
    }

}

