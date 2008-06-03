/*

Copyright (C) University of Oxford, 2008

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


#include "HeartConfig.hpp"

HeartConfig* HeartConfig::mpInstance = NULL;

HeartConfig* HeartConfig::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new HeartConfig;
    }
    return mpInstance;
}

HeartConfig::HeartConfig()
: mDefaultsFile("ChasteDefaults.xml")
{
    assert(mpInstance == NULL);
}

void HeartConfig::SetDefaultsFile(std::string fileName)
{
    mDefaultsFile = fileName;
}

void HeartConfig::SetParametersFile(std::string fileName)
{
    // get the parameters using the method 'ChasteParameters(filename)',
    // which returns a std::auto_ptr. We don't want to use a std::auto_ptr because
    // it will delete memory when out of scope, or no longer point when it is copied,
    // so we reallocate memory using a normal pointer and copy the data to there
    try
    {   
        std::auto_ptr<chaste_parameters_type> p_user(ChasteParameters(fileName.c_str()));
        mpUserParameters = new chaste_parameters_type(*p_user);
        assert(mpUserParameters);
    }
    catch (const xml_schema::exception& e)
    {
         std::cerr << e << std::endl;
         EXCEPTION("XML parsing error in user provided configuration file");
    }
    

    // Read default values
    try
    {
        std::auto_ptr<chaste_parameters_type> p_default(ChasteParameters(mDefaultsFile));
        mpDefaultParameters = new chaste_parameters_type(*p_default);
        assert(mpDefaultParameters);
    }
    catch (const xml_schema::exception& e)
    {
         std::cerr << e << std::endl;
         EXCEPTION("XML parsing error in default configuration file");
    }

}

chaste_parameters_type* HeartConfig::UserParameters()
{
    return mpUserParameters;
}

chaste_parameters_type* HeartConfig::DefaultParameters()
{
    return mpDefaultParameters;
}

void HeartConfig::Destroy()
{
    delete mpUserParameters;
    delete mpDefaultParameters;    
}

ionic_model_type HeartConfig::GetIonicModel()
{
    if (mpUserParameters->Simulation().IonicModel().present())
    {
        return mpUserParameters->Simulation().IonicModel().get();
    }
    else
    {
        if (mpDefaultParameters->Simulation().IonicModel().present())
        {
            return mpDefaultParameters->Simulation().IonicModel().get();            
        }
        else
        {
            EXCEPTION("No IonicModel provided (neither default of user defined)");
        }             
    }
}

c_vector<double, 3> HeartConfig::GetIntracellularConductivities()
{
    double intra_x_cond;
    double intra_y_cond;
    double intra_z_cond;                    
    
    if (mpUserParameters->Simulation().IonicModel().present())
    {
        intra_x_cond = mpUserParameters->Physiological().IntracellularConductivities().get().longi();
        intra_y_cond = mpUserParameters->Physiological().IntracellularConductivities().get().trans();
        intra_z_cond = mpUserParameters->Physiological().IntracellularConductivities().get().normal();                
    }
    else
    {
        if (mpDefaultParameters->Simulation().IonicModel().present())
        {
            intra_x_cond = mpDefaultParameters->Physiological().IntracellularConductivities().get().longi();
            intra_y_cond = mpDefaultParameters->Physiological().IntracellularConductivities().get().trans();
            intra_z_cond = mpDefaultParameters->Physiological().IntracellularConductivities().get().normal();                
        }
        else
        {
            EXCEPTION("No IntracellularConductivities provided (neither default of user defined)");
        }             
    }        

    return Create_c_vector(intra_x_cond, intra_y_cond, intra_z_cond);   
}
