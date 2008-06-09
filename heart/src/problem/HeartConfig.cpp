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

using namespace xsd::cxx::tree;

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
{
    assert(mpInstance == NULL);
    
    mpDefaultParameters = ReadFile("ChasteDefaults.xml");
}

HeartConfig::~HeartConfig()
{
    delete mpUserParameters;
    delete mpDefaultParameters;  	
}

void HeartConfig::SetDefaultsFile(std::string fileName)
{
    mpDefaultParameters = ReadFile(fileName);
}
std::string GetOutputDirectory();
chaste_parameters_type* HeartConfig::ReadFile(std::string fileName)
{
	// get the parameters using the method 'ChasteParameters(filename)',
    // which returns a std::auto_ptr. We don't want to use a std::auto_ptr because
    // it will delete memory when out of scope, or no longer point when it is copied,
    // so we reallocate memory using a normal pointer and copy the data to thereS
 	try
    {
        std::auto_ptr<chaste_parameters_type> p_default(ChasteParameters(fileName));
        return new chaste_parameters_type(*p_default);
    }
    catch (const xml_schema::exception& e)
    {
         std::cerr << e << std::endl;
         EXCEPTION("XML parsing error in configuration file: " + fileName);
    }
}

void HeartConfig::SetParametersFile(std::string fileName)
{
    mpUserParameters = ReadFile(fileName);
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
	delete mpInstance;
	mpInstance = NULL;  
}

template<class TYPE>
TYPE* HeartConfig::DecideLocation(TYPE* ptr1, TYPE* ptr2, std::string nameParameter)
{
    if (ptr1->present())
    {
        return ptr1;
    }
    else
    {
        if (ptr2->present())
        {
            return ptr2;            
        }
        else
        {
            EXCEPTION("No " + nameParameter + " provided (neither default nor user defined)");
        }             
    }

}

double HeartConfig::GetSimulationDuration()
{
	return DecideLocation( & mpUserParameters->Simulation().SimulationDuration(), 
						   & mpDefaultParameters->Simulation().SimulationDuration(), 
						   "SimulationDuration")->get(); 
}

domain_type HeartConfig::GetDomain()
{
	return DecideLocation( & mpUserParameters->Simulation().Domain(), 
						   & mpDefaultParameters->Simulation().Domain(), 
						   "Domain")->get();
}

ionic_model_type HeartConfig::GetIonicModel()
{
	return DecideLocation( & mpUserParameters->Simulation().IonicModel(), 
	                       & mpDefaultParameters->Simulation().IonicModel(), 
	                       "IonicModel")->get(); 
}

void HeartConfig::GetStimuli(std::vector<SimpleStimulus>& stimuliApplied, std::vector<ChasteCuboid>& stimulatedAreas)
{

    simulation_type::Stimuli::_xsd_Stimuli_::Stimuli::Stimulus::container&
         stimuli = DecideLocation( & mpUserParameters->Simulation().Stimuli(), 
	                       & mpDefaultParameters->Simulation().Stimuli(), 
	                       "Stimuli")->get().Stimulus();
    for (simulation_type::Stimuli::_xsd_Stimuli_::Stimuli::Stimulus::iterator i = stimuli.begin();
         i != stimuli.end();
         ++i)
    {                     
        stimulus_type stimulus(*i);           
        point_type point_a = stimulus.Location().CornerA();
        point_type point_b = stimulus.Location().CornerB();
        
        ChastePoint<3> chaste_point_a ( point_a.x(), 
                                        point_a.y(),
                                        point_a.z());

        ChastePoint<3> chaste_point_b ( point_b.x(),
                                        point_b.y(),
                                        point_b.z());
                    
        stimuliApplied.push_back( SimpleStimulus(stimulus.Strength(), stimulus.Duration(), stimulus.Delay() ) );
        stimulatedAreas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );
    }	
}

void HeartConfig::GetCellHeterogeneities(std::vector<ChasteCuboid>& cellHeterogeneityAreas,
    								     std::vector<double>& scaleFactorGks,
    							         std::vector<double>& scaleFactorIto)
{
	simulation_type::CellHeterogeneities::_xsd_CellHeterogeneities_::CellHeterogeneities::CellHeterogeneity::container&
         cell_heterogeneity = DecideLocation( & mpUserParameters->Simulation().CellHeterogeneities(), 
	                       					  & mpDefaultParameters->Simulation().CellHeterogeneities(), 
	                       					  "CellHeterogeneities")->get().CellHeterogeneity();
	
    for (simulation_type::CellHeterogeneities::_xsd_CellHeterogeneities_::CellHeterogeneities::CellHeterogeneity::iterator i = cell_heterogeneity.begin();
         i != cell_heterogeneity.end();
         ++i)
    {                     
        cell_heterogeneity_type ht(*i);           
        point_type point_a = ht.Location().CornerA();
        point_type point_b = ht.Location().CornerB();
        
        ChastePoint<3> chaste_point_a (point_a.x(), 
                                       point_a.y(),
                                       point_a.z());

        ChastePoint<3> chaste_point_b (point_b.x(),
                                       point_b.y(),
                                       point_b.z());
                    
        scaleFactorGks.push_back (ht.ScaleFactorGks());
        scaleFactorIto.push_back (ht.ScaleFactorIto());                                    
        cellHeterogeneityAreas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );
    }        
}

void HeartConfig::GetConductivityHeterogeneities(std::vector<ChasteCuboid>& conductivitiesHeterogeneityAreas,
					  			 	std::vector< c_vector<double,3> >& intraConductivities,
									std::vector< c_vector<double,3> >& extraConductivities)
{
	simulation_type::ConductivityHeterogeneities::_xsd_ConductivityHeterogeneities_::ConductivityHeterogeneities::ConductivityHeterogeneity::container&
         conductivity_heterogeneity = DecideLocation( & mpUserParameters->Simulation().ConductivityHeterogeneities(), 
	                       							  & mpDefaultParameters->Simulation().ConductivityHeterogeneities(), 
	                       					  		  "CellHeterogeneities")->get().ConductivityHeterogeneity();
	
    for (simulation_type::ConductivityHeterogeneities::_xsd_ConductivityHeterogeneities_::ConductivityHeterogeneities::ConductivityHeterogeneity::iterator i = conductivity_heterogeneity.begin();
         i != conductivity_heterogeneity.end();
         ++i)
    {                     
        conductivity_heterogeneity_type ht(*i);           
        point_type point_a = ht.Location().CornerA();
        point_type point_b = ht.Location().CornerB();
        
        ChastePoint<3> chaste_point_a (point_a.x(), 
                                       point_a.y(),
                                       point_a.z());

        ChastePoint<3> chaste_point_b (point_b.x(),
                                       point_b.y(),
                                       point_b.z());
                    
        conductivitiesHeterogeneityAreas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );
        
        if (ht.IntracellularConductivities().present())
        {
        	double intra_x = ht.IntracellularConductivities().get().longi();
        	double intra_y = ht.IntracellularConductivities().get().trans();
        	double intra_z = ht.IntracellularConductivities().get().normal();
        	
        	intraConductivities.push_back( Create_c_vector(intra_x, intra_y, intra_z) );
        }
        else
        {
        	intraConductivities.push_back( GetIntracellularConductivities() );
        }

        if (ht.ExtracellularConductivities().present())
        {
        	double extra_x = ht.ExtracellularConductivities().get().longi();
        	double extra_y = ht.ExtracellularConductivities().get().trans();
        	double extra_z = ht.ExtracellularConductivities().get().normal();
        	
        	extraConductivities.push_back( Create_c_vector(extra_x, extra_y, extra_z) );
        }
        else
        {
        	extraConductivities.push_back( GetExtracellularConductivities() );
        }

    }        
}

std::string HeartConfig::GetOutputDirectory()
{
    return DecideLocation( & mpUserParameters->Simulation().OutputDirectory(), 
                           & mpDefaultParameters->Simulation().OutputDirectory(), 
                           "OutputDirectory")->get();        
}

c_vector<double, 3> HeartConfig::GetIntracellularConductivities()
{
    optional<conductivities_type, false>* intra_conductivities  = DecideLocation( & mpUserParameters->Physiological().IntracellularConductivities(), 
                                                                                  & mpDefaultParameters->Physiological().IntracellularConductivities(), 
                                                                                  "IntracellularConductivities");           
    double intra_x_cond = intra_conductivities->get().longi();
    double intra_y_cond = intra_conductivities->get().trans();
    double intra_z_cond = intra_conductivities->get().normal();;

    return Create_c_vector(intra_x_cond, intra_y_cond, intra_z_cond);       
}

c_vector<double, 3> HeartConfig::GetExtracellularConductivities()
{
    optional<conductivities_type, false>* extra_conductivities  = DecideLocation( & mpUserParameters->Physiological().ExtracellularConductivities(), 
                                                                                  & mpDefaultParameters->Physiological().ExtracellularConductivities(), 
                                                                                  "ExtracellularConductivities");
    double extra_x_cond = extra_conductivities->get().longi();
    double extra_y_cond = extra_conductivities->get().trans();
    double extra_z_cond = extra_conductivities->get().normal();;

    return Create_c_vector(extra_x_cond, extra_y_cond, extra_z_cond);   
}

double HeartConfig::GetSurfaceAreaToVolumeRatio()
{
    return DecideLocation( & mpUserParameters->Physiological().SurfaceAreaToVolumeRatio(), 
                           & mpDefaultParameters->Physiological().SurfaceAreaToVolumeRatio(), 
                           "SurfaceAreaToVolumeRatio")->get();            
}

double HeartConfig::GetCapacitance()
{    
    return DecideLocation( & mpUserParameters->Physiological().Capacitance(), 
                           & mpDefaultParameters->Physiological().Capacitance(), 
                           "Capacitance")->get();            
}

double HeartConfig::GetOdeTimestep()
{
    return DecideLocation( & mpUserParameters->Numerical().Timesteps(), 
                           & mpDefaultParameters->Numerical().Timesteps(), 
                           "Timesteps")->get().ode();                
}

double HeartConfig::GetPdeTimestep()
{
    return DecideLocation( & mpUserParameters->Numerical().Timesteps(), 
                           & mpDefaultParameters->Numerical().Timesteps(), 
                           "Timesteps")->get().pde();                
}

double HeartConfig::GetPrintingTimestep()
{
    return DecideLocation( & mpUserParameters->Numerical().Timesteps(), 
                           & mpDefaultParameters->Numerical().Timesteps(), 
                           "Timesteps")->get().printing();                
}

bool HeartConfig::GetUseAbsoluteTolerance()
{
    ksp_use_type use_value = DecideLocation( & mpUserParameters->Numerical().KSPTolerances(), 
                                             & mpDefaultParameters->Numerical().KSPTolerances(), 
                                             "KSPTolerances")->get().use();
                                   
    return (use_value == ksp_use_type::absolute);                                   
}

double HeartConfig::GetAbsoluteTolerance()
{
    return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(), 
                           & mpDefaultParameters->Numerical().KSPTolerances(), 
                           "KSPTolerances")->get().KSPAbsolute();        
}

bool HeartConfig::GetUseRelativeTolerance()
{
    ksp_use_type use_value = DecideLocation( & mpUserParameters->Numerical().KSPTolerances(), 
                                             & mpDefaultParameters->Numerical().KSPTolerances(), 
                                             "KSPTolerances")->get().use();
                                   
    return (use_value == ksp_use_type::relative);                                   
}

double HeartConfig::GetRelativeTolerance()
{
    return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(), 
                           & mpDefaultParameters->Numerical().KSPTolerances(), 
                           "KSPTolerances")->get().KSPRelative();        
}

ksp_solver_type HeartConfig::GetKSPSolver()
{
    return DecideLocation( & mpUserParameters->Numerical().KSPSolver(), 
                           & mpDefaultParameters->Numerical().KSPSolver(), 
                           "KSPSolver")->get();    
}

ksp_preconditioner_type HeartConfig::GetKSPPreconditioner()
{
    return DecideLocation( & mpUserParameters->Numerical().KSPPreconditioner(), 
                           & mpDefaultParameters->Numerical().KSPPreconditioner(), 
                           "KSPPreconditioner")->get();
}
