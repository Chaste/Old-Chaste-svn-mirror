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


#ifndef HEARTCONFIG_HPP_
#define HEARTCONFIG_HPP_

#include "ChasteParameters.hpp"
#include <iostream>
#include "UblasCustomFunctions.hpp"
#include "Exception.hpp"
#include <vector>
#include "AbstractStimulusFunction.hpp"
#include "SimpleStimulus.hpp"
#include "ChasteCuboid.hpp"

class HeartConfig
{
public:
    /**
     * Call this method to access the global parameters holder.
     * 
     * @return a single instance of the class
     */
    static HeartConfig* Instance();
    
    chaste_parameters_type* UserParameters();
    chaste_parameters_type* DefaultParameters();
    
    void SetDefaultsFile(std::string fileName);
    void SetParametersFile(std::string fileName);
    static void Destroy();
    
    // Simulation
    double GetSimulationDuration();
    domain_type GetDomain();
    ionic_model_type GetIonicModel();
    void GetStimuli(std::vector<SimpleStimulus>& stimuliApplied, std::vector<ChasteCuboid>& stimulatedAreas);
    void GetCellHeterogeneities(std::vector<ChasteCuboid>& cellHeterogeneityAreas,
    							std::vector<double>& scaleFactorGks,
    							std::vector<double>& scaleFactorIto);
    void GetConductivityHeterogeneities(std::vector<ChasteCuboid>& conductivitiesHeterogeneityAreas,
				  					 	std::vector< c_vector<double,3> >& intraConductivities,
										std::vector< c_vector<double,3> >& extraConductivities);
    std::string GetOutputDirectory();
    
    // Physiological
    c_vector<double, 3> GetIntracellularConductivities();
    c_vector<double, 3> GetExtracellularConductivities();
    double GetSurfaceAreaToVolumeRatio();
    double GetCapacitance();
    
    // Numerical
    double GetOdeTimestep();
    double GetPdeTimestep();
    double GetPrintingTimestep();
    
    bool GetUseAbsoluteTolerance();
    double GetAbsoluteTolerance();

    bool GetUseRelativeTolerance();
    double GetRelativeTolerance();

    ksp_solver_type GetKSPSolver();
    ksp_preconditioner_type GetKSPPreconditioner();
    
private:
    HeartConfig();
    ~HeartConfig();
    
    
    /** The single instance of the class */
    static HeartConfig* mpInstance;
    
    chaste_parameters_type* mpUserParameters;
    chaste_parameters_type* mpDefaultParameters;
    
    // Misc
    template<class TYPE> 
    TYPE* DecideLocation(TYPE* ptr1, TYPE* ptr2, std::string nameParameter);
    //Utility method to parse an XML parameters file
    chaste_parameters_type* ReadFile(std::string fileName);
  
};

#endif /*HEARTCONFIG_HPP_*/
