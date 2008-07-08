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
private:
	void CheckTimeSteps() const;
	
public:
    /**
     * Call this method to access the global parameters holder.
     *
     * @return a single instance of the class
     */
    static HeartConfig* Instance();

    void SetDefaultsFile(std::string fileName);
    void SetParametersFile(std::string fileName);
    static void Destroy();

    /*
     *  Get methods
     */
    // Simulation
    double GetSimulationDuration() const;
    domain_type GetDomain() const;
    ionic_model_type GetIonicModel() const;
    
    bool GetCreateSlab() const;
    bool GetLoadMesh() const;
     
    void GetSlabDimensions(c_vector<double, 3>& slabDimensions) const;    
    double GetInterNodeSpace() const;
    
    std::string GetMeshName() const;
    
    void GetStimuli(std::vector<SimpleStimulus>& stimuliApplied, std::vector<ChasteCuboid>& stimulatedAreas) const;
    void GetCellHeterogeneities(std::vector<ChasteCuboid>& cellHeterogeneityAreas,
    							std::vector<double>& scaleFactorGks,
    							std::vector<double>& scaleFactorIto) const;
    void GetConductivityHeterogeneities(std::vector<ChasteCuboid>& conductivitiesHeterogeneityAreas,
				  					 	std::vector< c_vector<double,3> >& intraConductivities,
										std::vector< c_vector<double,3> >& extraConductivities) const;
    std::string GetOutputDirectory() const;

    // Physiological
    void GetIntracellularConductivities(c_vector<double, 3>& intraConductivities) const;
    void GetIntracellularConductivities(c_vector<double, 2>& intraConductivities) const;    
    void GetIntracellularConductivities(c_vector<double, 1>& intraConductivities) const;
        
    void GetExtracellularConductivities(c_vector<double, 3>& extraConductivities) const;
    void GetExtracellularConductivities(c_vector<double, 2>& extraConductivities) const;    
    void GetExtracellularConductivities(c_vector<double, 1>& extraConductivities) const;

    bool GetIsMediaOrthotropic() const;
    double GetSurfaceAreaToVolumeRatio() const;
    double GetCapacitance() const;

    // Numerical
    double GetOdeTimeStep() const;
    double GetPdeTimeStep() const;
    double GetPrintingTimeStep() const;

    bool GetUseAbsoluteTolerance() const;
    double GetAbsoluteTolerance() const;

    bool GetUseRelativeTolerance() const;
    double GetRelativeTolerance() const;

    ksp_solver_type GetKSPSolver() const;
    ksp_preconditioner_type GetKSPPreconditioner() const;


    /*
     *  Set methods
     */
    // Simulation
    void SetSimulationDuration(double simulationDuration);
    void SetDomain(domain_type domain);
    void SetIonicModel(ionic_model_type ionicModel);
    void SetMeshFileName(std::string meshPrefix);
    void SetOutputDirectory(std::string outputDirectory);

    // Physiological
    void SetIntracellularConductivities(const c_vector<double, 3>& intraConductivities);
    void SetIntracellularConductivities(const c_vector<double, 2>& intraConductivities);
    void SetIntracellularConductivities(const c_vector<double, 1>& intraConductivities);
        
    void SetExtracellularConductivities(const c_vector<double, 3>& extraConductivities);
    void SetExtracellularConductivities(const c_vector<double, 2>& extraConductivities);
    void SetExtracellularConductivities(const c_vector<double, 1>& extraConductivities);
    
    void SetMediaIsOrthotropic();
    void SetMediaIsAxisymmetric();
    void SetSurfaceAreaToVolumeRatio(double ratio);
    void SetCapacitance(double capacitance);

    // Numerical
    void SetOdePdeAndPrintingTimeSteps(double odeTimeStep, double pdeTimeStep, double printingTimeStep);
    void SetOdeTimeStep(double odeTimeStep);
    void SetPdeTimeStep(double pdeTimeStep);
    void SetPrintingTimeStep(double printingTimeStep);

    void SetUseRelativeTolerance(double relativeTolerance);
    void SetUseAbsoluteTolerance(double absoluteTolerance);

    void SetKSPSolver(ksp_solver_type kspSolver);
    void SetKSPPreconditioner(ksp_preconditioner_type kspPreconditioner);

protected:
    // Only to be accesed by the tests
    friend class TestHeartConfig;

    chaste_parameters_type* UserParameters();
    chaste_parameters_type* DefaultParameters();


private:
    HeartConfig();
    ~HeartConfig();

    chaste_parameters_type* mpUserParameters;
    chaste_parameters_type* mpDefaultParameters;

    /** The single instance of the class */
    static HeartConfig* mpInstance;

    // Misc
    template<class TYPE>
    TYPE* DecideLocation(TYPE* ptr1, TYPE* ptr2, const std::string& nameParameter) const;
    //Utility method to parse an XML parameters file
    chaste_parameters_type* ReadFile(std::string fileName);

};

#endif /*HEARTCONFIG_HPP_*/
