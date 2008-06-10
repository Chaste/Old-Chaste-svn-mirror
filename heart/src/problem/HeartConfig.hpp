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
    void GetStimuli(std::vector<SimpleStimulus>& stimuliApplied, std::vector<ChasteCuboid>& stimulatedAreas) const;
    void GetCellHeterogeneities(std::vector<ChasteCuboid>& cellHeterogeneityAreas,
    							std::vector<double>& scaleFactorGks,
    							std::vector<double>& scaleFactorIto) const;
    void GetConductivityHeterogeneities(std::vector<ChasteCuboid>& conductivitiesHeterogeneityAreas,
				  					 	std::vector< c_vector<double,3> >& intraConductivities,
										std::vector< c_vector<double,3> >& extraConductivities) const;
    std::string GetOutputDirectory() const;

    // Physiological
    c_vector<double, 3> GetIntracellularConductivities() const;
    c_vector<double, 3> GetExtracellularConductivities() const;
    double GetSurfaceAreaToVolumeRatio() const;
    double GetCapacitance() const;

    // Numerical
    double GetOdeTimestep() const;
    double GetPdeTimestep() const;
    double GetPrintingTimestep() const;

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
    void SetOutputDirectory(std::string outputDirectory);

    // Physiological
    void SetIntracellularConductivities(const c_vector<double, 3>& intraConductivities);
    void SetExtracellularConductivities(const c_vector<double, 3>& extraConductivities);
    void SetSurfaceAreaToVolumeRatio(double ratio);
    void SetCapacitance(double capacitance);

    // Numerical
    void SetTimesteps(double odeTimestep, double pdeTimestep, double printingTimestep);
    void SetOdeTimestep(double odeTimestep);
    void SetPdeTimestep(double pdeTimestep);
    void SetPrintingTimestep(double printingTimestep);

    void SetTolerances(double relativeTolerance, double absoluteTolerance, ksp_use_type use);
    void SetUseRelativeTolerance(void);
    void SetUseAbsoluteTolerance(void);
    void SetRelativeTolerance(double relativeTolerance);
    void SetAbsoluteTolerance(double absoluteTolerance);

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
