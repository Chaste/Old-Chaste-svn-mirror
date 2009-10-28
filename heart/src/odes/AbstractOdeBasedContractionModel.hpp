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

#ifndef ABSTRACTODEBASEDCONTRACTIONMODEL_
#define ABSTRACTODEBASEDCONTRACTIONMODEL_

#include "AbstractOdeSystem.hpp"

/**
 *  Options for the different contraction models (both stretch-dependent and independent, ditto stretch-rate)
 *  that has been implemented
 */
typedef enum ContractionModel_
{
    KERCHOFFS2003,
    NHS
} ContractionModel;


/**
 *  Struct storing the input parameters that might be used by a contraction model (excl stretch and stretch-rate,
 *  as these may be set several times using the current deformation guess by the implicit assembler).
 */ 
typedef struct ContractionModelInputParameters_
{
    double voltage;                           /**< Input voltage (mV)*/
    double intracellularCalciumConcentration; /**< Input calcium concentration (mMol) */
    double time;                              /**< Input time (ms) (for time dependent contraction models) */
} ContractionModelInputParameters;
    

/**
 *  Abstract base class for ODE-based contraction models. Inherits from AbstractOdeSystem and defines
 *  a contraction-model interface.
 */
class AbstractOdeBasedContractionModel : public AbstractOdeSystem
{
public:
    /**
     *  Constructor does nothing except pass through the number of state variables
     *  @param numStateVariables Number of state variables in the ODEs
     */
    AbstractOdeBasedContractionModel(unsigned numStateVariables)
        : AbstractOdeSystem(numStateVariables)
    {
    }
    
    /** 
     *  Does the model depend on the stretch. (Pure, to be implemented in the concrete class).
     */
    virtual bool IsStretchDependent()=0;

    /** 
     *  Does the model depend on the stretch-rate. (Pure, to be implemented in the concrete class).
     */
    virtual bool IsStretchRateDependent()=0;

    /** 
     *  Set any input parameters (excl stretch and stretch rate). (Pure, to be implemented in the concrete class).
     *
     *  @param rInputParameters  contains various parameters: voltage, intracellular calcium concentration and 
     *  time (at next timestep)
     */
    virtual void SetInputParameters(ContractionModelInputParameters& rInputParameters)=0;

    /** 
     *  Set the stretch and stretch rate. (Pure, to be implemented in the concrete class).
     * 
     *  @param stretch  fibre stretch (dimensionless)
     *  @param stretchRate  fibre stretch rate (1/ms)
     */
    virtual void SetStretchAndStretchRate(double stretch, double stretchRate)=0;
    
    /** Safe setting of stretch-only, for stretch-rate independent models ONLY
     *  @oaram stretch Stretch in fibre direction
     */
    void SetStretch(double stretch)
    {
        assert(!IsStretchRateDependent());
        SetStretchAndStretchRate(stretch, 0.0);
    }

    /** 
     *  Get the current active tension (note, actually a stress). (Pure, to be implemented in the concrete class).
     */
    virtual double GetActiveTension()=0;
};


#endif /*ABSTRACTODEBASEDSTRETCHINDEPENDENTCONTRACTIONMODEL_*/
