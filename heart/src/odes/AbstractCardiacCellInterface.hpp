/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTCARDIACCELLINTERFACE_HPP_
#define ABSTRACTCARDIACCELLINTERFACE_HPP_

#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"
#include "OdeSolution.hpp"

/**
 * This class defines a common interface to AbstractCardiacCell and AbstractCvodeCell,
 * primarily for the benefit of single-cell cardiac simulations, so that calling code
 * doesn't need to know which solver is being used.
 * 
 * Strictly speaking this isn't an interface, since some methods have implementations
 * defined.  But the name AbstractCardiacCell was already taken.
 * 
 * Note that serialization is not defined for this class.  AbstractCvodeCell does not
 * have (or need) serialization at present, so adding serialization here would break
 * archive backwards compatibility for little gain.
 */
class AbstractCardiacCellInterface
{
public:
    /**
     * Create a new cardiac cell. The state variables of the cell will be 
     * set to AbstractOdeSystemInformation::GetInitialConditions(). Note that
     * calls to SetDefaultInitialConditions() on a particular instance of this class
     * will not modify its state variables. You can modify them directly with 
     * rGetStateVariables(). 
     *
     * @param pOdeSolver  the ODE solver to use when simulating this cell
     *    (ignored for some subclasses; can be an empty pointer in these cases)
     * @param voltageIndex  the index of the transmembrane potential within the vector of state variables
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    AbstractCardiacCellInterface(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                                 unsigned voltageIndex,
                                 boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /** Virtual destructor */
    virtual ~AbstractCardiacCellInterface();
    
    /**
     * Reset the model's state variables to the default initial conditions.
     */
    virtual void ResetToInitialConditions()=0;
    
    /**
     * Set the timestep (or maximum timestep) to use for simulating this cell.
     *
     * @param dt  the timestep
     */
    virtual void SetTimestep(double dt)=0;

    /**
     * Simulate this cell's behaviour between the time interval [tStart, tEnd],
     * updating the internal state variable values.
     * The timestep used will depend on the subclass implementation.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void SolveAndUpdateState(double tStart, double tEnd)=0;
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * and return state variable values.  The timestep used will depend on the
     * subclass implementation.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     * @param tSamp  sampling interval for returned results (defaults to dt)
     */
    virtual OdeSolution Compute(double tStart, double tEnd, double tSamp=0.0)=0;

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * but does not update the voltage.  The timestep used will depend on the
     * subclass implementation.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void ComputeExceptVoltage(double tStart, double tEnd)=0;

    /**
     * Computes the total current flowing through the cell membrane, using the current
     * values of the state variables.
     *
     * Should return a value in units of microamps/cm^2.  Note that many cell models
     * do not use these dimensions (let alone these units) and so a complex conversion
     * is required.  There are 2 main cases:
     *   - Cell model uses amps per unit capacitance.  Often in this case the units used
     *     for the cell capacitance don't make sense (e.g. uF/cm^2 is used, and dV/dt=I/C_m).
     *     Hence we suggest examining the equation for dV/dt given in the model to determine
     *     what the cell model really considers the value for C_m to be, and scaling by
     *     Chaste's C_m / cell model C_m (the latter implicitly being dimensionless).
     *   - Cell model uses amps.  In this case you need to divide by an estimate of the cell
     *     surface area.  Assuming the model represents a single cell, and gives C_m in farads,
     *     then scaling by Chaste's C_m / model C_m seems reasonable.  If the 'cell model'
     *     doesn't actually represent a single whole cell, then you'll have to be more careful.
     * In both cases additional scaling may be required to obtain correct units once the
     * dimensions have been sorted out.
     *
     * Chaste's value for C_m can be obtained from HeartConfig::Instance()->GetCapacitance()
     * and is measured in uF/cm^2.
     * 
     * @param pStateVariables  optionally can be supplied to evaluate the ionic current at the
     *     given state; by default the cell's internal state will be used.
     */
    virtual double GetIIonic(const std::vector<double>* pStateVariables=NULL)=0;

    /** Set the transmembrane potential
     * @param voltage  new value
     */
    virtual void SetVoltage(double voltage)=0;

    /**
     * Get the current value of the transmembrane potential, as given
     * in our state variable vector.
     */
    virtual double GetVoltage()=0;

    /** Get the index of the transmembrane potential within our state variable vector. */
    unsigned GetVoltageIndex();

    /**
     * Set the intracellular stimulus.
     * Shorthand for SetIntracellularStimulusFunction.
     * @param pStimulus  new stimulus function
     */
    void SetStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus);

    /**
     * Get the value of the intracellular stimulus.
     * Shorthand for GetIntracellularStimulus.
     * @param time  the time at which to evaluate the stimulus
     */
    double GetStimulus(double time);

    /**
     * Set the intracellular stimulus.
     * This should have units of uA/cm^2 for single-cell problems,
     * or uA/cm^3 in a tissue simulation.
     * @param pStimulus  new stimulus function
     */
    void SetIntracellularStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus);

    /**
     * Get the value of the intracellular stimulus.
     * This will have units of uA/cm^2 for single-cell problems,
     * or uA/cm^3 in a tissue simulation.
     *
     * @param time  the time at which to evaluate the stimulus
     */
    double GetIntracellularStimulus(double time);

    /**
     * Get the value of the intracellular stimulus.
     * This will always be in units of uA/cm^2.
     *
     * @param time  the time at which to evaluate the stimulus
     */
    double GetIntracellularAreaStimulus(double time);

    /**
     * Set whether this cell object exists in the context of a tissue simulation,
     * or can be used for single cell simulations.  This affects the units of the
     * intracellular stimulus (see GetIntracellularStimulus) and so is used by
     * GetIntracellularAreaStimulus to perform a units conversion if necessary.
     *
     * @param tissue  true if cell is in a tissue
     */
    void SetUsedInTissueSimulation(bool tissue=true);

    /**
     * Use CellML metadata to set up the default stimulus for this cell.
     * By default this method will always throw an exception.  For suitably annotated
     * models, PyCml will override this to provide a RegularStimulus as defined in
     * the CellML.
     */
    virtual void UseCellMLDefaultStimulus();
    
    /**
     * @return Whether the cell was generated from a CellML file with stimulus metadata.
     */
    bool HasCellMLDefaultStimulus();

    /**
     *  Empty method which can be over-ridden in concrete cell class which should
     *  go through the current state vector and go range checking on the values
     *  (eg check that concentrations are positive and gating variables are between
     *  zero and one). This method is called in the ComputeExceptVoltage method.
     */
    virtual void VerifyStateVariables()
    {
        // See also #794.
    }
    
    /**
     * @return The Intracellular stimulus function pointer
     */
    boost::shared_ptr<AbstractStimulusFunction> GetStimulusFunction();
    
    /**
     * For boost archiving use only
     * (that's why the 'consts' are required)
     *
     * @return The Intracellular stimulus function pointer
     */
    const boost::shared_ptr<AbstractStimulusFunction> GetStimulusFunction() const;

    /**
     * For boost archiving use only
     * (that's why the 'consts' are required)
     *
     * @return pointer to the ODE solver being used
     */
    const boost::shared_ptr<AbstractIvpOdeSolver> GetSolver() const;
    
protected:
    /** The index of the voltage within our state variable vector. */
    unsigned mVoltageIndex;

    /** Pointer to the solver used to simulate this cell. */
    boost::shared_ptr<AbstractIvpOdeSolver> mpOdeSolver;

    /** The intracellular stimulus current. */
    boost::shared_ptr<AbstractStimulusFunction> mpIntracellularStimulus;

    /**
     * Flag set to true if ComputeExceptVoltage is called, to indicate
     * to subclass EvaluateYDerivatives methods that V should be
     * considered fixed, and hence dV/dt set to zero.
     */
    bool mSetVoltageDerivativeToZero;

    /** Whether this cell exists in a tissue, or is an isolated cell. */
    bool mIsUsedInTissue;

    /** Whether this cell has a default stimulus specified by CellML metadata */
    bool mHasDefaultStimulusFromCellML;

private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * This is needed to register the base-derived relationship between this class
     * and AbstractCardiacCell; however, serialization of our members is done by
     * AbstractCardiacCell in order to maintain archive backwards compatibility more
     * easily, since we never archive AbstractCvodeCell s.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
    }
};

CLASS_IS_ABSTRACT(AbstractCardiacCellInterface)


#endif /*ABSTRACTCARDIACCELLINTERFACE_HPP_*/
