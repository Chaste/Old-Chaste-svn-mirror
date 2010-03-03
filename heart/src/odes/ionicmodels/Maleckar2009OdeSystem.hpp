#ifndef _MALECKAR2009ODESYSTEM_
#define _MALECKAR2009ODESYSTEM_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: Maleckar
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translate: 7385, pycml: 7152)
//! on Mon Dec 14 15:22:44 2009
//! 
//! <autogenerated>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <cmath>
#include <cassert>
#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "OdeSystemInformation.hpp"
#include "AbstractStimulusFunction.hpp"
/**
 * 
 * This class implements the Maleckar 2009 atrial cell model
 * The code has been translated by PyCml.
 * Some scale factors have been manually added to be used to simulate the effect of some drug.
 * By default, it behaves in a 'control' situation i.e., like the file downloaded from the CellML repository.
 * State variables are initialised as to be at steady state for 1Hz.
 * 
 */
class Maleckar2009OdeSystem : public AbstractCardiacCell
{
    
private:
    /** Scale factor for Gks, set to 1 by the constructor*/
    double mScaleFactorGks;
    /** Scale factor for Gto, set to 1 by the constructor*/
    double mScaleFactorIto;
    /** Scale factor for Gkr, set to 1 by the constructor*/
    double mScaleFactorGkr;
    /** Scale factor for Gna, set to 1 by the constructor*/
    double mScaleFactorGna;
    /** Scale factor for Ach, set to 1e-24 by the constructor*/
    double mScaleFactorAch;
    /** Scale factor for GNaK,  set to 1 by the constructor*/
    double mScaleFactorGNaK;
    /** Scale factor for GNaCa, set to 1 by the constructor*/
    double mScaleFactorGNaCa;
    /** Scale factor for GCaL, set to 1 by the constructor*/
    double mScaleFactorGCaL;
    /** Scale factor for GKur, set to 1 by the constructor*/
    double mScaleFactorGKur;
    /** Scale factor for GK1, set to 1 by the constructor*/
    double mScaleFactorGK1;
    /** Scale factor for AZD, set to 0 by the constructor*/
    double mScaleFactorAZD;
    
    /**
     *  Range-checking on the current values of the state variables. Make sure
     *  all gating variables have are within zero and one, and all concentrations
     *  are positive. 
     * 
     *  ///\ TODO : implement this!
     */
    void VerifyStateVariables();
    
    /** Needed for serialization. */
    friend class boost::serialization::access;
    
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
    }
    
public:

    /**
     * Constructor
     * 
     * @param pSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    Maleckar2009OdeSystem(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, 
                 boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
                 
    /**
     * Destructor
     */   
    ~Maleckar2009OdeSystem(void);

    /**
     * Set the scale factor for Gks
     * @param sfgks
     */
    void SetScaleFactorGks(double sfgks);
    
    /**
     * Set the scale factor for Ito
     * @param sfito
     */
    void SetScaleFactorIto(double sfito);
    
    /**
     * Set the scale factor for GKr
     * @param sfgkr
     */
    void SetScaleFactorGkr(double sfgkr);
    
    /**
     * Set the scale factor for GNa
     * @param sfgna
     */
    void SetScaleFactorGna(double sfgna);
    
    /**
     * Set the scale factor for Ach, acetilcholine concentration
     * @param sfach
     */
    void SetScaleFactorAch(double sfach);
    
    /**
     * Set the scale factor for Na/K exchanger
     * @param sfgnak
     */
    void SetScaleFactorGNaK(double sfgnak);

    /**
     * Set the scale factor for Na/Ca exchanger
     * @param sfgnaca
     */
    void SetScaleFactorGNaCa(double sfgnaca);

    /**
     * Set the scale factor for L-type current
     * @param sfgcal
     */
    void SetScaleFactorGCaL(double sfgcal);

    /**
     * Set the scale factor for the Kur current
     * @param sfgkur
     */   
    void SetScaleFactorGKur(double sfgkur);
    
    /**
     * Set the scale factor for the K1 current
     * @param sfgk1
     */       
    void SetScaleFactorGK1(double sfgk1);

    /**
     * Set the scale factor for concentration of AZD compound (fom AstraZeneca)
     * @param sfazd
     */   
    void SetScaleFactorAZD(double sfazd);

    /**
     * Returns the ionic current
     * 
     * @return the total ionic current
     */
    double GetIIonic();

    /**
     * Fill in a vector representing the RHS of the TenTusscher2006 system
     * of Odes at each time step, y' = [y1' ... yn'].
     * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */    
    void EvaluateYDerivatives(
            double var_environment__time,
            const std::vector<double> &rY,
            std::vector<double> &rDY);
    
};

// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Maleckar2009OdeSystem)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const Maleckar2009OdeSystem * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }
        
        template<class Archive>
        inline void load_construct_data(
            Archive & ar, Maleckar2009OdeSystem * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)Maleckar2009OdeSystem(p_solver, p_stimulus);
        }
        
    }
    
}

#endif
