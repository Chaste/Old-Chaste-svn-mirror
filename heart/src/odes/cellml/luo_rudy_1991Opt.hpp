#ifndef CELLLUO_RUDY_1991FROMCELLMLOPT_HPP_
#define CELLLUO_RUDY_1991FROMCELLMLOPT_HPP_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: luo_rudy_1991
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translate: 7848, pycml: 7844)
//! on Sat Feb  6 09:54:39 2010
//! 
//! <autogenerated>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"

class Cellluo_rudy_1991FromCellMLOpt : public AbstractCardiacCell
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacCell >(*this);
    }
    
    // 
    // Settable parameters and readable variables
    // 
    double var_membrane__I_stim;
    double var_membrane__i_Na;
    double var_membrane__i_si;
    double var_membrane__i_K;
    double var_membrane__i_K1;
    double var_membrane__i_Kp;
    double var_membrane__i_b;
    
public:
    double Get_membrane__I_stim();
    double Get_membrane__i_Na();
    double Get_membrane__i_si();
    double Get_membrane__i_K();
    double Get_membrane__i_K1();
    double Get_membrane__i_Kp();
    double Get_membrane__i_b();
    Cellluo_rudy_1991FromCellMLOpt(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~Cellluo_rudy_1991FromCellMLOpt();
    void VerifyStateVariables();
    double GetIIonic();
    void EvaluateYDerivatives(double var_environment__time, const std::vector<double>& rY, std::vector<double>& rDY);
};


// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Cellluo_rudy_1991FromCellMLOpt)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const Cellluo_rudy_1991FromCellMLOpt * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }
        
        template<class Archive>
        inline void load_construct_data(
            Archive & ar, Cellluo_rudy_1991FromCellMLOpt * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)Cellluo_rudy_1991FromCellMLOpt(p_solver, p_stimulus);
        }
        
    }
    
}

#endif // CELLLUO_RUDY_1991FROMCELLMLOPT_HPP_
