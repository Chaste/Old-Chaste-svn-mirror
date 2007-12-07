#ifndef CANCERPARAMETERS_HPP_
#define CANCERPARAMETERS_HPP_

#include <boost/serialization/access.hpp>
#include <cassert>

class CancerParameters
{
public:
    /**
     * Call this method to access the global parameters holder.
     * 
     * @return a single instance of the class
     */
    static CancerParameters* Instance();
    
    double GetStemCellG1Duration();
    double GetTransitCellG1Duration();
    double GetHepaOneCellG1Duration();
    double GetSG2MDuration();    
    double GetSDuration();
    double GetG2Duration();
    double GetMDuration();    
    unsigned GetMaxTransitGenerations();
    double GetCryptLength();
    double GetCryptWidth();
    double GetSpringStiffness();
    double GetDampingConstantNormal();
    double GetDampingConstantMutant();
    double GetBetaCatSpringScaler();
    double GetApoptosisTime();
    double GetDivisionRestingSpringLength();
    double GetDivisionSeparation();
    double GetHepaOneCellHypoxicConcentration();
    double GetRadialWntThreshold();
    double GetCriticalHypoxicDuration();
    
    void SetStemCellG1Duration(double);
    void SetTransitCellG1Duration(double);
    void SetHepaOneCellG1Duration(double);    
    void SetSDuration(double);
    void SetG2Duration(double);
    void SetMDuration(double);    
    void SetMaxTransitGenerations(unsigned);
    void SetCryptLength(double);
    void SetCryptWidth(double);
    void SetSpringStiffness(double);
    void SetDampingConstantNormal(double);
    void SetDampingConstantMutant(double);
    void SetBetaCatSpringScaler(double);
    void SetApoptosisTime(double);
    void SetDivisionRestingSpringLength(double);
    void SetDivisionSeparation(double);
    void SetHepaOneCellHypoxicConcentration(double);
    void SetRadialWntThreshold(double); 
    void SetCriticalHypoxicDuration(double);
    void SetHepaOneParameters();
      
    /** 
     *  Reset all parameters to their defaults
     */
    void Reset();
    
protected:
    CancerParameters();
    CancerParameters(const CancerParameters&);
    CancerParameters& operator= (const CancerParameters&);
    
private:
    /** The single instance of the class */
    static CancerParameters *mpInstance;
    
    /**
     * Stem cell cycle time, used to non-dimensionalise the problem
     */
    double mStemCellG1Duration;
    /**
     * Transit cell cycle time.
     * May be used as a mean time for stochastic cell cycle models.
     * Should probably be non-dimensionalised with stem cell cycle time (ticket:204)
     */
    double mTransitCellG1Duration;    
    /**
     * HEPA-1 cell cycle time. 
     * For use in monolayer/spheroid simulations.
     * May be used as a mean time for stochastic cell cycle models.
     */
    double mHepaOneCellG1Duration;  
      
    /**
     * S Phase Duration, currently for all cell cycle models except T&N.
     */
    double mSDuration;
    
    /**
     * G2 Phase Duration.
     * Used by the cell cycle models.
     */
    double mG2Duration;
    
    /**
     * M Phase Duration.
     * Used by the cell cycle models, and the mDivisionPairs methods.
     */
    double mMDuration;
    
    /**
     * How many generations a transit cell lives for before becoming fully differentiated.
     */
    unsigned mMaxTransitGenerations;
    
    /**
     * The non-dimensionalised (with cell length) length of the crypt.
     * This determines when cells are sloughed from the crypt.
     */
    double mCryptLength;
    
    /**
    * The non-dimensionalised (with cell length) width of the crypt.
    * This determines when cells are sloughed from the crypt. in 2D
    */
    double mCryptWidth;
    
    /**
     * Spring stiffness represented by mu in Meineke
     */
    double mSpringStiffness;
    
    /**
     * Damping constant for normal cells, eta in Meineke
     */
    double mDampingConstantNormal;
    
    /**
     * Damping constant for mutant cells, eta in Meineke
     */
    double mDampingConstantMutant;
    
    /**
     * Scaling factor for beta catenin to spring strength
     */
    double mBetaCatSpringScaler;
    
    /**
     * The time it takes to fully undergo apoptosis
     */
    double mApoptosisTime;
    
    /**
     * Initial separation placement of mother/daughter at birth
     */
    double mDivisionSeparation;
    
    /**
     * Initial resting spring length after division (should be longer than 
     * mDivisionSeparation because of pressure from neighbouring springs)
     */
    double mDivisionRestingSpringLength;
    
    /**
     * Non-dimensionalized oxygen concentration below which HEPA-1 are hypoxic 
     */
    double mHepaOneCellHypoxicConcentration;
      
    /**
     * Non-dimensionalized radial Wnt threshold, below which cells become TRANSIT
     */
    double mRadialWntThreshold;
    
    /**
     * Non-dimensionalized critical hypoxic duration
     */
    double mCriticalHypoxicDuration;
      
    friend class boost::serialization::access;
    /**
     * As with other singleton classes, ensure the instance of this
     * class is serialized directly before being serialized via a
     * pointer.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mStemCellG1Duration;
        archive & mTransitCellG1Duration;
        archive & mHepaOneCellG1Duration;       
        archive & mSDuration;
        archive & mG2Duration;
        archive & mMDuration;                       
        archive & mMaxTransitGenerations;
        archive & mCryptLength;
        archive & mCryptWidth;
        archive & mSpringStiffness;
        archive & mDampingConstantNormal;
        archive & mDampingConstantMutant;
        archive & mBetaCatSpringScaler;
        archive & mApoptosisTime;
        archive & mHepaOneCellHypoxicConcentration;
        archive & mRadialWntThreshold;
        archive & mCriticalHypoxicDuration;
    }
};


#endif /*CANCERPARAMETERS_HPP_*/
