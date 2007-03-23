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
    
    double GetStemCellCycleTime();
    double GetTransitCellCycleTime();
    double GetSG2MDuration();
    unsigned GetMaxTransitGenerations();
    double GetCryptLength();
    double GetCryptWidth();
    double GetMeinekeLambda();
    double GetApoptosisTime();
    
    void SetStemCellCycleTime(double);
    void SetTransitCellCycleTime(double);
    void SetSG2MDuration(double);
    void SetMaxTransitGenerations(unsigned);
    void SetCryptLength(double);
    void SetCryptWidth(double);
    void SetMeinekeLambda(double);
    void SetApoptosisTime(double);
    
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
    double mStemCellCycleTime;
    /**
     * Transit cell cycle time.
     * May be used as a mean time for stochastic cell cycle models.
     * Should probably be non-dimensionalised with stem cell cycle time (ticket:204)
     */
    double mTransitCellCycleTime;
    
    /**
     * S-G2-M Phase Duration.
     * Used by the Wnt signalling model which only models the G1 phase.
     */
    double mSG2MDuration;
    
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
     * The ratio of mu (spring constant) to eta (damping constant).
     * Has dimensions 1/time.
     */
    double mMeinekeLambda;
    
    /**
     * The time it takes to fully undergo apoptosis
     */
    double mApoptosisTime;
    
    friend class boost::serialization::access;
    /**
     * As with other singleton classes, ensure the instance of this
     * class is serialized directly before being serialized via a
     * pointer.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mStemCellCycleTime;
        archive & mTransitCellCycleTime;
        archive & mSG2MDuration;
        archive & mMaxTransitGenerations;
        archive & mCryptLength;
        archive & mCryptWidth;
        archive & mMeinekeLambda;
        archive & mApoptosisTime;
    }
};


#endif /*CANCERPARAMETERS_HPP_*/
