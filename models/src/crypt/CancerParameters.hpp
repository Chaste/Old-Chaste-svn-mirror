#ifndef CANCERPARAMETERS_HPP_
#define CANCERPARAMETERS_HPP_

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
    unsigned GetMaxTransitGenerations();
    double GetCryptLength();
    double GetCryptWidth();
    double GetMeinekeLambda();
    double GetAlpha();
    double GetNaturalSpringLength();
    
    void SetStemCellCycleTime(double);
    void SetTransitCellCycleTime(double);
    void SetMaxTransitGenerations(unsigned);
    void SetCryptLength(double);
    void SetCryptWidth(double);
    void SetMeinekeLambda(double);
    void SetNaturalSpringLength(double);
    
protected:
    CancerParameters();
    CancerParameters(const CancerParameters&);
    CancerParameters& operator= (const CancerParameters&);
    
private:

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
     * The non-dimensionalised parameter.
     * alpha = stem cell cycle time * Meineke lambda.
     */
    double mAlpha;
    
    /**
     * The resting length of springs connected mature cells.
     * (i.e. the size of mature cells)
     */
    double mNaturalSpringLength;
};

#endif /*CANCERPARAMETERS_HPP_*/
