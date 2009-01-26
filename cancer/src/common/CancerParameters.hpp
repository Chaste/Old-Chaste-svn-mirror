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
#ifndef CANCERPARAMETERS_HPP_
#define CANCERPARAMETERS_HPP_

#include <boost/serialization/access.hpp>
#include <cassert>

/**
 * A special singleton class which holds all of the parameters used in the cancer simulations
 *
 * Because this is a singleton class it can be called from whichever part of the code needs
 * to find out a parameter value, the structure is quite simple with default values given
 * upon initialisation and Set() and Get() methods for each parameter.
 *
 * For details of each parameter refer to the member variable documentation for this class
 * rather than the Get() and Set() function descriptions.
 */
class CancerParameters
{
public:

    /**
     * Call this method to access the global parameters holder.
     *
     * @return a single instance of the class
     */
    static CancerParameters* Instance();

    /**
     * Get methods
     */
    double GetStemCellG1Duration();
    double GetTransitCellG1Duration();
    double GetHepaOneCellG1Duration();
    double GetMinimumGapDuration();
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
    double GetHepaOneCellQuiescentConcentration();
    double GetWntTransitThreshold();
    double GetWntStemThreshold();
    double GetTopOfLinearWntConcentration();
    double GetCriticalHypoxicDuration();
    double GetCryptProjectionParameterA();
    double GetCryptProjectionParameterB();
    double GetApoptoticSpringTensionStiffness();
    double GetApoptoticSpringCompressionStiffness();
    double GetWntChemotaxisStrength();
    double GetSymmetricDivisionProbability();
    double GetAreaBasedDampingConstantParameter();

    /**
     * Set methods
     */
    void SetStemCellG1Duration(double);
    void SetTransitCellG1Duration(double);
    void SetHepaOneCellG1Duration(double);
    void SetMinimumGapDuration(double);
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
    void SetHepaOneCellQuiescentConcentration(double);
    void SetWntTransitThreshold(double);
    void SetWntStemThreshold(double);
    void SetTopOfLinearWntConcentration(double);
    void SetCriticalHypoxicDuration(double);
    void SetHepaOneParameters();
    void SetCryptProjectionParameterA(double);
    void SetCryptProjectionParameterB(double);
    void SetApoptoticSpringTensionStiffness(double);
    void SetApoptoticSpringCompressionStiffness(double);
    void SetWntChemotaxisStrength(double);
    void SetSymmetricDivisionProbability(double);
    void SetAreaBasedDampingConstantParameter(double);

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
     * Duration of G1 phase for stem cells.
     * May be used as a mean duration for stochastic cell cycle models.
     *
     */
    double mStemCellG1Duration;

    /**
     * Duration of G1 phase for transit cells.
     * May be used as a mean duration for stochastic cell cycle models.
     */
    double mTransitCellG1Duration;

    /**
     * Duration of G1 phase for HEPA-1 cells, for use in monolayer/spheroid simulations.
     * May be used as a mean duration for stochastic cell cycle models.
     */
    double mHepaOneCellG1Duration;

    /**
     * Minimum possbile duration of either of the gap phases (G1 or G2).
     * Used to guarantee a strictly positive duration in cell cycle models that
     * use normal random deviates for G1 or G2 phases.
     */
    double mMinimumGapDuration;

    /**
     * Duration of S phase for all cell types.
     */
    double mSDuration;

    /**
     * Duration of G2 phase for all cell types.
     */
    double mG2Duration;

    /**
     * Duration of M phase for all cell types.
     */
    double mMDuration;

    /**
     * How many generations a transit cell lives for before becoming fully differentiated.
     */
    unsigned mMaxTransitGenerations;

    /**
     * The length of the crypt, non-dimensionalised with cell length.
     * This parameter determines when cells are sloughed from the crypt.
     */
    double mCryptLength;

    /**
    * The width of the crypt, non-dimensionalised with cell length.
    * This determines when cells are sloughed from the crypt in 2D.
    */
    double mCryptWidth;

    /**
     * Spring stiffness.
     * Represented by the parameter mu in the model by Meineke et al (2001).
     */
    double mSpringStiffness;

    /**
     * Damping constant for normal cells.
     * Represented by the parameter eta in the model by Meineke et al (2001).
     */
    double mDampingConstantNormal;

    /**
     * Damping constant for mutant cells.
     */
    double mDampingConstantMutant;

    /**
     * Scaling factor for beta catenin to spring strength
     */
    double mBetaCatSpringScaler;

    /**
     * The time it takes for a cell to fully undergo apoptosis
     */
    double mApoptosisTime;

    /**
     * Initial separation placement of mother/daughter cells at birth
     */
    double mDivisionSeparation;

    /**
     * Initial resting spring length after cell division.
     * The value of thiis parameter should be larger than mDivisionSeparation,
     * because of pressure from neighbouring springs.
     */
    double mDivisionRestingSpringLength;

    /**
     * Non-dimensionalized oxygen concentration below which HEPA-1 cells are
     * considered to be hypoxic.
     * A prolonged period of hypoxia causes the cell to become apoptotic.
     */
    double mHepaOneCellHypoxicConcentration;

    /**
     * Non-dimensionalized oxygen concentration below which HEPA-1 cells are
     * considered to be quiescent and slow their progress through the G1 phase
     * of the cell cycle.
     */
    double mHepaOneCellQuiescentConcentration;

    /**
     * Non-dimensionalized Wnt threshold, above which cells progress through the cell cycle.
     */
    double mWntTransitThreshold;

    /**
     * Non-dimensionalized Wnt threshold, above which cells behave as stem cells.
     */
    double mWntStemThreshold;

    /**
     * The proportion of the crypt that has a Wnt gradient.
     * The Wnt concentration goes to zero at this height up the crypt.
     */
    double mTopOfLinearWntConcentration;

    /**
     * Non-dimensionalized critical hypoxic duration.
     */
    double mCriticalHypoxicDuration;

    /**
     * Parameter a, for use in crypt projection simulations, in which the crypt
     * surface is given in cylindrical polar coordinates by z = a*r^b.
     */
    double mCryptProjectionParameterA;

    /**
     * Parameter b, for use in crypt projection simulations, in which the crypt
     * surface is given in cylindrical polar coordinates by z = a*r^b.
     */
    double mCryptProjectionParameterB;

    /**
     * Non-dimensionalized 'stiffness' of a apoptotic cell under tension.
     */
    double mApoptoticSpringTensionStiffness;

    /**
     * Non-dimensionalized 'stiffness' of a apoptotic cell under compression.
     */
    double mApoptoticSpringCompressionStiffness;

    /**
     * Strength of Wnt-based chemotactic force.
     */
    double mWntChemotaxisStrength;

    /**
     * Probability of symmetric division.
     */
    double mSymmetricDivisionProbability;

    /**
     * Non-dimensional parameter d0 for use in area-based damping constant calculations.
     */
    double mAreaBasedDampingConstantParameter;

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
        archive & mMinimumGapDuration;
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
        archive & mHepaOneCellQuiescentConcentration;
        archive & mWntTransitThreshold;
        archive & mWntStemThreshold;
        archive & mTopOfLinearWntConcentration;
        archive & mCriticalHypoxicDuration;
        archive & mCryptProjectionParameterA;
        archive & mCryptProjectionParameterB;
        archive & mApoptoticSpringTensionStiffness;
        archive & mApoptoticSpringCompressionStiffness;
        archive & mWntChemotaxisStrength;
        archive & mSymmetricDivisionProbability;
        archive & mAreaBasedDampingConstantParameter;
    }
};


#endif /*CANCERPARAMETERS_HPP_*/
