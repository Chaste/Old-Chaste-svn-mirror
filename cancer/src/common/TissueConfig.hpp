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
#ifndef TISSUECONFIG_HPP_
#define TISSUECONFIG_HPP_

#include <boost/serialization/access.hpp>
#include <cassert>
#include "Exception.hpp"

/**
 * A special singleton class which holds all of the parameters used in the cancer simulations.
 *
 * Because this is a singleton class it can be called from whichever part of the code needs
 * to find out a parameter value, the structure is quite simple with default values given
 * upon initialisation and Set() and Get() methods for each parameter.
 *
 * For details of each parameter refer to the member variable documentation for this class
 * rather than the Get() and Set() function descriptions.
 */
class TissueConfig
{
public:

    /**
     * Call this method to access the global parameters holder.
     *
     * @return a single instance of the class
     */
    static TissueConfig* Instance();

    /**
     * @return mStemCellG1Duration
     */
    double GetStemCellG1Duration();
    /**
     * @return mTransitCellG1Duration
     */
    double GetTransitCellG1Duration();
    /**
     * @return mHepaOneCellG1Duration
     */
    double GetHepaOneCellG1Duration();
    /**
     * @return mMinimumGapDuration
     */
    double GetMinimumGapDuration();
    /**
     * @return mSDuration + mG2Duration + mMDuration
     */
    double GetSG2MDuration();
    /**
     * @return mSDuration
     */
    double GetSDuration();
    /**
     * @return mG2Duration
     */
    double GetG2Duration();
    /**
     * @return mMDuration
     */
    double GetMDuration();
    /**
     * @return mMaxTransitGenerations
     */
    unsigned GetMaxTransitGenerations();
    /**
     * @return mCryptLength
     */
    double GetCryptLength();
    /**
     * @return mCryptWidth
     */
    double GetCryptWidth();
    /**
     * @return mSpringStiffness
     */
    double GetSpringStiffness();
    /**
     * @return mMechanicsCutOffLength
     */
    double GetMechanicsCutOffLength();
    /**
     * @return mDampingConstantNormal
     */
    double GetDampingConstantNormal();
    /**
     * @return mDampingConstantMutant
     */
    double GetDampingConstantMutant();
    /**
     * @return mBetaCatSpringScaler
     */
    double GetBetaCatSpringScaler();
    /**
     * @return mApoptosisTime
     */
    double GetApoptosisTime();
    /**
     * @return mDivisionRestingSpringLength
     */
    double GetDivisionRestingSpringLength();
    /**
     * @return mDivisionSeparation
     */
    double GetDivisionSeparation();
    /**
     * @return mHepaOneCellHypoxicConcentration
     */
    double GetHepaOneCellHypoxicConcentration();
    /**
     * @return mHepaOneCellQuiescentConcentration
     */
    double GetHepaOneCellQuiescentConcentration();
    /**
     * @return mWntTransitThreshold
     */
    double GetWntTransitThreshold();
    /**
     * @return mWntStemThreshold
     */
    double GetWntStemThreshold();
    /**
     * @return mTopOfLinearWntConcentration
     */
    double GetTopOfLinearWntConcentration();
    /**
     * @return mCriticalHypoxicDuration
     */
    double GetCriticalHypoxicDuration();
    /**
     * @return mCryptProjectionParameterA
     */
    double GetCryptProjectionParameterA();
    /**
     * @return mCryptProjectionParameterB
     */
    double GetCryptProjectionParameterB();
    /**
     * @return mApoptoticSpringTensionStiffness
     */
    double GetApoptoticSpringTensionStiffness();
    /**
     * @return mApoptoticSpringCompressionStiffness
     */
    double GetApoptoticSpringCompressionStiffness();
    /**
     * @return mWntChemotaxisStrength
     */
    double GetWntChemotaxisStrength();
    /**
     * @return mSymmetricDivisionProbability
     */
    double GetSymmetricDivisionProbability();
    /**
     * @return mAreaBasedDampingConstantParameter
     */
    double GetAreaBasedDampingConstantParameter();
    /**
     * @return mMatureCellTargetArea
     */
    double GetMatureCellTargetArea();
    /**
     * @return mDeformationEnergyParameter
     */
    double GetDeformationEnergyParameter();
    /**
     * @return mMembraneSurfaceEnergyParameter
     */
    double GetMembraneSurfaceEnergyParameter();
    /**
     * @return mCellCellAdhesionEnergyParameter
     */
    double GetCellCellAdhesionEnergyParameter();
    /**
     * @return mCellBoundaryAdhesionEnergyParameter
     */
    double GetCellBoundaryAdhesionEnergyParameter();
	/**
	 * @return mOutputCellIdData
	 */
	bool GetOutputCellIdData();
	/**
	 * @return mOutputCellMutationStates
	 */
	bool GetOutputCellMutationStates();
	/**
	 * @return mOutputCellAncestors
	 */
    bool GetOutputCellAncestors();
    /**
     * @return mOutputCellTypes
     */
    bool GetOutputCellTypes();
    /**
     * @return mOutputCellVariables
     */
    bool GetOutputCellVariables();
    /**
     * @return mOutputCellCyclePhases
     */
    bool GetOutputCellCyclePhases();
    /**
     * @return mOutputCellAges
     */
    bool GetOutputCellAges();
    /**
     * @return mOutputCellAreas
     */
    bool GetOutputCellAreas();

    /**
     * Set mStemCellG1Duration.
     */
    void SetStemCellG1Duration(double);
    /**
     * Set mTransitCellG1Duration.
     */
    void SetTransitCellG1Duration(double);
    /**
     * Set mHepaOneCellG1Duration.
     */
    void SetHepaOneCellG1Duration(double);
    /**
     * Set mMinimumGapDuration.
     */
    void SetMinimumGapDuration(double);
    /**
     * Set mSDuration.
     */
    void SetSDuration(double);
    /**
     * Set mG2Duration.
     */
    void SetG2Duration(double);
    /**
     * Set mMDuration.
     */
    void SetMDuration(double);
    /**
     * Set mMaxTransitGenerations.
     */
    void SetMaxTransitGenerations(unsigned);
    /**
     * Set mCryptLength.
     */
    void SetCryptLength(double);
    /**
     * Set mCryptWidth.
     */
    void SetCryptWidth(double);
    /**
     * Set mSpringStiffness.
     */
    void SetSpringStiffness(double);
    /**
     * Set mMechanicsCutOffLength.
     */
    void SetMechanicsCutOffLength(double);
    /**
     * Set mDampingConstantNormal.
     */
    void SetDampingConstantNormal(double);
    /**
     * Set mDampingConstantMutant.
     */
    void SetDampingConstantMutant(double);
    /**
     * Set mBetaCatSpringScaler.
     */
    void SetBetaCatSpringScaler(double);
    /**
     * Set mApoptosisTime.
     */
    void SetApoptosisTime(double);
    /**
     * Set mDivisionRestingSpringLength.
     */
    void SetDivisionRestingSpringLength(double);
    /**
     * Set mDivisionSeparation.
     */
    void SetDivisionSeparation(double);
    /**
     * Set mHepaOneCellHypoxicConcentration.
     */
    void SetHepaOneCellHypoxicConcentration(double);
    /**
     * Set mHepaOneCellQuiescentConcentration.
     */
    void SetHepaOneCellQuiescentConcentration(double);
    /**
     * Set mWntTransitThreshold.
     */
    void SetWntTransitThreshold(double);
    /**
     * Set mWntStemThreshold.
     */
    void SetWntStemThreshold(double);
    /**
     * Set mTopOfLinearWntConcentration.
     */
    void SetTopOfLinearWntConcentration(double);
    /**
     * Set mCriticalHypoxicDuration.
     */
    void SetCriticalHypoxicDuration(double);
    /**
     * Set mHepaOneParameters.
     */
    void SetHepaOneParameters();
    /**
     * Set mCryptProjectionParameterA.
     */
    void SetCryptProjectionParameterA(double);
    /**
     * Set mCryptProjectionParameterB.
     */
    void SetCryptProjectionParameterB(double);
    /**
     * Set mApoptoticSpringTensionStiffness.
     */
    void SetApoptoticSpringTensionStiffness(double);
    /**
     * Set mApoptoticSpringCompressionStiffness.
     */
    void SetApoptoticSpringCompressionStiffness(double);
    /**
     * Set mWntChemotaxisStrength.
     */
    void SetWntChemotaxisStrength(double);
    /**
     * Set mSymmetricDivisionProbability.
     */
    void SetSymmetricDivisionProbability(double);
    /**
     * Set mAreaBasedDampingConstantParameter.
     */
    void SetAreaBasedDampingConstantParameter(double);
    /**
     * Set mMatureCellTargetArea.
     */
    void SetMatureCellTargetArea(double);
    /**
     * Set mDeformationEnergyParameter.
     */
    void SetDeformationEnergyParameter(double);
    /**
     * Set mMembraneSurfaceEnergyParameter.
     */
    void SetMembraneSurfaceEnergyParameter(double);
    /**
     * Set mCellCellAdhesionEnergyParameter.
     */
    void SetCellCellAdhesionEnergyParameter(double);
    /**
     * Set mCellBoundaryAdhesionEnergyParameter.
     */
    void SetCellBoundaryAdhesionEnergyParameter(double);
	/**
	 * Set mOutputCellIdData.
	 */
	void SetOutputCellIdData(bool);
	/**
	 * Set mOutputCellMutationStates.
	 */
    void SetOutputCellMutationStates(bool);
	/**
	 * Set mOutputCellAncestors.
	 */
    void SetOutputCellAncestors(bool);
	/**
	 * Set mOutputCellTypes.
	 */
    void SetOutputCellTypes(bool);
	/**
	 * Set mOutputCellVariables.
	 */
    void SetOutputCellVariables(bool);
	/**
	 * Set mOutputCellCyclePhases.
	 */
    void SetOutputCellCyclePhases(bool);
	/**
	 * Set mOutputCellAges.
	 */
    void SetOutputCellAges(bool);
	/**
	 * Set mOutputCellAreas.
	 */
    void SetOutputCellAreas(bool);

    /**
     *  Reset all parameters to their defaults
     */
    void Reset();

protected:

    /**
     * Default constructor.
     */
    TissueConfig();

    /**
     * Copy constructor.
     */
    TissueConfig(const TissueConfig&);

    /**
     * Overloaded assignement operator.
     */
    TissueConfig& operator= (const TissueConfig&);

private:

    /** The single instance of the class */
    static TissueConfig *mpInstance;

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
     * 
     * Represented by the parameter mu in the model by Meineke et al (2001) in 
     * their off-lattice model of the intestinal crypt 
     * (doi:10.1046/j.0960-7722.2001.00216.x).
     */
    double mSpringStiffness;
    
    /**
     * Mechanics cut off length.
     * Used in NodeBasedTissue.
     */
    double mMechanicsCutOffLength;

    /**
     * Damping constant for normal cells.
     * Represented by the parameter eta in the model by Meineke et al (2001) in 
     * their off-lattice model of the intestinal crypt 
     * (doi:10.1046/j.0960-7722.2001.00216.x).
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

    /**
     * Non-dimensional target area of a mature (fully-grown) TissueCell.
     * For use in vertex-based models.
     */
    double mMatureCellTargetArea;

    /**
     * Cell deformation energy parameter.
     * For use in vertex-based models.
     */
    double mDeformationEnergyParameter;

    /**
     * Cell membrane energy parameter.
     * For use in vertex-based models.
     */
    double mMembraneSurfaceEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter.
     * For use in vertex-based models.
     */
    double mCellCellAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter.
     * For use in vertex-based models.
     */
    double mCellBoundaryAdhesionEnergyParameter;

    /** Whether to write cell ID data to file. */
    bool mOutputCellIdData;

    /**
     * Whether to count the number of each cell mutation 
     * state and output to file.
     */
    bool mOutputCellMutationStates;

    /**
     * Whether to output the ancestor of each cell to a 
     * visualizer file.
     */
    bool mOutputCellAncestors;

    /**
     * Whether to count the number of each cell type and output 
     * to file.
     */
    bool mOutputCellTypes;

    /**
     * Whether to write the cell variables to a file.
     */
    bool mOutputCellVariables;

    /** Whether to write the cell cycle phases to a file. */
    bool mOutputCellCyclePhases;

    /** Whether to write the cell ages to a file. */
    bool mOutputCellAges;

	/** Whether to write the cell areas to a file. */
    bool mOutputCellAreas;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * As with other singleton classes, ensure the instance of this
     * class is serialized directly before being serialized via a
     * pointer.
     *
     * @param archive the archive
     * @param version the current version of this class
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
        archive & mMechanicsCutOffLength;
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
        archive & mMatureCellTargetArea;
        archive & mDeformationEnergyParameter;
        archive & mMembraneSurfaceEnergyParameter;
        archive & mCellCellAdhesionEnergyParameter;
        archive & mCellBoundaryAdhesionEnergyParameter;
        archive & mOutputCellIdData;
        archive & mOutputCellMutationStates;
        archive & mOutputCellAncestors;
        archive & mOutputCellTypes;
        archive & mOutputCellVariables;
        archive & mOutputCellCyclePhases;
        archive & mOutputCellAges;
        archive & mOutputCellAreas;
    }
};


#endif /*TISSUECONFIG_HPP_*/
