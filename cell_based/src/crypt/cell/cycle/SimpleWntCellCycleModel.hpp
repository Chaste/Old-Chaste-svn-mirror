/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef SIMPLEWNTCELLCYCLEMODEL_HPP_
#define SIMPLEWNTCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "WntConcentration.hpp"


/**
 *  Simple Wnt-dependent cell cycle model
 */
class SimpleWntCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;

        archive & mUseCellProliferativeTypeDependentG1Duration;
        archive & mDimension;
    }

    /**
     * Whether to use different mean G1 durations for different cell types.
     * For use in SetG1Duration().
     */
    bool mUseCellProliferativeTypeDependentG1Duration;

    /**
     * The spatial dimension (needed by the templated class WntConcentration).
     */
    unsigned mDimension;

    /**
     * Get the Wnt level experienced by the cell.
     */
    double GetWntLevel();

    /**
     * Get the type of Wnt concentration (LINEAR, RADIAL, or NONE).
     * This affects how the cell cycle phase is updated.
     */
    WntConcentrationType GetWntType();

protected:

    /**
     * Stochastically set the G1 duration. The G1 duration is taken
     * from a normal distribution whose mean is the G1 duration given
     * in TissueConfig for the cell type and whose standard deviation
     * is 1.
     *
     * Called on cell creation at the start of a simulation, and for both
     * parent and daughter cells at cell division.
     */
    void SetG1Duration();

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called.
     *
     * @param dimension the spatial dimension (needed by the templated class WntConcentration)
     * @param useCellProliferativeTypeDependentG1Duration  Whether the duration of the G1 phase is dependent on cell type
     */
    SimpleWntCellCycleModel(unsigned dimension, bool useCellProliferativeTypeDependentG1Duration=false);

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    void UpdateCellCyclePhase();

    /**
     * Overridden ResetForDivision() method.
     */
    void ResetForDivision();

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    void InitialiseDaughterCell();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Get the spatial dimension.
     *
     * @return mDimension
     */
    unsigned GetDimension();
};

#include "TemplatedExport.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SimpleWntCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SimpleWntCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const SimpleWntCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SimpleWntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SimpleWntCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of SimpleWntCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    unsigned dimension = UINT_MAX;
    ::new(t)SimpleWntCellCycleModel(dimension);
}
}
} // namespace ...

#endif /*SIMPLEWNTCELLCYCLEMODEL_HPP_*/
