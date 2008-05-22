/*

Copyright (C) University of Oxford, 2008

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
#ifndef OXYGENBASEDCELLKILLER_HPP_
#define OXYGENBASEDCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/** 
 *  Kills cells that have experienced a prolonged continuous period of hypoxia.
 * 
 *  The non-dimensionalised oxygen concentration at which cells become 
 *  hypoxic is optionally passed into the constructor. 
 */

template <unsigned SPACE_DIM>
class OxygenBasedCellKiller : public AbstractCellKiller<SPACE_DIM>
{
private:

    // The oxygen concentration below which cells become hypoxic
    double mHypoxicConcentration;
     
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<SPACE_DIM> >(*this);        
    }

public:

    /**
     * Constructor.
     * 
     * @param pointer to the tissue.
     * @param oxygen concentration below which cells become hypoxic.
     */  
    OxygenBasedCellKiller(AbstractTissue<SPACE_DIM>* pTissue, 
                          double concentration=CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration());
    
    /**
     * Set method for mHypoxicConcentration.
     * 
     * @param hypoxicConcentration the oxygen concentration below which cells become hypoxic.
     */
    void SetHypoxicConcentration(double hypoxicConcentration);
    
    /**
     * Get method for mHypoxicConcentration.
     */
    double GetHypoxicConcentration() const;
    
    /**
     *  Starts apoptosis if the cell has has been hypoxic for longer than 
     *  some critical period, and  it is currently hypoxic, and a random number 
     *  is less than some probability of death (which scales linearly with the 
     *  local oxygen concentration).
     * 
     *  @param rCell reference to the cell to test for apoptosis.
     */  
    void TestAndLabelSingleCellForApoptosis(TissueCell& rCell);
    
    /**
     * Loop over cells and start apoptosis if the cell has been undergone
     * a prolonged period of hypoxia.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath();

};

template <unsigned SPACE_DIM>
OxygenBasedCellKiller<SPACE_DIM>::OxygenBasedCellKiller(AbstractTissue<SPACE_DIM>* pTissue, 
                                                        double concentration)
        : AbstractCellKiller<SPACE_DIM>(pTissue),
          mHypoxicConcentration(concentration) 
{
}
          
template <unsigned SPACE_DIM>
void OxygenBasedCellKiller<SPACE_DIM>::SetHypoxicConcentration(double hypoxicConcentration)
{
    mHypoxicConcentration = hypoxicConcentration;    
}

template <unsigned SPACE_DIM>
double OxygenBasedCellKiller<SPACE_DIM>::GetHypoxicConcentration() const
{
    return mHypoxicConcentration;
}

template <unsigned SPACE_DIM>
void OxygenBasedCellKiller<SPACE_DIM>::TestAndLabelSingleCellForApoptosis(TissueCell& rCell)
{   
    if (rCell.GetCellType()==NECROTIC && !(rCell.HasApoptosisBegun()))
    {
        rCell.StartApoptosis();
    }          
}

template <unsigned SPACE_DIM>
void OxygenBasedCellKiller<SPACE_DIM>::TestAndLabelCellsForApoptosisOrDeath()
{      
    for (typename AbstractTissue<SPACE_DIM>::Iterator cell_iter = this->mpTissue->Begin();
        cell_iter != this->mpTissue->End();
        ++cell_iter)
    {               
        TestAndLabelSingleCellForApoptosis(*cell_iter);
    }
}     
    
#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(OxygenBasedCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OxygenBasedCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const AbstractTissue<DIM>* const p_tissue = t->GetTissue();
    ar << p_tissue;
    double conc = t->GetHypoxicConcentration();
    ar << conc;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OxygenBasedCellKiller<DIM> * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    double conc;
    ar >> conc;
    // invoke inplace constructor to initialize instance
    ::new(t)OxygenBasedCellKiller<DIM>(p_tissue, conc);
}
}
} // namespace ...

#endif /*OXYGENBASEDCELLKILLER_HPP_*/
