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


#ifndef BIDOMAINCELLCOLLECTION_HPP_
#define BIDOMAINCELLCOLLECTION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

#include "AbstractCardiacCellCollection.hpp"
#include "AbstractConductivityTensors.hpp"




/**
 * #
 * 
 * 
 * BidomainCellCollection class.
 *
 * The bidmain equation is of the form:
 *
 * A_m ( C_m d(V_m)/dt + I_ion ) = div ( sigma_i grad( V_m + phi_e ) ) + I_intra_stim
 *   and
 * div ( (sigma_i + sigma_e) grad phi_e    +   sigma_i (grad V_m) )   = 0
 *
 * where V_m is the trans-membrane potential = phi_i - phi_e            (mV),
 *       phi_i is the intracellular potential                           (mV),
 *       phi_e is the intracellular potential                           (mV),
 * and   A_m is the surface area to volume ratio of the cell membrane   (1/cm),
 *       C_m is the membrane capacitance                                (uF/cm^2),
 *       sigma_i is the intracellular conductivity tensor               (mS/cm),
 *       sigma_e is the intracellular conductivity tensor               (mS/cm),
 * and   I_ion is the ionic current                                     (uA/cm^2),
 *       I_intra_stim is the internal stimulus                          (uA/cm^3),
 *
 * An extracellular stimulus can only be applied as a boundary condition through the Electrodes class, 
 * in which case the units are uA/cm^2.
 *
 */
template <unsigned SPACE_DIM>
class BidomainCellCollection : public virtual AbstractCardiacCellCollection<SPACE_DIM>
{
private:
    friend class TestBidomainCellCollection; // for testing.

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
        archive & boost::serialization::base_object<AbstractCardiacCellCollection<SPACE_DIM> >(*this);
        // Conductivity tensors are dealt with by HeartConfig, and the caches get regenerated.
    }

    /** Extracellular conductivity tensors. */
    AbstractConductivityTensors<SPACE_DIM> *mpExtracellularConductivityTensors;

    /**
     * Convenience method for extracellular conductivity tensors creation
     */
    void CreateExtracellularConductivityTensors();

public:
    /**
     * Constructor sets up extracellular conductivity tensors.
     * @param pCellFactory factory to pass on to the base class constructor
     */
    BidomainCellCollection(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory);

    /**
     *  Archiving constructor
     * @param rCellsDistributed  local cell models (recovered from archive)
     * @param pMesh  a pointer to the AbstractTetrahedral mesh (recovered from archive).
     */
    BidomainCellCollection(std::vector<AbstractCardiacCell*> & rCellsDistributed,AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh);

    /**
     * Destructor
     */
    ~BidomainCellCollection();

    /**
     * Get the extracellular conductivity tensor for the given element
     * @param elementIndex  index of the element of interest
     */
     const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetExtracellularConductivityTensor(unsigned elementIndex);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BidomainCellCollection)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const BidomainCellCollection<SPACE_DIM> * t, const unsigned int file_version)
{
    const AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* p_mesh = t->pGetMesh();
    ar & p_mesh;

    // Don't use the std::vector serialization for cardiac cells, so that we can load them
    // more cleverly when migrating checkpoints.
    t->SaveCardiacCells(ar, file_version);

    // CreateIntracellularConductivityTensor() is called by constructor and uses HeartConfig. So make sure that it is
    // archived too (needs doing before construction so appears here instead of usual archive location).
    HeartConfig* p_config = HeartConfig::Instance();
    ar & *p_config;
    ar & p_config;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, BidomainCellCollection<SPACE_DIM> * t, const unsigned int file_version)
{
    std::vector<AbstractCardiacCell*> cells_distributed;
    AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* p_mesh;

    ar & p_mesh;
    // Load only the cells we actually own
    AbstractCardiacCellCollection<SPACE_DIM,SPACE_DIM>::LoadCardiacCells(
            *ProcessSpecificArchive<Archive>::Get(), file_version, cells_distributed, p_mesh);

    // CreateIntracellularConductivityTensor() is called by AbstractCardiacCellCollection constructor and uses HeartConfig.
    // (as does CreateExtracellularConductivityTensor). So make sure that it is
    // archived too (needs doing before construction so appears here instead of usual archive location).
    HeartConfig* p_config = HeartConfig::Instance();
    ar & *p_config;
    ar & p_config;

    ::new(t)BidomainCellCollection<SPACE_DIM>(cells_distributed, p_mesh);
}
}
} // namespace ...


#endif /*BIDOMAINCELLCOLLECTION_HPP_*/
