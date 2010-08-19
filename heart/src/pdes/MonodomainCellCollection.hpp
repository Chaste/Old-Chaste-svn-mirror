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


#ifndef MONODOMAINCELLCOLLECTION_HPP_
#define MONODOMAINCELLCOLLECTION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include <vector>
#include "AbstractCardiacCellCollection.hpp"


/**
 * MonodomainCellCollection class.
 *
 * The monodomain equation is of the form:
 * A (C dV/dt + Iionic) +Istim = Div( sigma_i Grad(V) )
 *
 * where A is the surface area to volume ratio         (1/cm),
 *       C is the capacitance                          (uF/cm^2),
 *       sigma_i is the intracellular conductivity     (mS/cm),
 *       I_ionic is the ionic current                  (uA/cm^2),
 *       I_stim is the intracellular stimulus current  (uA/cm^3).
 *
 * Note that default values of A, C and sigma_i are stored in the parent class
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class MonodomainCellCollection : public virtual AbstractCardiacCellCollection<ELEMENT_DIM,SPACE_DIM>
{
private:
    friend class TestMonodomainCellCollection;

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
        archive & boost::serialization::base_object<AbstractCardiacCellCollection<ELEMENT_DIM, SPACE_DIM> >(*this);
    }


public:
    /**
     *  Constructor
     *
     *  @param pCellFactory  Provides the mesh and cells
     */
    MonodomainCellCollection(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory);


    /**
     * Another constructor (for archiving)
     *
     * @param rCellsDistributed  local cell models (recovered from archive)
     * @param pMesh the mesh (also recovered from archive)
     */
    MonodomainCellCollection(std::vector<AbstractCardiacCell*> & rCellsDistributed,
                  AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);


};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp" // Must be last
EXPORT_TEMPLATE_CLASS2(MonodomainCellCollection, 1, 1)
EXPORT_TEMPLATE_CLASS2(MonodomainCellCollection, 1, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainCellCollection, 1, 3)
EXPORT_TEMPLATE_CLASS2(MonodomainCellCollection, 2, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainCellCollection, 3, 3)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const MonodomainCellCollection<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    const AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh = t->pGetMesh();
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
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, MonodomainCellCollection<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    std::vector<AbstractCardiacCell*> cells_distributed;
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh;

    ar & p_mesh;
    // Load only the cells we actually own
    AbstractCardiacCellCollection<ELEMENT_DIM,SPACE_DIM>::LoadCardiacCells(
            *ProcessSpecificArchive<Archive>::Get(), file_version, cells_distributed, p_mesh);

    // CreateIntracellularConductivityTensor() is called by AbstractCardiacCellCollection constructor and uses HeartConfig.
    // So make sure that it is archived too (needs doing before construction so appears here instead of usual archive location).
    HeartConfig* p_config = HeartConfig::Instance();
    ar & *p_config;
    ar & p_config;

    ::new(t)MonodomainCellCollection<ELEMENT_DIM, SPACE_DIM>(cells_distributed, p_mesh);
}
}
} // namespace ...


#endif /*MONODOMAINCELLCOLLECTION_HPP_*/
