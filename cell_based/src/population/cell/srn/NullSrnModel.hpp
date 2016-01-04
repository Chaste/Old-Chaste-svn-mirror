/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef NULLSRN_HPP_
#define NULLSRN_HPP_

#include <vector>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractSrnModel.hpp"
#include "SimulationTime.hpp"

class AbstractSrnModel;

/**
 * This class contains a dummy/null SRN model that can be used for any cell models
 * that do not have an SRN, or where the SRN is combined with the main cell-cycle
 * model.
 */

class NullSrnModel : public AbstractSrnModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the srn model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSrnModel>(*this);
    }

protected:

public:
    /**
     * Creates an NullSrnModel, calls SetBirthTime on the
     * AbstractSrnModel to make sure that can be set 'back in time' for
     * cells which did not divide at the current time.
     *
     */
    NullSrnModel();

    /**
     * Destructor.
     */
    virtual ~NullSrnModel();

    /*
     * Overridden  SimulateToCurrentTime method. Here we dont do anything as this is a
     * Null model.
     */
    void SimulateToCurrentTime();

    /*
     * Overridden CreateSrnModel method.
     *
     * Builder method to create new instances of the SRN model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     *
     * This method is called by Cell::Divide() to create a SRN
     * model for the daughter cell.  Note that the parent SRN
     * model will have had ResetForDivision() called just before
     * CreateSrnModel() is called, so performing an exact copy of the
     * parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @return new SRN model
     */
    AbstractSrnModel* CreateSrnModel();
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(NullSrnModel)

#endif /* NULLSRN_HPP_ */
