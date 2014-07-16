/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef ACTIVATIONOUTPUTMODIFIER_HPP_
#define ACTIVATIONOUTPUTMODIFIER_HPP_

#include "AbstractOutputModifier.hpp"

/**
 * Specialised class for on-the-fly calculation of activation and recovery times.
 * This is hard-coded for the first 2 rounds of activation/recovery.
 *
 * Times are relative to simulation start and are listed in node order and in 4 columns comma separated (.CSV)
 *
 * node_0_activation, node_0_recovery, node_0_later_activation, node_0_later_recovery
 * node_1_activation, node_1_recovery, node_1_later_activation, node_1_later_recovery
 * ...
 *
 * Any missing data (node was not activated a second time or did not recover) is marked with -1.
 */
class ActivationOutputModifier : public AbstractOutputModifier
{
private:
    double mThreshold; /**< The user-defined threshold at which activation and recovery is to be measured */
    unsigned mLocalSize; /**< The number of nodes on this process (calculated in #InitialiseAtStart)*/
    std::vector<double> mFirstActivitationTimes; /**< The first activation (first time above threshold) for all local nodes on this process*/
    std::vector<double> mFirstRecoveryTimes; /**< The first recovery (first time subsequent time below threshold) for all local nodes on this process*/
    std::vector<double> mSecondActivitationTimes; /**< The second activation time for local nodes */
    std::vector<double> mSecondRecoveryTimes; /**< The second recovery time for local nodes */
public:
    /**
     * Constructor
     *
     * @param rFilename  The file which is eventually produced by this modifier
     * @param threshold  The transmembrane voltage threshold (in mV) at which activation is deemed to have been trigged.
     *                   This is also used for calculating relaxation time (this is not so sophisticated as an APD90 calculation).
     *
     */
    ActivationOutputModifier(const std::string& rFilename, double threshold)
        : AbstractOutputModifier(rFilename),
          mThreshold(threshold)
    {
    }
    /**
     * Initialise the modifier (make space for local activation times) when the solve loop is starting.
     *
     * Note the problem passes parameters in a non-templated fashion in order to keep the interface as lightweight as
     * possible.
     * @param pVectorFactory  The vector factory which is associated with the calling problem's mesh
     */
    virtual void InitialiseAtStart(DistributedVectorFactory* pVectorFactory);

    /**
     * Finalise the modifier (write all results to the file)
     */
    virtual void FinaliseAtEnd();

    /**
     * Process a solution time-step (memoize all new activations)
     * @param time  The current simulation time
     * @param solution  A working copy of the solution at the current time-step.  This is the PETSc vector which is distributed across the processes.
     * @param problemDim  The calling problem dimension. Used here to avoid probing the size of the solution vector
     */
    virtual void ProcessSolutionAtTimeStep(double time, Vec solution, unsigned problemDim);
};

#endif /* ACTIVATIONOUTPUTMODIFIER_HPP_ */
