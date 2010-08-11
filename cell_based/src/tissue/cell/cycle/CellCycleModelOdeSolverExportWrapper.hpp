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

#ifndef CELLCYCLEMODELODESOLVEREXPORTWRAPPER_HPP_
#define CELLCYCLEMODELODESOLVEREXPORTWRAPPER_HPP_

#include "CellCycleModelOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "HeunIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

#define EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(CLASS) \
    EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, CLASS, BackwardEulerIvpOdeSolver) \
    EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, CLASS, EulerIvpOdeSolver) \
    EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, CLASS, HeunIvpOdeSolver) \
    EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, CLASS, RungeKutta2IvpOdeSolver) \
    EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, CLASS, RungeKutta4IvpOdeSolver)

#endif /*CELLCYCLEMODELODESOLVEREXPORTWRAPPERHPP_*/
