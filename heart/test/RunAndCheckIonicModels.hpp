/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _RUNANDCHECKIONICMODELS_HPP_
#define _RUNANDCHECKIONICMODELS_HPP_

#include <vector>
#include <string>

#include "OdeSolution.hpp"

#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"

#include "AbstractCardiacCell.hpp"

void RunOdeSolverWithIonicModel(AbstractCardiacCell *pOdeSystem,
                                double endTime,
                                std::string filename,
                                int stepPerRow=100,
                                bool doComputeExceptVoltage=true);
                  
void CheckCellModelResults(std::string baseResultsFilename);

void CompareCellModelResults(std::string baseResultsFilename1, std::string baseResultsFilename2, double tolerance);


#endif //_RUNANDCHECKIONICMODELS_HPP_
