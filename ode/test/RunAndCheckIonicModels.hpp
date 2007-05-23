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

void SaveSolution(std::string baseResultsFilename, AbstractCardiacCell *pOdeSystem,
                  OdeSolution& rSolution, int stepPerRow);
                  
void CheckCellModelResults(std::string baseResultsFilename);

void CompareCellModelResults(std::string baseResultsFilename1, std::string baseResultsFilename2, double tolerance);


#endif //_RUNANDCHECKIONICMODELS_HPP_
