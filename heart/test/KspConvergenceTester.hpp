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

#ifndef KSPCONVERGENCETESTER_HPP_
#define KSPCONVERGENCETESTER_HPP_

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class KspConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM, PROBLEM_DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->SetKspRelativeTolerance(1e-2);
    }
    void UpdateConvergenceParameters()
    {
        this->SetKspRelativeTolerance(this->GetKspRelativeTolerance()*0.1);
    
    }
    bool GiveUpConvergence()
    {
        return this->GetKspRelativeTolerance()<1e-9;
    }
    double Abscissa()
    {
        return this->GetKspRelativeTolerance();
    }
};

#endif /*KSPCONVERGENCETESTER_HPP_*/
