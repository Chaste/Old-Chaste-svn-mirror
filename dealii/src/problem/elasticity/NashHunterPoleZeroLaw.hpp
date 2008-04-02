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

#ifndef NASHHUNTERPOLEZEROLAW_HPP_
#define NASHHUNTERPOLEZEROLAW_HPP_

#include "PoleZeroMaterialLaw.hpp"

/**
 *  The pole-zero law but with parameters set using values from Chapter 41 of [??], 
 *  Remme, Nash and Hunter
 * 
 *  The stiffness are in KILOPASCALS
 */
template<unsigned DIM>
class NashHunterPoleZeroLaw : public PoleZeroMaterialLaw<DIM>
{
friend class TestMaterialLaws;

public :
    NashHunterPoleZeroLaw()
        : PoleZeroMaterialLaw<DIM>()
    {
        assert(DIM==2 || DIM==3);
        if (DIM==2)
        {
            SetUpFor2d();
        }
        else
        {
            SetUpFor3d();
        }
    }
    
    void SetUpFor3d()
    {        
        assert(DIM==3);
        
        std::vector<std::vector<double> > k(3),a(3),b(3);
        for(unsigned i=0; i<3; i++)
        {
            k[i].resize(3);
            a[i].resize(3);
            b[i].resize(3);
        }
        
        /////////////////////////////////////////////////////////////////
        // Everything here has been entered in kPa!!
        // Currently the NHS cellular model returns the active tension
        // in kPa so no scaling is needed
        /////////////////////////////////////////////////////////////////
        k[0][0] = 2; //ff
        k[1][0] = k[0][1] = 1; //fs 
        k[0][2] = k[2][0] = 1; //fn
        k[1][1] = 2; //ss
        k[1][2] = k[2][1] = 1; //sn 
        k[2][2] = 2; //nn
        
        // dimensionless
        a[0][0] = 0.475; //ff
        a[1][0] = a[0][1] = 0.8; //fs 
        a[2][0] = a[0][2] = 0.8; //fn
        a[1][1] = 0.619; //ss
        a[2][1] = a[1][2] = 0.8; //sn 
        a[2][2] = 0.943; //nn
        
        // dimensionless
        b[0][0] = 1.5;
        b[1][0] = b[0][1] = 1.2; 
        b[2][0] = b[0][2] = 1.2; 
        b[1][1] = 1.5;
        b[2][1] = b[1][2] = 1.2; 
        b[2][2] = 0.442;
        
        this->SetParameters(k,a,b);
    }    

    void SetUpFor2d()
    {        
        assert(DIM==2);
        
        std::vector<std::vector<double> > k(2),a(2),b(2);
        for(unsigned i=0; i<2; i++)
        {
            k[i].resize(2);
            a[i].resize(2);
            b[i].resize(2);
        }
        
        k[0][0] = 2; //ff
        k[1][0] = k[0][1] = 1; //fs 
        k[1][1] = 2; //ss
        
        a[0][0] = 0.475; //ff
        a[1][0] = a[0][1] = 0.8; //fs 
        a[1][1] = 0.619; //ss
        
        b[0][0] = 1.5;
        b[1][0] = b[0][1] = 1.2; 
        b[1][1] = 1.5;
        
        this->SetParameters(k,a,b);
    }    
};

#endif /*NASHHUNTERPOLEZEROLAW_HPP_*/
