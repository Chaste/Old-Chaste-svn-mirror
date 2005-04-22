#ifndef _TESTLR91_HPP_
#define _TESTLR91_HPP_


#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "SodiumCurrentLR91.hpp"
#include "BackgroundCurrentLR91.hpp"
#include "PlateauPotassiumCurrentLR91.hpp"
#include "CalciumConcentrationLR91.hpp"
#include "SlowInwardCurrentLR91.hpp"
#include "PotassiumTimeDependentCurrentLR91.hpp"
#include "TransmembranePotentialLR91.hpp"
//#include "TimeDependentPotassiumCurrent.hpp"
//#include "LuoRudyModel.hpp" 
#include "ConstantsLR91.hpp"

class TestLR91 : public CxxTest::TestSuite
{
    public:
    
    // Tests that currents INa and IK are calculated correctly
    void donttestCurrentsMagnitude( void )
    {
//        SodiumCurrentLR91 *pINa;
//        //SodiumCurrent(double m, double h);
//        double m = 0.5;
//        double h = 0.01;
//        double j = 0.1;
//        double V = -10.0;
//        
//        pINa = new SodiumCurrentLR91();
//        
//        pINa->UpdateMagnitudeOfCurrent(V, m, h, j);
//        double iNa = pINa->GetMagnitudeOfCurrent();
//        
//        double mPrime = pINa->ComputeMPrime(V, m, h, j);
//        double hPrime = pINa-> ComputeHPrime(V, m, h, j);
//        double jPrime = pINa->ComputeJPrime(V, m, h, j);
//        
//        
//        std::cout << "\n";
//        std::cout << "INa is equal to " << iNa << "\n";
//        TS_ASSERT_DELTA( iNa, gNa * pow(m,3) * h * j * (V- eNa),0.0001);
//        std::cout << "\n";
//        std::cout << "mPrime is equal to " << mPrime << "\n";
//        std::cout << "hPrime is equal to " << hPrime<< "\n";
//        std::cout << "jPrime is equal to " << jPrime << "\n";   
//    
//    // Tests that currents Ib and IKp are calculated correctly
//        V = -10.0;
//        BackgroundCurrentLR91 *pIB;
//        pIB = new BackgroundCurrentLR91();
//        pIB->UpdateMagnitudeOfCurrent(V);
//        double iB = pIB->GetMagnitudeOfCurrent();
//        
//        std::cout << "\n"  <<"Ib is equal to " << iB << "\n";
//        TS_ASSERT_DELTA( iB, gB * (V + eB),0.0001);
//        
//        PlateauPotassiumCurrentLR91 * pIKp;
//        pIKp = new PlateauPotassiumCurrentLR91();
//        pIKp->UpdateMagnitudeOfCurrent(V);
//        double iKp = pIKp->GetMagnitudeOfCurrent();
//        
//        std::cout << "\n" << "IKp is equal to " <<  iKp << "\n";
//        TS_ASSERT_DELTA(  iKp , gKp *  1/ ( 1 + exp(7.488 - V)/5.98) * (V + eKp),0.0001);
//        
// 
//   // Tests that currents ICa compiles
//        double caI = 4; 
//        CalciumConcentrationLR91 *pCaI;
//        pCaI = new CalciumConcentrationLR91();
//        pCaI->SetMagnitudeOfIonicConcentration(caI);
//        double iCa = pCaI->GetMagnitudeOfIonicConcentration();
//         
//        std::cout << "\n" << "Intracellular calcium is " <<  iCa << "\n";      
//        
//    // test that Isi computes correctly
//        double d = 0.3;
//        double f = 0.6;
//        SlowInwardCurrentLR91 *pISi;
//        pISi = new SlowInwardCurrentLR91();
//        pISi->UpdateMagnitudeOfCurrent(V,d,f,caI);
//        double iSi = pISi->GetMagnitudeOfCurrent();
//        double eSi = 7.7 - 13.0287 * log(caI);
//        
//        std::cout << "\n" << "Slow inward current is " <<  iSi << "\n";  
//        TS_ASSERT_DELTA(  iSi , gSi * d *f * (V - eSi), 0.0001);    
//        
//         
//        
//        //test transmembrane potentail works;
//        double iTotal = 10;
//        double iStim = 5;
//        TransmembranePotentialLR91 *pV;
//        pV = new TransmembranePotentialLR91();
//        double VPrime = pV->ComputeVPrime(iStim, iTotal);
//        
//        std::cout << "\n" << "Transmembrane potential is " <<  VPrime << "\n";  
        
    }
    
};

#endif //_TESTLR91_HPP_
