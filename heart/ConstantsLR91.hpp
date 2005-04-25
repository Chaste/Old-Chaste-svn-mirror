#ifndef _CONSTANTSLR91_HPP_
#define _CONSTANTSLR91_HPP_

// Sodium Current Constants
const double eNa = 54.4;          // milivolts
const double gNa = 23.0;          // mS/cm^2

// Slow Inward Current Constant
const double gSi = 0.09;          // mS/cm^2

// Potassium Current Constants
const double eK = -77.0;          // millivolts
const double gK = 0.282;          // mS/cm^2
const double eK1 = eK;
const double gK1 = 0.6047;        // mS/cm^2
const double eKp = eK;            // millivolts
const double gKp = 0.0183;        // mS/cm^2

// Background Current Constants
const double eB = 59.87;          // millivolts
const double gB = 0.03921;        // mS/cm^2
// -59.87 in original paper, formula for Ib changed from:
//
// mMagnitudeOfCurrent = gB * (voltage - eB)
//
// to:
//
// mMagnitudeOfCurrent = gB * (voltage + eB), so actually equivalent

// Other Constants
const double membraneCapacitance = 1.0; // microFarads/ mm^-2

#endif //_CONSTANTSLR91_HPP_
