#ifndef _CONSTANTSLR91_HPP_
#define _CONSTANTSLR91_HPP_

// Sodium Dependent Constants
const double eNa = 54.4;        // milivolts
const double gNa = 23.0;          // mS/cm^2

// Slow Inward Current
const double gSi = 0.09;        // mS/cm^2

// Potassium Dependent Constants
const double eK = -77.0;          // millivolts
const double gK = 0.282;        // mS/cm^2
const double eK1 = eK;
const double gK1 = 0.6047;       // mS/cm^2
const double gKp = 0.0183;      // mS/cm^2
const double eKp = eK;

// Background Currents
const double gB = 0.03921;      // mS/cm^2
const double eB = 59.87;       // millivolts

// Other Constants
const double membraneCapacitance = 1.0; // microFarads/ mm^-2

 
  



#endif //_CONSTANTSLR91_HPP_
