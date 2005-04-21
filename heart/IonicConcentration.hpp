#ifndef _IONICCONCENTRATION_HPP_
#define _IONICCONCENTRATION_HPP_

class IonicConcentration
{
	protected:
		double mMagnitudeOfConcentration;
	
	public:
		// Constructor
		IonicConcentration(void);
		// Overloaded Constructor
		IonicConcentration(double transmembranePotential);
		// Destructor
		~IonicConcentration();
		// Get methods
		double GetConcentration(void);
};

#endif //_IONICCONCENTRATION_HPP_
