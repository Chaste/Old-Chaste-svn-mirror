#ifndef _IONICCONCENTRATION_HPP_
#define _IONICCONCENTRATION_HPP_
/**
 * An ionic concentration. Base class for specific ionic concentrations. 
 * In Lou Rudy model only the intracellular calcium concentration varies
 */
class IonicConcentration
{
	protected:
		double mMagnitudeOfIonicConcentration;
	
	public:
		// Constructor
		IonicConcentration(void);
		// Overloaded Constructor
		IonicConcentration(double transmembranePotential);
		// Destructor
		~IonicConcentration();
		// Get & Set methods
		double GetMagnitudeOfIonicConcentration(void);
		void SetMagnitudeOfIonicConcentration(const double &rIonicConcentration);
};

#endif //_IONICCONCENTRATION_HPP_
