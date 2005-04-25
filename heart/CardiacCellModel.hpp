#ifndef _CARDIACCELLMODEL_HPP_
#define _CARDIACCELLMODEL_HPP_

/**
 * A basic cardiac cell model.
 * 
 * *** This class is currently redundant, but may be used in the future  ***
 */
class CardiacCellModel
{
	protected:
		double mTransmembranePotential;
		double mStimulusCurrent;
		// The total current crossing the membrane, Ii in:
		// dV/dt=-(1/C)(Ii+Istim)
		double mTotalIonicCurrent;
	
	public:
		// Constructor
		CardiacCellModel();
		// Overloaded Constructor
		CardiacCellModel(double transmembranePotential);
		// Destructor
		~CardiacCellModel();
		// Get methods
		double GetTransmembranePotential(void);
		double GetStimulusCurrent(void);
		double GetTotalIonicCurrent(void);
		// Set methods
		void SetStimulusCurrent(double magnitude);
};

#endif //_CARDIACCELLMODEL_HPP_
