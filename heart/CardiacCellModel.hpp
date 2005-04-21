#ifndef _CARDIACCELLMODEL_HPP_
#define _CARDIACCELLMODEL_HPP_

class CardiacCellModel
{
	protected:
		double mTransmembranePotential;
		double mStimulusCurrent;
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
