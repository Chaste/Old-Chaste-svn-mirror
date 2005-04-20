#ifndef _IONICCURRENT_HPP_
#define _IONICCURRENT_HPP_

class IonicCurrent
{
	protected:
		double mMagnitudeOfCurrent;
	
	public:
		// Constructor
		IonicCurrent(void);
		IonicCurrent(const double &rMagnitudeOfCurrent);
		// Destructor
		~IonicCurrent(void);
		// Set and get methods
		void SetMagnitudeOfCurrent(const double &rMagnitudeOfCurrent);
		double GetMagnitudeOfCurrent(void);
};

#endif //_IONICCURRENT_HPP_
