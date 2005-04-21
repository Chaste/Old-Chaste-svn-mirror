#ifndef _VECTORDOUBLE_HPP_
#define _VECTORDOUBLE_HPP_

class VectorDouble
{
	private:
		int mSize;
		double *mElementArray;
	public:
		VectorDouble(int Size);
		VectorDouble(const VectorDouble& rOtherVector);
		~VectorDouble();
		VectorDouble& operator=(const VectorDouble& rOtherVector);
		double &VectorDouble::operator()(int Entry) const;
		int Size( void ) const;
		
};

#endif //_VECTORDOUBLE_HPP_
