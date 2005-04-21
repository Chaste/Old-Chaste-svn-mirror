#ifndef _MATRIXDOUBLE_HPP_
#define _MATRIXDOUBLE_HPP_
#include "VectorDouble.hpp"

class MatrixDouble
{
	private:
		int mRows;
		int mColumns;
		int mNumberOfElements;
		double *mElementArray;
	public:
		MatrixDouble(int Rows, int Columns);
		MatrixDouble(const MatrixDouble& rOtherMatrix);
		~MatrixDouble();
		MatrixDouble& operator=(const MatrixDouble& rOtherMatrix);
		double &MatrixDouble::operator()(int Row, int Column);
		static MatrixDouble Identity(int Size);
		int Rows( void );
		int Columns( void );
		MatrixDouble Inverse( void );
		double Determinant( void );
		
		MatrixDouble& operator*(double scalar);
		VectorDouble operator*(VectorDouble& rSomeVector);
};

#endif //_MATRIXDOUBLE_HPP_
