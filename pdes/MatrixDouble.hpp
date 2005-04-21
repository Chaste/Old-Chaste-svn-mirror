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
		double &MatrixDouble::operator()(int Row, int Column) const;
		static MatrixDouble Identity(int Size);
		int Rows( void ) const;
		int Columns( void ) const;
		MatrixDouble Inverse( void ) const;
		double Determinant( void ) const;
		
		MatrixDouble& operator*(double scalar);
		VectorDouble operator*(VectorDouble& rSomeVector);
};

#endif //_MATRIXDOUBLE_HPP_
