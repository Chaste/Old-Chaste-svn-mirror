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
		void ResetToZero( void );
		
		MatrixDouble& operator*(double scalar);
		friend MatrixDouble operator*(const double scalar, const MatrixDouble& rMatrix);
		friend MatrixDouble operator*(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix);
		friend MatrixDouble operator+(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix);
		friend MatrixDouble operator-(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix);
		
		friend VectorDouble operator*(const VectorDouble& rSomeVector, const MatrixDouble& rSomeMatrix);
		VectorDouble operator*(const VectorDouble& rSomeVector) const;
		MatrixDouble Transpose() const;
		bool IsSquare( void ) const;
		double GetTrace() const;
		double GetFirstInvariant() const;
		double GetSecondInvariant() const;
		double GetThirdInvariant() const;


};


#endif //_MATRIXDOUBLE_HPP_
