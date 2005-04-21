#ifndef _TINYMATRIX_HPP_
#define _TINYMATRIX_HPP_

class TinyMatrix
{
private:
		int mRows;
		int mColumns;
		int mNumberOfElements;
		double *mElementArray;
	public:
		TinyMatrix(int Rows, int Columns);
		TinyMatrix(const TinyMatrix& rOtherMatrix);
		~TinyMatrix();
		TinyMatrix& operator=(const TinyMatrix& rOtherMatrix);
		double &TinyMatrix::operator()(int Row, int Column);
		static TinyMatrix Identity(int Size);
		double Determinant( void );
		TinyMatrix Inverse( void );
		int Rows( void );
		int Columns( void );		
};

#endif //_TINYMATRIX_HPP_
