/* Copyright dyj (c) 2018
 * any codes cannot be used for business
 * template for project special
 * updated on 2018/1/2
 * Translation:
 *	矩阵			matrix
 *	向量			vector
 *	行				row
 *	列				column
 *	秩				rank
 *
 *	方阵			square
 *	单位矩阵		identity
 *
 *	行列式			determinant
 *	逆				inversion
 *	求逆			inverse
 *	转置矩阵		transpose
 *	伴随矩阵		adjoint
 *
 *	行主元			row_pivot
 *	列主元			column_pivot
 *	行阶梯形		row_echelon
 *	列阶梯形		column_echelon
 *	行最简形		row_canonical
 *	列最简形		column_canonical
 *	化为行阶梯形	row_eliminate
 *	化为列阶梯形	column_eliminate
 *	化为行最简形	row_reduce
 *	化为列最简形	column_reduce
 *
 * Introductions:
 *	class la_vector:
 *		private members:
 *			_vec		where the values are stored
 *			_dimention	storing the dimention of the vector
 *			is_row		false	column vector
 *						true	row vector
 *		public members:
 *			ROW,COLUMN	constant flag for constructing vectors
 *			transpose()	to transpose this function and return its reference
 *			row()		to return if the vector is a row vector
 *			dimention()	to return the dimention of the vector
 *
 *	class matrix:
 *		private members:
 *			step_by_step	flag to sign if print the steps during solving, 
 *							which will be set to false defaultly
 *			_row,_column	storing the number of rows/columns of the matrix
 *			_mat			storing the matrix itself
 *		public members:
 *			transposition()	to return the copy of it's transposition
 *			inversion()	to return the copy of it's inversion
 *			column()
 *			row()		to return the number of columns/rows of the matrix
 *			identity()	to check if this matrix is an identity matrix
 *			square()	to check if this matrix is a square matrix
 *			rank()				to calculate the rank of the matrix
 *			print(_os = cout)	to print the matrix to _os, in a formated way
 *			set_step_flags(_f)	to set the flags of step_by_step. Returning 1 if
 *								the flag is changed; 0 if unchanged
 *			block(_u, _d, _l, _r)
 *						returning the block matrix from [_u][_l] to [_d][_r]
 *
 *			row_pivot(_r)		to find the first non-zero element in row _r, 
 *								returning -1 if not found
 *			row_echelon()		to check if the matrix is row_echelon_form matrix
 *			row_exchange(_i, _j)		to exchange the row of _i and _j
 *			row_time(_i, _k)			to make each element on row _i time _k
 *			row_time_plus(_i, _k, _j)	to add _k times of row _i to row _j
 *			column_pivot(_r)
 *			column_echelon()
 *			column_exchange(_i, _j)
 *			column_time(_i, _k)
 *			column_time_plus(_i, _k, _j)
 *
 *			operator []	to return the certian row in the matrix, range: [1, _row]
 *						WARNING never access the _mat directly, that will make
 *						the program insanely chaotic, for the index of this []
 *						is actually different from that of _mat
 *		functions related:
 *			TODO
 *			solve
 *			equivallent
 *			determinant
 *			resize
 *
 *			row_eliminate(_mat, _row_start, _column_start)
 *				to make _mat a echelon form row equivalent to the original
 *				one. Directly modify on _mat
 *			row_reduce(_mat)
 *				to make _mat to a reduced echelon form
 *			inverse(_mat)
 *				to return the inversion of _mat, if exists;
 *				empty matrix, if non-exists;
 *			transpose(_mat)
 *				to return the transposition of _mat
 *			concat(_mat1, _mat2)
 *				to attach _mat2 to _mat1 horizontally
 *			vconcat(_mat1, _mat2)
 *				to attach _mat2 to _mat1 vertically
 *			equivallent(_mat1, _mat2)
 *				to check if _mat1 and _mat2 are equivallent
 */
#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cassert>
#include "rational.h"

using namespace std;
using namespace rat;

//linear algebra namespace
namespace lnr {

	//key class: linear algebra vector
	//typedef rat::rational T;
	typedef rat::rational T;
	class la_vector {
		private:
			vector<T> _vec;	//vector<int> to store the value
			int _dimention;		//dimention of the vector
			bool is_row;		//flag of column/row
		public:
			//flags for direction
			const bool ROW = true, COLUMN = false;
			//constructors
			la_vector();
			la_vector(const la_vector&) = default;
			explicit la_vector(const int &_dimention);
			explicit la_vector(const int &_dimention, const bool &_direction);
			la_vector(const int &_dimention, T *_array);
			la_vector(const int &_dimention, const bool &_direction, T *_array);

			//operator for assignment
			la_vector &operator = (const la_vector&);
			la_vector &operator += (const la_vector&);
			la_vector &operator -= (const la_vector&);
			la_vector &operator *= (const T&);	//times a number
			T &operator [] (const int&);
			const T &operator [] (const int&) const;	//this version is for const ones

			//operators
			friend la_vector operator + (const la_vector&, const la_vector&);
			friend la_vector operator - (const la_vector&, const la_vector&);
			friend la_vector operator + (const la_vector&);
			friend la_vector operator - (const la_vector&);
			friend la_vector operator * (const T&, const la_vector&);	//dot multiplication
			friend la_vector operator * (const la_vector&, const T&);	//dot multiplication
			friend T operator * (const la_vector&, const la_vector&);	//times a number

			//functions
			//la_vector transpose() const;	//return the copy of the transformred vector
			la_vector &transpose();		//transpose this vector
			bool row() const;
			int dimention() const;
	};

	class matrix_vector {
		private:
			vector<T> _vec;
		public:
			matrix_vector() = default;
			
			vector<T> &vec() {
				return _vec;
			}
			const vector<T> &vec() const{
				return _vec;
			}
			T &operator [] (const int &i) {
				return _vec[i-1];
			}
			const T &operator [] (const int &i) const {
				return _vec[i-1];
			}
	};

	class matrix{
		private:
			//static member flag
			static bool step_by_step;

			int _row, _column;
			vector<matrix_vector> _mat;

		public:

			static bool vertical, horizontal;

			//constructors
			matrix();
			explicit matrix(int n);	//elemental TEST
			matrix(const matrix &);
			matrix(const matrix &_matl, const matrix &_matr, bool _ver = false);
			matrix(const la_vector &);
			matrix(const int &_row, const int &_column);
			matrix(const int &_row, const int &_column, T **_array);

			//assignments
			matrix &operator = (const matrix &) = default;
			matrix &operator += (const matrix &);
			matrix &operator -= (const matrix &);
			matrix &operator *= (const T &);
			matrix_vector &operator [] (const int &);
			const matrix_vector &operator [] (const int &) const;

			//operators
			friend matrix operator + (const matrix &, const matrix &);
			friend matrix operator + (const matrix &);
			friend matrix operator - (const matrix &, const matrix &);
			friend matrix operator - (const matrix &);
			friend matrix operator * (const matrix &, const matrix &);
			friend matrix operator * (const T &, const matrix &);
			friend matrix operator * (const matrix &, const T &);

			//functions
			matrix transposition() const;
			matrix inversion() const;
			matrix block(int _u, int _d, int _l, int _r) const;	//TEST
			int column() const;
			int row() const;
			int rank() const;	//returning the rank of this matrix
			void print(ostream &_os = cout) const;
			bool identity() const;	//TEST
			bool square() const;	//TODO
			static int set_step_flags(bool _flag);

			//row transformation
			int row_pivot(int _row) const;
			bool row_echelon() const;	//returning if this is a row_echelon matrix
			matrix &row_exchange(const int &_i, const int &_j);
			matrix &row_time(const int &_i, const T &_k);
			matrix &row_time_plus(const int &_i, const T &_k, const int &_j);
			matrix row_echelon_form() const;	//returning a copy of row_echelon matrix

			//column transformation
			int column_pivot(int _column) const;
			bool column_echelon() const;
			matrix &column_exchange(const int &_i, const int &j);
			matrix &column_time(const int &_i, const T &_k);
			matrix &column_time_plus(const int &_i, const T &_k, const int &_j);
			matrix column_echelon_form() const;
	};
	bool matrix::step_by_step;
	bool matrix::vertical = true;
	bool matrix::horizontal = false;

	la_vector dimention_error;
	matrix matrix_dimention_error;

	//declaration of functions related
	matrix &transpose(matrix &_mat);
	matrix transpose(const matrix &_mat);
	matrix &row_eliminate(matrix &_mat, int _row_start = 1, int _column_start = 1);
	matrix row_eliminate(const matrix &_mat);
	matrix &column_eliminate(matrix &_mat, int _row_start = 1, int _column_start = 1);
	matrix column_eliminate(const matrix &_mat);
	matrix &row_reduce(matrix &_mat);
	matrix row_reduce(const matrix &_mat);
	matrix &column_reduce(matrix &_mat);
	matrix volumn_reduce(const matrix &_mat);
	matrix &inverse(matrix &_mat);
	matrix inverse(const matrix &_mat);
	matrix concat(const matrix &_matl, const matrix &_matr);	//TODO
	matrix vconcat(const matrix &_matu, const matrix &_matd);	//TODO
	la_vector solve(const matrix &_mat);

	//definition of the menber functions of la_vector
	//constructors
	la_vector::la_vector() {
		_dimention = 0;
		is_row = false;
	}

	la_vector::la_vector(const int &d) {
		_dimention = d;
		is_row = false;
		_vec.resize(d);
	}

	la_vector::la_vector(const int &d, const bool &r) {
		_dimention = d;
		is_row = r;
		_vec.resize(d);
	}

	la_vector::la_vector(const int &d, T *a) {
		_dimention = d;
		is_row = false;
		_vec.resize(d);
		for(int i = 0; i < d; ++i) {
			_vec[i] = a[i];
		}
	}

	la_vector::la_vector(const int &d, const bool &r, T *a) {
		_dimention = d;
		is_row = r;
		_vec.resize(d);
		for(int i = 0; i < d; ++i) {
			_vec[i] = a[i];
		}
	}
	//assignments
	la_vector &la_vector::operator = (const la_vector &other) {
		_dimention = other.dimention();
		is_row = other.row();
		_vec.resize(_dimention);
		for(int i = 0; i < _dimention; ++i) {
			_vec[i] = other[i+1];
		}
		return *this;
	}

	la_vector &la_vector::operator += (const la_vector& b) {
		if(_dimention != b.dimention() || is_row != b.row()) {
			return dimention_error;
		}
		for(int i = 0; i < _dimention; ++i) {
			_vec[i] += b[i+1];
		}
		
		return *this;
	}

	la_vector &la_vector::operator -= (const la_vector& b) {
		if(_dimention != b.dimention() || is_row != b.row()) {
			return dimention_error;
		}
		for(int i = 0; i < _dimention; ++i) {
			_vec[i] -= b[i+1];
		}
		return *this;
	}

	la_vector &la_vector::operator *= (const T& b) {	//times a number
		for(int i = 0; i < _dimention; ++i) {
			_vec[i] *= b;
		}
		return *this;
	}

	T &la_vector::operator [] (const int &i) {
		return _vec[i-1];
	}

	const T &la_vector::operator [] (const int &i) const {
		return _vec[i-1];
	}

	//operators
	la_vector operator + (const la_vector &a, const la_vector &b){
		if(a.dimention() != b.dimention() || a.row() != b.row()) {
			return dimention_error;
		}
		la_vector ans(a);
		ans += b;
		return ans;
	}

	la_vector operator - (const la_vector &a, const la_vector &b) {
		if(a.dimention() != b.dimention() || a.row() != b.row()) {
			return dimention_error;
		}
		la_vector ans(a);
		ans -= b;
		return ans;
	}

	la_vector operator + (const la_vector &a) {
		return a;
	}

	la_vector operator - (const la_vector &a) {
		la_vector ret = a*(-1);
		return ret;
	}

	la_vector operator * (const T &lamb, const la_vector &a) {
		la_vector ans(a);
		ans *= lamb;
		return ans;
	}

	la_vector operator * (const la_vector &a, const T &lamb) {
		la_vector ans(a);
		ans *= lamb;
		return ans;
	}

	T operator * (const la_vector &a, const la_vector &b) {
		T ans = 0;
		for(int i = 1; i <= a.dimention(); ++i) {
			ans += a[i]*b[i];
		}
		return ans;
	}

	//functions
	la_vector &la_vector::transpose() {		//transpose this vector
		is_row = !is_row;
		return *this;
	}

	bool la_vector::row() const {
		return is_row;
	}

	int la_vector::dimention() const {
		return _dimention;
	}

	//definition of the member functions of matrix
	//constructors
	matrix::matrix() {
		_row = _column = 0;
	}

	matrix::matrix(int n) {
		assert(n > 0);
		_row = _column = n;
		_mat.resize(n);
		for(int i = 1; i <= n; ++i) {
			(*this)[i].vec().resize(n);
			for(int j = 1; j <= n; ++j) {
				(*this)[i][j] = (i == j ? 1 : 0);
			}
		}
	}

	matrix::matrix(const matrix &other) {
		_row = other.row();
		_column = other.column();
		_mat.resize(_row);
		for(int i = 1; i <= _row; ++i) {
			(*this)[i].vec().resize(_column);
			for(int j = 1; j <= _column; ++j) {
				(*this)[i][j] = other[i][j];
			}
		}
	}

	matrix::matrix(const matrix &matl, const matrix &matr, bool ver) {
		if(!ver) {
			if(matl.row() != matr.row()) {
				(*this) = matrix();
				return;
			} else {
				int r = _row = matl.row();
				int cl = matl.column(), cr = matr.column();
				_column = cl+cr;
				_mat.resize(_row);
				for(int i = 1; i <= r; ++i) {
					(*this)[i].vec().resize(_column);
					for(int j = 1; j <= cl; ++j) {
						(*this)[i][j] = matl[i][j];
					}
					for(int j = 1; j <= cr; ++j) {
						(*this)[i][j+cl] = matr[i][j];
					}
				}
			}
		} else {
			if(matl.column() != matr.column()) {
				(*this) = matrix();
				return;
			} else {
				int c = _column = matl.column();
				int ru = matl.row(), rd = matr.row();
				_row = ru+rd;
				_mat.resize(_row);
				for(int i = 1; i <= _row; ++i) {
					(*this)[i].vec().resize(_column);
				}
				for(int i = 1; i <= c; ++i) {
					for(int j = 1; j <= ru; ++j) {
						(*this)[j][i] = matl[j][i];
					}
					for(int j = 1; j <= rd; ++j) {
						(*this)[j+ru][i] = matr[j][i];
					}
				}
			}
		}
	}

	matrix::matrix(const la_vector &vec) {		//copy from a la_vector
		if(vec.row()) {
			_row = 1;
			_column = vec.dimention();
			_mat.resize(1);
			(*this)[1].vec().resize(_column);
			for(int i = 1; i <= _column; ++i) {
				(*this)[1][i] = vec[i];
			}
		} else {
			_row = vec.dimention();
			_column = 1;
			_mat.resize(_row);
			for(int i = 1; i <= _row; ++i) {
				(*this)[i].vec().resize(1);
				(*this)[i][1] = vec[i];
			}
		}
	}

	matrix::matrix(const int &row, const int &column) {
		_row = row;
		_column = column;
		_mat.resize(row);
		for(int i = 1; i <= row; ++i) {
			(*this)[i].vec().resize(column);
			for(int j = 1; j <= column; ++j) {
				(*this)[i][j] = 0;
			}
		}
	}

	/*
	matrix::matrix(const int &row, const int &column, T **array) {
		_row = row;
		_column = column;
		_mat.resize(row);
		for(int i = 1; i <= row; ++i) {
			(*this)[i].vec().resize(column);
			for(int j = 1; j <= column; ++j) {
				(*this)[i][j] = array[i][j];
			}
		}
	}
	*/

	//assignments
	matrix &matrix::operator += (const matrix &b) {
		if(_row != b.row() || _column != b.column()) {
			return matrix_dimention_error;	//error
		}
		for(int i = 1; i <= _row; ++i) {
			for(int j = 1; j <= _column; ++j) {
				(*this)[i][j] += b[i][j];
			}
		}
		return *this;
	}

	matrix &matrix::operator -= (const matrix &b) {
		if(_row != b.row() || _column != b.column()) {
			return matrix_dimention_error;	//error
		}
		for(int i = 1; i <= _row; ++i) {
			for(int j = 1; j <= _column; ++j) {
				(*this)[i][j] -= b[i][j];
			}
		}
		return *this;
	}

	matrix &matrix::operator *= (const T &k) {
		for(int i = 1; i <= _row; ++i) {
			for(int j = 1; j <= _column; ++j) {
				(*this)[i][j] *= k;
			}
		}
		return *this;
	}

	matrix_vector &matrix::operator [] (const int &i) {
		return _mat[i-1];
	}

	const matrix_vector &matrix::operator [] (const int &i) const {
		return _mat[i-1];
	}

	//operators
	matrix operator + (const matrix &a, const matrix &b) {
		matrix ret(a);
		ret += b;
		return ret;
	}

	matrix operator + (const matrix &a) {
		matrix ret(a);
		return ret;
	}

	matrix operator - (const matrix &a, const matrix &b) {
		matrix ret(a);
		ret -= b;
		return ret;
	}

	matrix operator - (const matrix &a) {
		matrix ret(a);
		ret *= -1;
		return ret;
	}

	matrix operator * (const matrix &a, const matrix &b) {
		matrix ret(a.row(), b.column());
		if(a.column() == b.row()) {
			for(int i = 1; i <= a.row(); ++i) {
				for(int j = 1; j <= b.column(); ++j) {
					for(int k = 1; k <= a.column(); ++k) {
						ret[i][j] += a[i][k]*b[k][j];
					}
				}
			}
		}
		return ret;
	}
	
	matrix operator * (const T &k, const matrix &a) {
		matrix ret(a);
		ret *= k;
		return ret;
	}

	matrix operator * (const matrix &a, const T &k) {
		matrix ret(a);
		ret *= k;
		return ret;
	}

	//functions
	matrix matrix::transposition() const {
		matrix ret(*this);
		transpose(ret);
		return ret;
	}

	matrix matrix::inversion() const {
		matrix ret(*this);
		inverse(ret);
		return ret;
	}

	matrix matrix::block(int u, int d, int l, int r) const {
		matrix ret(d-u+1, r-l+1);
		for(int i = u; i <= d; ++i) {
			for(int j = l; j <= r; ++j) {
				ret[i-u+1][j-l+1] = (*this)[i][j];
			}
		}
		return ret;
	}

	int matrix::column() const {
		return _column;
	}

	int matrix::row() const {
		return _row;
	}

	int matrix::row_pivot(int r) const {
		if(r <= 0 || r > _row) {
			return -1;	//error
		}
		for(int i = 1; i <= _column; ++i) {
			if((*this)[r][i] != 0) {
				return i;
			}
		}
		return -1;
	}

	int matrix::column_pivot(int c) const {
		if(c <= 0 || c > _column) {
			return -1;	//error
		}
		for(int i = 1; i <= _row; ++i) {
			if((*this)[i][c] != 0) {
				return i;
			}
		}
		return -1;
	}

	bool matrix::row_echelon() const {
		int now = 0, pre = 0;
		for(int i = 1; i <= _row; ++i) {
			pre = now;
			now = row_pivot(i);
			if(now != -1 && now <= pre) {
				return false;
			}
		}
		return true;
	}

	bool matrix::column_echelon() const {
		int now = 0, pre = 0;
		for(int i = 1; i <= _column; ++i) {
			pre = now;
			now = column_pivot(i);
			if(now != -1 && now <= pre) {
				return false;
			}
		}
		return true;
	}

	int matrix::rank() const {
		int ret = 0, row = _row;
		matrix cpy(*this);
		if(!row_echelon()) {
			row_eliminate(cpy);
		}
		for(int i = 1; i <= row; ++i) {
			if(cpy.row_pivot(i) != -1) {
				ret++;
			}
		}
		return ret;
	}

	void matrix::print(ostream &os) const {
		int row = _row, column = _column, width = -1;
		string ans[row][column];
		for(int i = 1; i <= row; ++i) {
			for(int j = 1; j <= column; ++j) {
				ostringstream is;
				is << (*this)[i][j];
				//is >> ans[i][j];
				ans[i-1][j-1] = is.str();
				width = max(width, (int)ans[i-1][j-1].size());
			}
		}
		width++;
		for(int i = 1; i <= row; ++i) {
			for(int j = 1; j <= column; ++j) {
				os << setw(width) << ans[i-1][j-1];
			}
			os << '\n';
		}
		os << flush;
	}

	bool matrix::identity() const {
		if(_row != _column) {
			return false;
		}
		for(int i = 1; i <= _row; ++i) {
			for(int j = 1; j <= _column; ++j) {
				if((i == j && (*this)[i][j] != 1) || (i != j && (*this)[i][j] != 0)) {
					return false;
				}
			}
		}
		return true;
	}

	bool matrix::square() const {
		return _row > 0 && _row == _column;
	}

	int matrix::set_step_flags(bool flag) {
		if(flag == step_by_step) {
			return 0;
		} else {
			step_by_step = flag;
			return 1;
		}
	}

	//row transformation
	matrix &matrix::row_exchange(const int &i, const int &j) {
		if(i < 0 || j < 0 || i > _row || j > _row) {
			return matrix_dimention_error;	//error
		}
		swap((*this)[i], (*this)[j]);
		if(step_by_step) {
			cout << "R:\t(r" << i << ", r" << j << ")" << endl;
			print();
		}
		return *this;
	}

	matrix &matrix::row_time(const int &r, const T &k) {
		if(r < 0 || r > _row) {
			return matrix_dimention_error;	//error
		}
		for(int i = 1; i <= _column; ++i) {
			(*this)[r][i] *= k;
		}
		if(step_by_step) {
			cout << "R:\t" << k << "*r" << r << endl;
			print();
		}
		return *this;
	}

	matrix &matrix::row_time_plus(const int &r, const T &k, const int &s) {
		if(r <= 0 || s <= 0 || r > _row || s > _row) {
			return matrix_dimention_error;	//error
		}
		for(int i = 1; i <= _column; ++i) {
			(*this)[s][i] += (*this)[r][i]*k;
		}
		if(step_by_step) {
			cout << "R:\t" << k << "*r" << r << "+r" << s << endl;
			print();
		}
		return *this;
	}

	matrix &matrix::column_exchange(const int &i, const int &j) {
		if(i < 0 || j < 0 || i > _column || j > _column) {
			return matrix_dimention_error;	//error
		}
		for(int k = 1; k <= _row; ++k) {
			swap((*this)[k][i], (*this)[k][j]);
		}
		if(step_by_step) {
			cout << "C:\t(c" << i << ", c" << j << ")" << endl;
			print();
		}
		return *this;
	}

	matrix &matrix::column_time(const int &c, const T &k) {
		if(c < 0 || c > _column) {
			return matrix_dimention_error; //error
		}
		for(int i = 1; i <= _row; ++i) {
			(*this)[i][c] *= k;
		}
		if(step_by_step) {
			cout << "C:\t" << k << "*c" << c << endl;
			print();
		}
		return *this;
	}

	matrix &matrix::column_time_plus(const int &c, const T &k, const int &s) {
		if(c <= 0 || s <= 0 || c > _column || s > _column) {
			return matrix_dimention_error; //error
		}
		for(int i = 1; i <= _row; ++i) {
			(*this)[i][s] += (*this)[i][c]*k;
		}
		if(step_by_step) {
			cout << "C:\t" << k << "*c" << c << "+c" << s << endl;
			print();
		}
		return *this;
	}

	matrix matrix::row_echelon_form() const {
		matrix ret(*this);
		row_eliminate(ret);
		return ret;
	}

	matrix matrix::column_echelon_form() const {
		matrix ret(*this);
		column_eliminate(ret);
		return ret;
	}

	//definition of the functions related
	matrix &transpose(matrix &mat) {
		int row = mat.row(), column = mat.column();
		matrix ret(column, row);
		for(int i = 1; i <= row; ++i) {
			for(int j = 1; j <= column; ++j) {
				ret[j][i] = mat[i][j];
			}
		}
		mat = ret;
		return mat;
	}

	matrix transpose(const matrix &mat) {
		matrix ret(mat);
		transpose(ret);
		return ret;
	}

	matrix &row_eliminate(matrix &mat, int rs, int cs) {
		if(rs > mat.row() || cs > mat.column()) {
			return mat;
		}
		bool done = true;
		int row = mat.row();
		for(int i = rs; i <= row; ++i) {
			if(mat[i][cs] != 0) {
				done = false;
				//swap(mat[i][cs], mat[rs][cs]);
				mat.row_exchange(i, rs);
				break;
			}
		}
		if(done) {
			return row_eliminate(mat, rs, cs+1);
		}
		mat.row_time(rs, 1/mat[rs][cs]);
		for(int i = rs+1; i <= row; ++i) {
			mat.row_time_plus(rs, -mat[i][cs], i);
		}
		return row_eliminate(mat, rs+1, cs+1);
	}

	matrix row_eliminate(const matrix &_mat) {
		matrix ret(_mat);
		row_eliminate(ret);
		return ret;
	}

	matrix &column_eliminate(matrix &mat, int rs, int cs) {
		if(rs > mat.row() || cs > mat.column()) {
			return mat;
		}
		bool done = true;
		int column = mat.column();
		for(int i = cs; i <= column; ++i) {
			if(mat[rs][i] != 0) {
				done = false;
				mat.column_exchange(i, cs);
				break;
			}
		}
		if(done) {
			return column_eliminate(mat, rs+1, cs);
		}
		mat.column_time(cs, 1/mat[rs][cs]);
		for(int i = cs+1; i <= column; ++i) {
			mat.column_time_plus(cs, -mat[rs][i], i);
		}
		return column_eliminate(mat, rs+1, cs+1);
	}

	matrix &row_reduce(matrix &mat) {
		if(!mat.row_echelon()) {
			row_eliminate(mat);
		}
		int head;
		for(int i = mat.rank(); i > 0; --i) {
			head = mat.row_pivot(i);
			for(int j = i-1; j > 0; --j) {
				mat.row_time_plus(i, mat[j][head]*(-1), j);
			}
		}
		return mat;
	}

	matrix row_reduce(const matrix &mat) {
		matrix ret(mat);
		row_reduce(ret);
		return ret;
	}

	matrix &column_reduce(matrix &mat) {
		if(!mat.column_echelon()) {
			column_eliminate(mat);
		}
		int head;
		for(int i = mat.rank(); i > 0; --i) {
			head = mat.column_pivot(i);
			for(int j = i-1; j > 0; --j) {
				mat.column_time_plus(i, mat[head][j]*(-1), j);
			}
		}
		return mat;
	}

	matrix column_reduce(const matrix &mat) {
		matrix ret(mat);
		column_reduce(ret);
		return ret;
	}

	matrix &inverse(matrix &mat) {
		if(!mat.square()) {
			return mat;	//error
		}
		matrix ret(mat, matrix(mat.row()));
		int n = ret.row();
#ifdef DEBUG
		cout << n << endl;
#endif
		row_reduce(ret);
		mat = ret.block(1, n, n+1, 2*n);
		return mat;
	}

	matrix inverse(const matrix &mat) {
		matrix ret(mat);
		inverse(ret);
		return ret;
	}

	matrix concat(const matrix &matl, const matrix &matr) {
		matrix ret(matl, matr);
		return ret;
	}

	matrix vconcat(const matrix &matu, const matrix &matd) {
		matrix ret(matu, matd, true);
		return ret;
	}

	la_vector solve(const matrix &mat) {
		int rank, column = mat.column();
		matrix str(mat);
		row_eliminate(str);
		rank = str.rank();
		if(str.row_pivot(rank) == str.column()) {
			return la_vector();	//no solution
		} else {
			for(int i = 1; i <= rank; ++i) {
				if(str.row_pivot(i) != i) {
					return la_vector();	//infinite solution
				}
			}
			if(str.row_pivot(rank) < str.column()-1) {
				return la_vector(); //infinite solution
			}
			row_reduce(str);
			la_vector ret(rank);
			for(int i = 1; i <= rank; ++i) {
				ret[i] = str[i][column];
			}
			return ret;
		}
	}
}
