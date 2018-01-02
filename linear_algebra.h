/* Copyright dyj (c) 2018
 * any codes cannot be used for business
 * template for project special
 * updated on 2018/1/2
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
 *			_row,_column	storing the number of rows/columns of the matrix
 *			_mat			storing the matrix itself
 *		public members:
 *			transpose()	to transpose the matrix and return its reference
 *			column()
 *			row()		to return the number of columns/rows of the matrix
 *			row_exchange(_i, _j)		to exchange the row of _i and _j
 *			row_exchange_step(_i, _j)	to print the step
 *			row_time(_i, _k)			to make each element on row _i time _k
 *			row_time_step(_i, _k)		...
 *			row_time_plus(_i, _k, _j)	to add _k times of row _i to row _j
 *			row_time_plus_step(_i, _k, _j) ...
 *			operator []	to return the certian row in the matrix, range: [1, _row]
 *						WARNING never access the _mat directly, that will make
 *						the program insanely chaotic, for the index of this []
 *						is actually different from that of _mat
 *		functions related:
 *			stairlize(_mat, _row_start, _column_start)
 *						to make _mat a echelon form row equivalent to the original
 *						one. Directly modify on _mat
 *			diagonalize(_mat)
 *						to make _mat to a reduced echelon form
 * */
#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
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
			int _row, _column;
			vector<matrix_vector> _mat;

		public:
			//constructors
			matrix() = default;
			matrix(const matrix &);
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
			matrix &transpose();
			int column() const;
			int row() const;
			int first_nzero(int _row) const;
			bool stairlized() const;	//returning if this is a stairlized matrix
			int rank() const;	//returning the rank of this matrix
			void print() const;

			//row transformation
			matrix &row_exchange(const int &_i, const int &_j);
			matrix &row_exchange_step(const int &_i, const int &_j);
			matrix &row_time(const int &_i, const T &_k);
			matrix &row_time_step(const &_i, const T &_k);
			matrix &row_time_plus(const int &_i, const T &_k, const int &_j);
			matrix &row_time_plus_step(const int &_i, const T &_k, const int &_j);
			matrix stair() const;	//returning a copy of stairlized matrix
	};

	la_vector dimention_error;
	matrix matrix_dimention_error;

	//declaration of functions related
	void stairlize(matrix &_mat, int _row_start = 1, int _column_start = 1);
	void stairlize_step(matrix &_mat, int _row_start = 1, int _column_start = 1);
	void diagonalize(matrix &_mat);
	void diagonalize_step(matrix &_mat);
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
	matrix &matrix::transpose() {
		int bigger = max(_row, _column);
		swap(_row, _column);
		_mat.resize(bigger);
		for(int i = 1; i <= bigger; ++i) {
			(*this)[i].vec().resize(bigger);
		}
		for(int i = 1; i <= bigger; ++i) {
			for(int j = i+1; j <= bigger; ++i) {
				swap((*this)[i][j], (*this)[j][i]);
			}
			(*this)[i].vec().resize(_column);
		}
		_mat.resize(_row);
		return *this;
	}

	int matrix::column() const {
		return _column;
	}

	int matrix::row() const {
		return _row;
	}

	int matrix::first_nzero(int r) const {
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

	bool matrix::stairlized() const {
		int now = 0, pre = 0;
		for(int i = 1; i <= _row; ++i) {
			pre = now;
			now = first_nzero(i);
			if(now != -1 && now <= pre) {
				return false;
			}
		}
		return true;
	}

	int matrix::rank() const {
		int ret = 0, row = _row;
		matrix cpy(*this);
		if(!stairlized()) {
			stairlize(cpy);
		}
		for(int i = 1; i <= row; ++i) {
			if(cpy.first_nzero(i) != -1) {
				ret++;
			}
		}
		return ret;
	}

	void matrix::print() const {
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
				cout << setw(width) << ans[i-1][j-1];
			}
			cout << endl;
		}
	}

	//row transformation
	matrix &matrix::row_exchange(const int &i, const int &j) {
		if(i < 0 || j < 0 || i > _row || j > _row) {
			return matrix_dimention_error;	//error
		}
		swap((*this)[i], (*this)[j]);
		return *this;
	}

	matrix &matrix::row_exchange_step(const int &i, const &j) {
		if(i < 0 || j < 0 || i > _row || j > _row) {
			return matrix_dimention_error;	//error
		}
		swap((*this)[i], (*this)[j]);
		cout << "R:\t(r" << i << ", r" << j << ")" << endl;
		print();
		return *this;
	}

	matrix &matrix::row_time(const int &r, const T &k) {
		if(r < 0 || r > _row) {
			return matrix_dimention_error;	//error
		}
		for(int i = 1; i <= _column; ++i) {
			(*this)[r][i] *= k;
		}
		return *this;
	}

	matrix &matrix::row_time_step(const int &r, const T &k) {
		if(r < 0 || r > _row) {
			return matrix_dimention_error;	//error
		}
		for(int i = 1; i <= _column; ++i) {
			(*this)[r][i] *= k;
		}
		cout << "R:\t" << k << "*r" << r << endl;
		print();
		return *this;
	}

	matrix &matrix::row_time_plus(const int &r, const T &k, const int &s) {
		if(r <= 0 || s <= 0 || r > _row || s > _row) {
			return matrix_dimention_error;	//error
		}
		for(int i = 1; i <= _column; ++i) {
			(*this)[s][i] += (*this)[r][i]*k;
		}
		return *this;
	}

	matrix &matrix::row_time_plus_step(const int &r, const T &k, const int &s) {
		if(r <= 0 || s <= 0 || r > _row || s > _row) {
			return matrix_dimention_error;	//error
		}
		for(int i = 1; i <= _column; ++i) {
			(*this)[s][i] += (*this)[r][i]*k;
		}
		cout << "R:\t" << k << "*r" << r << "+r" << s << endl;
		print();
		return *this;
	}

	matrix matrix::stair() const {
		matrix ret(*this);
		stairlize(ret);
		return ret;
	}


	//definition of the functions related
	void stairlize(matrix &mat, int rs, int cs) {
		if(rs > mat.row() || cs > mat.column()) {
			return;
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
			stairlize(mat, rs, cs+1);
			return;
		}
		mat.row_time(rs, 1/mat[rs][cs]);
		for(int i = rs+1; i <= row; ++i) {
			mat.row_time_plus(rs, -mat[i][cs], i);
		}
		stairlize(mat, rs+1, cs+1);
	}

	void stairlize_step(matrix &mat, int rs, int cs) {
		if(rs > mat.row() || cs > mat.column()) {
			return;
		}
		bool done = true;
		int row = mat.row();
		for(int i = rs; i <= row; ++i) {
			if(mat[i][cs] != 0) {
				done = false;
				//swap(mat[i][cs], mat[rs][cs]);
				mat.row_exchange_step(i, rs);
				break;
			}
		}
		if(done) {
			stairlize_step(mat, rs, cs+1);
			return;
		}
		mat.row_time_step(rs, 1/mat[rs][cs]);
		for(int i = rs+1; i <= row; ++i) {
			mat.row_time_plus_step(rs, -mat[i][cs], i);
		}
		stairlize_step(mat, rs+1, cs+1);
	}

	void diagonalize(matrix &mat) {
		if(!mat.stairlized()) {
			stairlize(mat);
		}
		int head;
		for(int i = mat.rank(); i > 0; --i) {
			head = mat.first_nzero(i);
			for(int j = i-1; j > 0; --j) {
				mat.row_time_plus(i, mat[j][head]*(-1), j);
			}
		}
	}

	void diagonalize_step(matrix &mat) {
		if(!mat.stairlized()) {
			stairlize_step(mat);
		}
		int head;
		for(int i = mat.rank(); i > 0; --i) {
			head = mat.first_nzero(i);
			for(int j = i-1; j > 0; --j) {
				mat.row_time_plus_step(i, mat[j][head]*(-1), j);
			}
		}
	}

	la_vector solve(const matrix &mat) {
		int rank, column = mat.column();
		matrix str(mat);
		stairlize(str);
		rank = str.rank();
		if(str.first_nzero(rank) == str.column()) {
			return la_vector();	//no solution
		} else {
			for(int i = 1; i <= rank; ++i) {
				if(str.first_nzero(i) != i) {
					return la_vector();	//infinite solution
				}
			}
			if(str.first_nzero(rank) < str.column()-1) {
				return la_vector(); //infinite solution
			}
			diagonalize(str);
			la_vector ret(rank);
			for(int i = 1; i <= rank; ++i) {
				ret[i] = str[i][column];
			}
			return ret;
		}
	}
}
