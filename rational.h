/* Code of Project Calculatica (R) by DYJ
 * Headfile : rational
 * usage : Basis of type value
 * namespace : rat
 * In year 2015
 * updated on 2015/8/14
 * updated on 2017/10/11
 */
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <cassert>

namespace rat{
	typedef int hint;
	typedef double hflo;

	class rational{
		private:
			hint mother = 0, son = 0;
			char sign = 1;
			void maintain();
		public:
			//All constructors
			rational() = default;
			rational(const rational&) = default;
			rational(const hint &, const hint & = 1);
			rational(const hflo &);
			//TODO some more constructors must be created.

			//operator for assignment
			rational &operator = (const rational&) = default;
			rational &operator += (const rational&);
			rational &operator -= (const rational&);
			rational &operator *= (const rational&);
			rational &operator /= (const rational&);
			rational &operator %= (const rational&);

			//operator for IO
			//TODO input functions must be made;
			friend std::ostream &operator << (std::ostream &, const rational&);
			friend std::istream &operator >> (std::istream &, const rational&);

			//other operator
			friend bool operator < (const rational&, const rational&);
			friend rational operator - (const rational&);
			friend rational operator + (const rational&, const rational&);
			friend rational operator * (const rational&, const rational&);
			friend rational operator % (const rational&, const rational&);

			//functions
			rational reciprocal() const;
			rational abs() const;
			hint floor() const;
			hint ceil() const;
			hflo float_t() const;
	};

	//declearation of all functions
	hint gcd(const hint &a, const hint &b);
	hint lcm(const hint &a, const hint &b);
	bool operator > (const rational&, const rational&);
	bool operator <= (const rational&, const rational&);
	bool operator >= (const rational&, const rational&);
	bool operator == (const rational&, const rational&);
	bool operator != (const rational&, const rational&);
	rational operator + (const rational&);
	rational operator - (const rational&, const rational&);
	rational operator / (const rational&, const rational&);
	rational reciprocal(const rational&);
	rational abs(const rational&);
	rational istr_rat(const std::string &);	//int_string to rational
	rational fstr_rat(const std::string &);	//float_string to rational
	rational fracstr_rat(const std::string &);	//fraction_string to rational
	rational str_rat(const std::string &);	//all types of strings above to rational
	hint floor(const rational&);
	hint ceil(const rational&);

	//defination of all functions
	hint gcd(const hint &a, const hint &b){
		return b == 0 ? a : gcd(b, a%b);
	}

	hint lcm(const hint &a, const hint &b) {
		return a*b/gcd(a, b);
	}

	bool operator < (const rational &a, const rational &b) {
		if(a.sign != b.sign) {
			return a.sign < b.sign;
		} else if(a.sign == 0) {
			return false;
		} else {
			hint l = lcm(a.mother, b.mother);
			return (a.son*l/a.mother < b.son*l/b.mother)^(a.sign == -1);
		}
	}

	bool operator > (const rational &a, const rational &b) {
		return b < a;
	}

	bool operator <= (const rational &a, const rational &b) {
		return !(b < a);
	}

	bool operator >= (const rational &a, const rational &b) {
		return !(a < b);
	}

	bool operator == (const rational &a, const rational &b) {
		return !((a < b) || (b < a));
	}

	bool operator != (const rational &a, const rational &b) {
		return (a < b) || (b < a);
	}

	rational operator + (const rational &a) {
		return a;
	}

	rational operator - (const rational &a) {
		rational ret = a;
		ret.sign *= -1;
		return ret;
	}

	rational operator + (const rational &a, const rational &b) {
		hint pro = a.sign*b.sign;
		if(pro == 0) {
			if(a == 0) {
				return b;
			} else {
				return a;
			}
		}
		rational ret;
		ret.mother = lcm(a.mother, b.mother);
		ret.son = a.sign*a.son*ret.mother/a.mother+b.sign*b.son*ret.mother/b.mother;
		ret.maintain();
		return ret;
	}

	rational operator - (const rational &a, const rational &b) {
		return a+(-b);
	}

	rational operator * (const rational &a, const rational &b) {
		rational ret;
		ret.son = a.son*b.son*a.sign*b.sign;
		ret.mother = a.mother*b.mother;
		ret.maintain();
		return ret;
	}

	rational operator / (const rational &a, const rational &b) {
		return a*b.reciprocal();
	}

	rational operator % (const rational &a, const rational &b) {
		assert(b != 0);
		hint ti = (a/b).floor();
		return a-b*ti;
	}

	rational &rational::operator += (const rational &a) {
		*this = *this+a;
		return *this;
	}

	rational &rational::operator -= (const rational &a) {
		*this = *this-a;
		return *this;
	}

	rational &rational::operator *= (const rational &a) {
		*this = *this*a;
		return *this;
	}
	
	rational &rational::operator /= (const rational &a) {
		*this = *this/a;
		return *this;
	}
	
	rational &rational::operator %= (const rational &a) {
		*this = *this%a;
		return *this;
	}

	rational reciprocal(const rational &a) {
		return a.reciprocal();
	}

	rational abs(const rational &a) {
		return a.abs();
	}

	rational istr_rat(const std::string &str) {
		rational ret = 0, bse = 1;
		for(auto i = str.rbegin(); i != str.rend() && isdigit(*i); ++i) {
			ret += bse*(*i-'0');
			bse *= 10;
		}
		if(str[0] == '-') {
			ret *= -1;
		}
		return ret;
	}

	rational fstr_rat(const std::string &str) {
		rational ret = 0, b1 = 1, b2(1, 10);
		size_t dot = str.find('.');
		for(int i = dot-1; i >= 0 && isdigit(str[i]); --i) {
			ret += b1*(str[i]-'0');
			b1 *= 10;
		}
		for(size_t i = dot+1; i < str.size(); ++i) {
			ret += b2*(str[i]-'0');
			b2 /= 10;
		}
		if(str[0] == '-') {
			ret *= -1;
		}
		return ret;
	}

	rational fracstr_rat(const std::string &str) {
		hint son = 0, mot = 0, bse = 1;
		size_t sls = str.find('/');
		for(int i = sls-1; i >= 0 && isdigit(str[i]); --i) {
			son += bse*(str[i]-'0');
			bse *= 10;
		}
		bse = 1;
		for(size_t i = str.size()-1; i > sls; --i) {
			mot += bse*(str[i]-'0');
			bse *= 10;
		}
		if(str[0] == '-') {
			son *= -1;
		}
		return rational(son, mot);
	}

	rational str_rat(const std::string &str) {
		rational ret = 0;
		if(str.find('.') != std::string::npos) {
			ret = fstr_rat(str);
		} else if(str.find('/') != std::string::npos) {
			ret = fracstr_rat(str);
		} else {
			ret = istr_rat(str);
		}
		return ret;
	}

	hint floor(const rational &a) {
		return a.floor();
	}

	hint ceil(const rational &a) {
		return a.ceil();
	}

	void rational::maintain(){
		if(son == 0){
			sign = mother = 0;
		} else {
			sign = mother*son > 0 ? 1 : -1;
			hint t = gcd(mother, son);
			mother = std::abs(mother/t);
			son = std::abs(son/t);
		}
	}

	rational::rational(const hint &devidend, const hint &divisor){
		assert(divisor != 0);
		son = devidend;
		mother = divisor;
		maintain();
	}

	std::ostream &operator << (std::ostream &os, const rational& r) {
		if(r.sign < 0){
			os << '-' << std::flush;
		} else if(r.sign == 0) {
			os << '0' << std::flush;
			return os;
		}
		os << r.son << std::flush;
		if(r.mother != 1) {
			os << '/' << r.mother << std::flush;
		}
		return os;
	}

	std::istream &operator >> (std::istream &is, rational &r) {
		std::string str;
		is >> str;
		r = str_rat(str);
		return is;
	}

	rational rational::reciprocal() const {
		rational ret = *this;
		std::swap(ret.son, ret.mother);
		return ret;
	}

	rational rational::abs() const {
		rational ret = *this;
		ret.sign = std::abs(ret.sign);
		return ret;
	}

	hint rational::floor() const {
		if(sign == 0) {
			return 0;
		} else if(sign == 1) {
			return son/mother;
		} else {
			if(son%mother == 0) {
				return sign*son/mother;
			} else {
				return sign*(son/mother+1);
			}
		}
	}

	hint rational::ceil() const {
		if(sign == 0) {
			return 0;
		} else if(sign == -1) {
			return sign*son/mother;
		} else {
			if(son%mother == 0) {
				return son/mother;
			} else {
				return son/mother+1;
			}
		}
	}

	hflo rational::float_t() const {
		return (hflo)son/mother;
	}
}
