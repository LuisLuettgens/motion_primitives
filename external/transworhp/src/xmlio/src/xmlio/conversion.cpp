#include "conversion.h"

#include <sstream>
#include <cstdlib>
#include <iostream>

std::string Replace(std::string &s, const std::string &r1, const std::string &r2) {

        if (r1 == r2) return s;

        std::string::size_type st;
        while ((st=s.find(r1))!=std::string::npos) {
                s.replace(st,r1.length(),r2);
        }

        return s;
}

std::string Join(const std::vector<std::string> &s, const std::string &j) {

	std::string ret = std::string("");
	std::vector<std::string>::const_iterator it = s.begin();

	for (;it!=s.end();it++) {
		if (it!=s.begin())
			ret = ret + j;
		ret = ret + *it;

	}

	return ret;
}

std::string Join(const std::vector<double> &s, const std::string &j) {

	std::string ret = std::string("");
	std::vector<double>::const_iterator it = s.begin();

	for (;it!=s.end();it++) {
		if (it!=s.begin())
			ret = ret + j;
		ret = ret + ToString(*it);

	}

	return ret;
}

std::string Join(const std::vector<int> &s, const std::string &j) {

	std::string ret;

	int ct = 0;
	for (int value : s) {
		if (ct)	ret = ret + j;
		ret = ret + ToString(value);
		ct++;
	}

	return ret;
}

std::vector<std::string> ToStringArray2(const std::string &s, const std::string &split) {

	int lasti=0;
	std::vector<std::string> ret;

	unsigned int i,j;

	for (i=0;i<s.length();i++) {
		for (j=0;j<split.length();j++) {
			if (s.at(i) == split.at(j)) {
				std::string s1 = s.substr(lasti,i-lasti);
				if (s1.length()) {
					eatwhitespace(s1); // Matthias: neu 9.5.14
					ret.push_back(s1);
				}
				lasti=i+1;
				break;
			}
		}
		if (ret.size()>0) break;
	}


	std::string s1 = s.substr(lasti,s.length());
	if (s1.length()) {
		eatwhitespace(s1); // Matthias: neu 9.5.14
		ret.push_back(s1);
	}
	return ret;
}


std::vector<std::string> ToStringArray(const std::string &s, const std::string &split) {

	int lasti=0;
	std::vector<std::string> ret;

	unsigned int i,j;

	for (i=0;i<s.length();i++) {
		for (j=0;j<split.length();j++) {
			if (s.at(i) == split.at(j)) {
				std::string s1 = s.substr(lasti,i-lasti);
				if (s1.length()) {
					ret.push_back(s1);
				}
				lasti=i+1;
				break;
			}
		}

	}

	std::string s1 = s.substr(lasti,i);
	if (s1.length()) {
		ret.push_back(s1);
	}
	return ret;
}

std::vector<std::string> ToStringArrayEmpty(const std::string &s, const std::string &split) {

	int lasti=0;
	std::vector<std::string> ret;

	unsigned int i,j;

	for (i=0;i<s.length();i++) {
		for (j=0;j<split.length();j++) {
			if (s.at(i) == split.at(j)) {
				std::string s1 = s.substr(lasti,i-lasti);
				/*if (s1.length())*/
				{
					ret.push_back(s1);
				}
				lasti=i+1;
				break;
			}
		}

	}

	std::string s1 = s.substr(lasti,i);
	ret.push_back(s1);

	return ret;
}
std::vector<std::string> ToStringArray(const std::string &s) {

	return ToStringArray(s,std::string(" \t"));
	/*
	int lasti=0;
	std::vector<std::string> ret;

	unsigned int i;

	for (i=0;i<s.length();i++) {
		if (s.at(i)==' ' || s.at(i)=='\t') {
			std::string s1 = s.substr(lasti,i-lasti);
			if (s1.length()) {
				ret.push_back(s1);
			}
			lasti=i+1;
		}

	}

	std::string s1 = s.substr(lasti,i);
	ret.push_back(s1);

	return ret;*/
}


std::vector<float> ToFloatArray(const std::string &s, const std::string &split) {

	int lasti=0;
	std::vector<float> ret;

	unsigned int i,j;

	for (i=0;i<s.length();i++) {
		for (j=0;j<split.length();j++) {
			if (s.at(i)==split.at(j)) {
				std::string s1 = s.substr(lasti,i-lasti);
				if (s1.length()) {
					float a = ToFloat(s1);
					ret.push_back(a);
				}
				lasti=i+1;
				break;
			}
		}

	}

	std::string s1 = s.substr(lasti,i);
	float a = (float) ToDouble(s1);
	ret.push_back(a);

	return ret;
}




std::vector<double> ToDoubleArray(const std::string &s, const std::string &split) {

	int lasti=0;
	std::vector<double> ret;

	unsigned int i,j;

	for (i=0;i<s.length();i++) {
		for (j=0;j<split.length();j++) {
			if (s.at(i)==split.at(j)) {
				std::string s1 = s.substr(lasti,i-lasti);
				if (s1.length()) {
					double a = ToDouble(s1);
					ret.push_back(a);
				}
				lasti=i+1;
				break;
			}
		}

	}

	std::string s1 = s.substr(lasti,i);
	if (s1.length()) {
	double a = ToDouble(s1);
	ret.push_back(a);
	}
	return ret;
}

std::vector<double> ToDoubleArrayEmpty(const std::string &s, const std::string &split) {

	int lasti=0;
	std::vector<double> ret;

	unsigned int i,j;

	for (i=0;i<s.length();i++) {
		for (j=0;j<split.length();j++) {
			if (s.at(i)==split.at(j)) {
				std::string s1 = s.substr(lasti,i-lasti);
				if (s1.length()) {
					double a = ToDouble(s1);
					ret.push_back(a);
				}
else
ret.push_back(0);
				lasti=i+1;
				break;
			}
		}

	}

	std::string s1 = s.substr(lasti,i);
	if (s1.length()) {
	double a = ToDouble(s1);
	ret.push_back(a);
	}
else
ret.push_back(0);

	return ret;
}

std::vector<float> ToFloatArray(const std::string &s) {

	return ToFloatArray(s,std::string(" \t,"));

}

std::vector<double> ToDoubleArray(const std::string &s) {

	return ToDoubleArray(s,std::string(" \t,"));
	/*
	int lasti=0;
	std::vector<double> ret;

	unsigned int i;

	for (i=0;i<s.length();i++) {
		if (s.at(i)==' ' || s.at(i)=='\t' || s.at(i)==',') {
			std::string s1 = s.substr(lasti,i-lasti);
			if (s1.length()) {
				double a = ToDouble(s1);
				ret.push_back(a);
			}
			lasti=i+1;
		}

	}

	std::string s1 = s.substr(lasti,i);
	double a = ToDouble(s1);
	ret.push_back(a);

	return ret;*/
}

std::vector<int> ToIntArray(const std::string &s, const std::string &split) {

	int lasti=0;
	std::vector<int> ret;

	unsigned int i,j;

	for (i=0;i<s.length();i++) {
		for (j=0;j<split.length();j++) {
			if (s.at(i)==split.at(j)) {
				std::string s1 = s.substr(lasti,i-lasti);
				if (s1.length()) {
					int a = ToInt(s1);
					ret.push_back(a);
				}
				lasti=i+1;
				break;
			}
		}

	}

	std::string s1 = s.substr(lasti,i);
	int a = ToInt(s1);
	ret.push_back(a);

	return ret;
}

std::vector<int> ToIntArray(const std::string &s) {

	return ToIntArray(s,std::string(" \t,"));
}


int ToInt(const std::string &s) {
	int a =atoi(s.c_str());
	return a;
}

bool ToBool(const std::string &s) {
	int a =atoi(s.c_str());
	return (a==0)?false:true;
}

float ToFloat(const std::string &s) {
	return (float) ToDouble(s);
}

int floatpoint=0 ; // 1=.   2=,


void InitLocale() {

	std::string spoint("0.2");
	std::string scomma("0,2");

	std::stringstream sspoint(spoint);
	std::stringstream sscomma(scomma);

	double a,b;
	sspoint >> a;
	sscomma >> b;

	//std::cout << a << " " << b << std::endl;

	if (a>=.1) {
	//	std::cout << "atof uses decimal point" << std::endl;
		floatpoint = 1;
	}
	else if (b>=.1) {
		std::cout << "@Jan: atof uses decimal comma" << std::endl;
		floatpoint = 2;
	}
	else {
		std::cout << "@Jan: Can't transform numbers." << std::endl;
	}
}


double ToDouble(const std::string &s) {

	if (floatpoint==0) InitLocale();

	std::string s2(s);

	if (floatpoint==2) {
		std::string::size_type st = s2.find('.');
		if (st!=std::string::npos) {
			s2.replace(st,1,std::string(","));
		}
	}

	if (floatpoint==1) {
		std::string::size_type st = s2.find(',');
		if (st!=std::string::npos) {
			s2.replace(st,1,std::string("."));
		}
	}

	std::stringstream ss(s2);
	double ret=0;
	ss>>ret;

	return ret;
}

#ifdef WIN32
std::string ToString(unsigned long a) {
	std::stringstream s;
	s << a;
	return s.str();
}
#endif

std::string ToString(int a) {
	std::stringstream s;
	s << a;
	return s.str();
}
std::string ToString(size_t a) {
	std::stringstream s;
	s << a;
	return s.str();
}
std::string ToString(double a) {
	std::stringstream s;
	s << a;
	return s.str();
}

std::string ToStringFakePoint(double a) {
	std::string s = ToString(a);
	size_t e = s.find('e');
	if (e == std::string::npos) {
		size_t i = s.find('.');
		if (i == std::string::npos) s += ".0";
	}
	return s;
}

std::string ToString(double a, int width, int prec) {
	std::stringstream s;
	s.setf(std::ios::fixed, std::ios::floatfield);
	s.precision(prec);
	s.width(width);
	s << a;
	return s.str();
}

std::string ToHexString(unsigned long a) {
	std::stringstream s;
	s << std::hex << a;
	return s.str();
}

std::string ToScientificString(double a) {
	std::stringstream s;
	s.setf(std::ios::scientific);
//	s.precision(prec);
//	s.width(width);
	s << a;
	return s.str();
}
std::string ToString(char a) {
	std::stringstream s;
	s << "'" << a << "'";
	return s.str();
}


std::vector<double> ToMatlabArray(const std::string &achse) {

	std::vector<double> a;

	std::vector<double> matlab = ToDoubleArray(achse,":");
	if (matlab.size()==3) {

		for (double i=matlab[0];i<=matlab[2];i+=matlab[1]) {
			a.push_back(i);
		}
	}
	else if (matlab.size() == 2) {

		for (double i = matlab[0]; i <= matlab[1]; i++) {
			a.push_back(i);
		}
	}
	else {
		a = ToDoubleArray(achse,",");
	}
	return a;
}


std::string eatwhitespace(std::string &c) {

	std::string::iterator it = c.begin();
	for (;it!=c.end();) {

		if (*it==' ' || *it==10 || *it==9 || *it=='\t') {
			c.erase(0,1);
			it=c.begin();
		} else
			break;
	}

	std::string::reverse_iterator rit = c.rbegin();
	for (;rit!=c.rend();) {
//		std::cout << c << ">>"<< (char) *rit << "<<" << (int)(*rit) << std::endl;

		if (*rit==' ' || *rit==10 || *rit==9 || *rit=='\t' || *rit==13) {
			c.erase(c.size()-1,1);
			rit=c.rbegin();
		} else
			break;
	}

	return c;
}
