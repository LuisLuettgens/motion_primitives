#include "exception.h"

#include <iostream>
#ifdef WIN32
#include "windows.h"
#endif

using namespace std;


/** Textcolor Implementation (1 Argument manipulator) */
class textc_int {
public:
	textc_int(std::ostream& (*f)(std::ostream&,int), int t) : _fp(f), _tp(t) {}

	std::ostream& (* _fp)(std::ostream&,int);
	int _tp;
	friend std::ostream& operator<<(std::ostream& s, const textc_int & sm);

};

textc_int textc( int l );
std::ostream& textc ( std::ostream& os, int i );



std::ostream& textc ( std::ostream& os, int i ) {

#ifdef WIN32
	switch (i) {
	case 0:    // White on Black
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
		                        FOREGROUND_INTENSITY | FOREGROUND_RED |
		                        FOREGROUND_GREEN | FOREGROUND_BLUE);
		break;
	case 1:    // Red on Black
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
		                        FOREGROUND_INTENSITY | FOREGROUND_RED);
		break;
	case 2:    // Blue on Black
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
		                        FOREGROUND_INTENSITY | FOREGROUND_BLUE);
		break;
	case 3:    // Green on Black
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
		                        FOREGROUND_INTENSITY | FOREGROUND_GREEN);
		break;
	case 4:    // Yellow on Black
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
		                        FOREGROUND_INTENSITY | FOREGROUND_RED |
		                        FOREGROUND_GREEN);
		break;
	case 5:    // Magenta on Black
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
		                        FOREGROUND_INTENSITY | FOREGROUND_RED |
		                        FOREGROUND_BLUE);
		break;

	}
#else
	switch(i) {
	case 0:
		os << "\033[30m";
		break;
	case 1: //rot
		os << "\033[31m";
		break;
	case 2: //blau
		os << "\033[32m";
		break;
	case 3: //grün
		os << "\033[34m";
		break;
	case 4: //grau
		os << "\033[37m";
		break;
	case 5: //dk red
		os << "\033[33m";
		break;
	}
#endif
	return os;
}

textc_int textc( int l ) {
	return textc_int ( textc, l );
}

std::ostream& operator<<(std::ostream& s, const textc_int & sm) {
	(*sm._fp)(s,sm._tp);
	return s;
}


std::ostream &operator <<(std::ostream &os, const Exception &e) {
	os << endl;
	os << textc(1);
	os << "--------------------< EXCEPTION >--------------------" << endl;
	os << e.what.str() << endl;
	os << "-----------------------------------------------------" << textc(0) << endl;
	return os;
}


void Exception::warn(int col) {
	cout << endl;
	cout << textc(col);
	cout << "---------------------< WARNING >---------------------" << endl;
	cout << what.str() << endl;
	cout << "-----------------------------------------------------" << textc(0) << endl;
	
}
