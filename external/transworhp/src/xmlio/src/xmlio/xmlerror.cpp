#include "xmlerror.h"
#include <iostream>
#include <sstream>
#include "textout.h"

std::string XMLError::GetErrorText() const {

	std::stringstream ret;
	
	switch (code) {
	case XML_ATTRIBUTE_NOT_QUOTED:
		ret << "Attribute '"<< info <<"' not quoted";
		break;
	case  XML_NO_ELEMENT_NAME:
		ret << "Element name missing";
		break;
	case XML_ATTRIBUTE_WITHOUT_ARGUMENTS:
		ret << "Attribute '"<< info <<"' requires an argument";
		break;
	case XML_EQUAL_EXPECTED:
		ret << "Missing '=' after attribute '"<< info <<"'";
		break;
	case XML_NO_ROOT:
		ret << "No root element";
		break;
	case XML_ELEMENT_NOT_CLOSED:
		ret << "Element '"<< info <<"' isn't closed";
		break;
	case XML_WRONG_ROOT_NAME:
		ret << "DOCTYPE says root element should be '"<< info <<"', but it isn't";
		break;
	case XML_UNEXPECTED_CHARACTER:
		ret << "Unexpected character '" << info << "'";
		break;
	case XML_FILE_NOT_FOUND:
		ret << "File '"<< info <<"' wasn't found, sorry.";
		break;


	default:
		ret << "No error code ";



	}

	return ret.str();
	

}

std::ostream &operator << (std::ostream &os, const XMLError &e) {

	os << textcolor(1) << "XMLError: " <<textcolor(0)
	<< e.GetErrorText()  << " (Line " << e.line << ")." << endline;
	return os;

}

