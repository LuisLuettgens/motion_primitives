#include "xmlio.h"
#include "textout.h"
#include "xmlerror.h"
#include <vector>
#include <string>
#include <sstream>
#include "conversion.h"

/*
char abc[][10] = {"false","true","no","yes","off","on"};

class XMLBase {

public:
	XMLBase() {
		XMLNode::BOOLEAN.push_back(abc[0]);
		XMLNode::BOOLEAN.push_back(abc[1]);

		XMLNode::LOGICAL.push_back(abc[2]);
		XMLNode::LOGICAL.push_back(abc[3]);

		XMLNode::SWITCH.push_back(abc[4]);
		XMLNode::SWITCH.push_back(abc[5]);
	}

	~XMLBase() {
		XMLNode::BOOLEAN.clear();
		XMLNode::LOGICAL.clear();
		XMLNode::SWITCH.clear();
	}

};

XMLBase xmlbase;*/




XMLParser::XMLParser() {

	root = nullptr;

}

XMLParser::~XMLParser() {

	std::vector<XMLError*>::iterator it = xmlerrors.begin();
	for (;it!=xmlerrors.end();it++) {
		delete *it;
	}
	xmlerrors.clear();

	delete root;
}

XMLNode *XMLParser::GetXMLDeclaration() const {
	return root ? root->GetTypedChild(1) : nullptr;
}

XMLNode *XMLParser::GetDoctype() const {
	return root ? root->GetTypedChild(3) : nullptr;
}

XMLNode *XMLParser::GetDocument() const {
	return root ? root->GetTypedChild(0) : nullptr;
}

XMLNode *XMLParser::ParseString(const std::string &text) {

	std::vector<XMLError*>::iterator it = xmlerrors.begin();
	for (;it!=xmlerrors.end();it++) {
		delete *it;
	}
	xmlerrors.clear();

	line = 1;
	std::stringstream is(text);

	root = new XMLNode(-1);
	if (!is) {
		AddError(XML_FILE_NOT_FOUND, "STRING");
		return nullptr;
	}

	return ParseStream(is);
}

XMLNode *XMLParser::Parse(const std::string &filename) {

	std::vector<XMLError*>::iterator it = xmlerrors.begin();
	for (;it!=xmlerrors.end();it++) {
		delete *it;
	}
	xmlerrors.clear();

	line = 1;
	std::ifstream is(filename.c_str());

	root = new XMLNode(-1);
	if (!is) {
		AddError(XML_FILE_NOT_FOUND, filename);
		return nullptr;
	}

	return ParseStream(is);
}



XMLNode *XMLParser::ParseStream(std::istream &is) {

	char c;
	while (is) {
		is.get(c);

		if (!is)
			break;

		switch (c) {
		case '<': {
				XMLNode *node = new XMLNode;
				node->Parse(is, this);
				root->AddChild(node);
				break;
			}
		case '\n':
			line++;
			break;
		case 13: // für Windows
			break;
		default:
			//std::cout << c << (int)c;
			std::string a = std::string("") + c;
			if (c<32)
				a = std::string("CHAR(")+ToString((int)c)+")";
			AddError(XML_UNEXPECTED_CHARACTER,a);
			break;
		}
	}

	XMLNode *xml = GetXMLDeclaration();
	if (xml) {
		/*std::cout << "Parser: XML        "
		<< xml->GetName() << std::endl;
		std::cout << "Parser: Version    "
		<< xml->GetAttribute(std::string("version")) << std::endl;
		std::cout << "Parser: encoding   "
		<< xml->GetAttribute(std::string("encoding")) << std::endl;
		std::cout << "Parser: standalone "
		<< xml->GetAttribute(std::string("standalone")) << std::endl;*/
	}
	XMLNode *doctype = GetDoctype();
	if (doctype) {
		/*std::cout << "Parser: Doctype    "
		<< doctype->GetName() << std::endl;
		std::cout << "Parser: id         "
		<< doctype->GetAttribute(std::string("id")) << std::endl;
		std::cout << "Parser: system     "
		<< doctype->GetAttribute(std::string("SYSTEM")) << std::endl;*/

		if (doctype->GetAttribute(std::string("id")) != GetDocument()->GetName()) {
			// Check if doctype->id == root->name; (only if doctype is used)
			AddError(XML_WRONG_ROOT_NAME,
			         doctype->GetAttribute(std::string("id")));
		}
	}

	// Open doctype->system;

	if (!xmlerrors.empty()) {
		return nullptr;
	}

	return root;
}

std::string XMLParser::getline(std::istream &is, char c) {

	std::string ret;

	//getline(is,ret,c);

	char e;

	int kommentar=0;
	int klammer=0;
	bool in_quote=false;
	bool in_comment = false;
	char commentstart[] = "!--";
	char commentend[] = "-->";

	while (is) {

		is.get(e);

		if (e=='\n') {
			line ++;
		}

		if (e=='[' && !in_quote)
			klammer++;

		if (e==']' && !in_quote)
			klammer--;

		if (e=='"')
			in_quote = !in_quote;

		if (!in_comment) {
			if (e==commentstart[kommentar] && kommentar<3) {
				kommentar++;
				if (kommentar==3) {
					kommentar=0;
					in_comment=true;
				}
			} else {
				kommentar=0;
			}
		} else {
			if (e==commentend[kommentar] && kommentar<3) {
				kommentar++;
				if (kommentar==3) {
					kommentar=0;
					in_comment=false;
				}
			} else {
				kommentar=0;
			}
		}


		if (e==c && klammer==0 && !in_comment)
			break;

		ret += e;

	}

	return ret;

}

bool XMLParser::EndOfFile(std::istream &is) const {
	return (!is);
}

void XMLParser::AddError(int c, const std::string &s) {

	XMLError *error = new XMLError(c,line,0,s);
	xmlerrors.push_back(error);

}

void XMLParser::GetError(std::ostream &os) const {

	std::vector<XMLError*>::const_iterator it = xmlerrors.begin();
	for (;it!=xmlerrors.end();it++) {
		os << **it;
	}
}
