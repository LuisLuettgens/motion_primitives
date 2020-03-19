#include "xmlio.h"
#include "textout.h"
#include "xmlerror.h"
#include "conversion.h"
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>

/*std::string *XMLNode::BOOLEAN[][6] = {"false", "true", ""};
std::string* XMLNode::LOGICAL[][6] = {"no", "yes", ""};
std::string *XMLNode::SWITCH[][6] = {"off", "on", ""};
*/
/*
 std::vector<char*>XMLNode::BOOLEAN;
 std::vector<char*>XMLNode::LOGICAL;
 std::vector<char*>XMLNode::SWITCH;*/
char XMLNode::BOOLEAN[3][10] = {"false", "true", ""};
char XMLNode::LOGICAL[3][10] = {"no", "yes", ""};
char XMLNode::SWITCH[3][10] = {"off", "on", ""};


XMLNode::XMLNode(int t) : type(t) {}

XMLNode::XMLNode(const std::string &str)
		: nodename(str), type(0) {}

XMLNode::XMLNode(const std::string &str, const std::string &t)
		: nodename(str), type(0), text(t) {}

XMLNode::~XMLNode() {
	std::vector<XMLNode*>::iterator it = child.begin();
	for (;it!=child.end(); it = child.begin()) {
		delete *it;
		child.erase(it);
	}
}

void XMLNode::RemoveAllChilds() {
	std::vector<XMLNode*>::iterator it = child.begin();
	for (;it!=child.end(); it = child.begin()) {
		delete *it;
		child.erase(it);
	}
}


void XMLNode::RemoveAllComments() {
	std::vector<XMLNode*>::iterator it = child.begin();
	for (; it != child.end();) {
		if ((*it)->IsComment()) {
			delete *it;
			child.erase(it);
			it = child.begin();
		}
		else {
			(*it)->RemoveAllComments();
			it++;
		}
	}
}


bool XMLNode::Parse(std::istream &is, XMLParser *parser) {

	std::string a = parser->getline(is, '>');

	if (a.at(0)=='/') {
		// Prüfen, ob start=endtag?
		a = a.erase(0,1);
		nodename = a;
		return false;
	}

	bool hasChilds=true;

	// Automatisches End-Tag
	if (a.at(a.size()-1) == '/') {
		hasChilds=false;
		a = a.erase(a.size()-1,std::string::npos);
	}

	// XML Deklaration erkennen
	if (a.at(0) == '?' && a.at(a.size()-1) == '?') {
		hasChilds=false;
		a = a.erase(a.size()-1,std::string::npos);
		a = a.erase(0,1);
		type = 1;
	}

	// Kommentar
	if (a.substr(0,3) == "!--") {
		if (a.substr(a.size()-2,2)=="--") {
			text = a;
			type = 2;
			return true;
		}
		return false;
	}

	// Doctype
	if (a.at(0) == '!') {
		type = 3;
		return doctypeparse(parser,a);
	}

	//std::cout << a << std::endl;

	std::vector<std::string> ss = attrsplit(a);
	if (!attrparse(parser, ss)) {
		return false;
	}

	if (!hasChilds)
		return true;

	while (!parser->EndOfFile(is)) {

		std::string c = parser->getline(is, '<');

		eatwhitespace(c);
		text += c;

		XMLNode *n = new XMLNode;
		if (n->Parse(is, parser)) {
			child.push_back(n);
		} else {
			if (nodename!=n->nodename) {
				parser->AddError(XML_ELEMENT_NOT_CLOSED,nodename);
				return false;
			}
			delete n;
			return true;
		}

	}

	return true;
}

void printattributepair(std::ostream &os, const std::string &first, const std::string &second) {
	os << " " << textcolor(3) << bold << first << normal << textcolor(0)
		<< "="
		<< textcolor(2)
		<< quote << second << quote << textcolor(0);
}

std::ostream &operator <<(std::ostream &os, const std::map<std::string, std::string> &m) {

	// version extra.
	// id extra

	std::map<std::string, std::string>::const_iterator it_version = m.find("version");
	if (it_version != m.end()) {
		printattributepair(os, it_version->first, it_version->second);
	}

	std::map<std::string, std::string>::const_iterator it_id = m.find("id");
	if (it_id != m.end()) {
		printattributepair(os, it_id->first, it_id->second);
	}

	std::map<std::string, std::string>::const_iterator it1 = m.begin();
	for (; it1 != m.end(); it1++) {
		if (it_version != it1 && it_id != it1) {
			printattributepair(os, it1->first, it1->second);
		}
	}

	return os;
}


void XMLNode::Debug(std::ostream &os, TextOutputType_e e) {
	DebugStream d(e);
	d << beginfile;
	d << *this << std::endl;
	d << endfile;
	os << d.GetString();
}



std::ostream &operator <<(std::ostream &os, const XMLNode &node) {

	os << tab;

	std::string text = node.text;

	if (node.type==-1) {
		os << removetab;
	}
	if (node.type==0) {

		if (text=="" && node.child.size()==0) {
			os << textcolor(3) << lt << node.nodename << node.attribute;

			DebugStream* deb = dynamic_cast<DebugStream*>(&os);
			if (deb->GetTextOutputType() == RTF) {
				std::string::size_type st = text.find('\\');
				while (st!=std::string::npos) {
					text.replace(st,1,std::string("\\\\"));
					st = text.find('\\',st+2);
				}
			}

			os << textcolor(3) << "/" << gt << textcolor(0)<< endline;


		} else {
			os << textcolor(3) << lt << node.nodename << node.attribute;

			DebugStream* deb = dynamic_cast<DebugStream*>(&os);
			if (deb->GetTextOutputType() == RTF) {
				std::string::size_type st = text.find('\\');
				while (st!=std::string::npos) {
					text.replace(st,1,std::string("\\\\"));
					st = text.find('\\',st+2);
				}
			}
			os << textcolor(3) << gt << textcolor(0) << text;
		}

	}
	if (node.type==1) {
		os << textcolor(4) << lt << "?" << textcolor(3)
		<< bold << node.nodename << normal << node.attribute
		<< textcolor(4) << "?" << gt << textcolor(0) << endline;
		return os;
	}
	if (node.type==2) {
		// node.text evtl. für HTML ohne <>><<><
		std::string a = node.text;
		DebugStream *deb = dynamic_cast<DebugStream*>(&os);
		if (deb->GetTextOutputType()==HTML) {
			replacebrackets(a);
		}

		if (deb->GetTextOutputType()==CPP) {
			replacequote(a);
		}

		if (deb->GetTextOutputType()==CONSOLE) {
			std::vector<std::string> aa = ToStringArray(a,std::string("\n"));
			//os << textcolor(4) << oblique << lt;

			std::vector<std::string>::iterator it = aa.begin();
			for (;it!=aa.end();it++) {
				if (it == aa.begin()) {
					os << textcolor(4) << oblique << lt;
				} else {
					os << textcolor(4) << oblique;
				}

				os << *it;

				if (it+1 == aa.end()) {
					os << gt << normal << textcolor(0) << endline;
				} else {
					os << normal << textcolor(0) << endline;
				}
			}

			//os << gt << normal << textcolor(0) << endline;

		} else {

			os << textcolor(4) << oblique << lt << a << gt << normal <<
			textcolor(0) << endline;

		}


		return os;
	}
	if (node.type==3) {
		// <!
		os << textcolor(3) << lt << "!" << textcolor(5) << node.nodename << textcolor(0);
		if (node.GetAttribute(std::string("id"))!="") {
			os << " " << node.GetAttribute(std::string("id"));
		}
		if (node.GetAttribute(std::string("SYSTEM"))!="") {
			std::string a = node.GetAttribute(std::string("SYSTEM"));

			os << textcolor(5) << " SYSTEM " << textcolor(2)
			<< quote << std::string(a,1,a.length()-2) << quote;
		}
		if (node.GetAttribute(std::string("PUBLIC"))!="") {
			os << textcolor(5) << " PUBLIC " << textcolor(2)
			<< node.GetAttribute(std::string("PUBLIC"));
			os << " " << node.GetAttribute(std::string("PUBLIC_URL"));
		}
		if (node.GetAttribute(std::string("LOCAL"))!="") {

			std::string a = node.GetAttribute(std::string("LOCAL"));
			DebugStream* deb = dynamic_cast<DebugStream*>(&os);
			if (deb->GetTextOutputType()==HTML) {
				replacebrackets(a);
			}
			if (deb->GetTextOutputType()==CPP) {
				replacequote(a);
			}
			os << textcolor(3) << " " << a;

		}

		os << textcolor(3) << gt << textcolor(0) << endline;

		return os;
	}

	if (node.child.size()) {
		if (node.type!=-1)
			os << endline;

		os << addtab;


		std::vector<XMLNode*>::const_iterator it = node.child.begin();
		for (;it!=node.child.end();it++) {
			os << **it;
		}

		os << removetab;

		os << tab;
	}

	if ((node.child.size() || text!="") && node.type==0 ) {
		os << textcolor(3) << lt << "/" << node.nodename
		<< gt << textcolor(0)<< endline;
	}

	return os;
}

XMLNode *XMLNode::GetFirstChild(const std::string &name)  {
	it_tmp = child.begin();
	for (;it_tmp!=child.end();it_tmp++) {
		if ((*it_tmp)->GetName()==name)
			return *it_tmp;
	}

	return nullptr;
}

XMLNode *XMLNode::GetChildWithType(const std::string &name, const std::string &type) const {

	std::vector<XMLNode*>::const_iterator tmp = child.begin();

	for (;tmp!=child.end();tmp++) {
		if ((*tmp)->GetName()==name) {
			if ((*tmp)->GetAttribute(std::string("type"))==type) {
				return *tmp;
			}
		}
	}

	return nullptr;
}

XMLNode *XMLNode::GetChildWithId(const std::string &name, const std::string &id) const {

	std::vector<XMLNode*>::const_iterator tmp = child.begin();

	for (;tmp!=child.end();tmp++) {
		if ((*tmp)->GetName()==name) {
			if ((*tmp)->GetAttribute(std::string("id"))==id) {
				return *tmp;
			}
		}
	}

	return nullptr;
}

XMLNode * XMLNode::AddChildAfter(XMLNode * n, XMLNode * par) {
	for (auto it = child.begin(); it != child.end(); ++it) {
		if ((*it) == par) {
			child.insert(++it, n);
			it_tmp = child.begin();
			return n;
		}
	}
	return nullptr;
}

bool XMLNode::RemoveChild(const std::string &name) {
	std::vector<XMLNode*>::iterator it = child.begin();
	for (;it!=child.end();it++) {
		if ((*it)->GetName()==name) {
			delete *it;
			child.erase(it);
			return true;
		}
	}

	return false;
}


XMLNode *XMLNode::GetNextChild(const std::string &name)  {

	++it_tmp;

	for (;it_tmp!=child.end();it_tmp++) {
		if ((*it_tmp)->GetName()==name)
			return *it_tmp;
	}

	return nullptr;
}

XMLNode *XMLNode::GetFirstChild()  {

	it_tmp = child.begin();
	if (it_tmp!=child.end())
		return *it_tmp;

	return nullptr;
}

XMLNode *XMLNode::GetNextChild()  {

	++it_tmp;
	if (it_tmp!=child.end())
		return *it_tmp;

	return nullptr;
}

XMLNode *XMLNode::GetTypedChild(int type) const {

	std::vector<XMLNode*>::const_iterator it = child.begin();
	while (it!=child.end()) {
		if ((*it)->type==type)
			return *it;
		it++;
	}
	return nullptr;
}

std::string XMLNode::GetText() const {
	std::string tmp = text;
	size_t i = tmp.find("&lt;");
	while (i != std::string::npos) {
		tmp = tmp.replace(i, 4, "<");
		i = tmp.find("&lt;");
	}
	i = tmp.find("&gt;");
	while (i != std::string::npos) {
		tmp = tmp.replace(i, 4, ">");
		i = tmp.find("&gt;");
	}
	return tmp;
}

std::vector<std::string> XMLNode::attrsplit(const std::string &s) const {

	unsigned int lasti=0;
	std::vector<std::string> ret;

	unsigned int i;
	bool cite = false;
	//std::cout << s <<":"<< std::endl;


	for (i=0;i<s.length();i++) {
		if (s.at(i) == '"' || s.at(i) == '\'')
			cite = !cite;

		if (!cite) {
			if (isSplitCharacter(s.at(i))) {
				std::string s1 = s.substr(lasti,i-lasti);
				if (s1.length()) {
					ret.push_back(s1);
				}
				lasti=i+1;
			}
			if (s.at(i)=='=') {
				ret.push_back(std::string("="));
				lasti=i+1;
			}
		}
	}

	if (lasti!=i) {
		std::string s1 = s.substr(lasti,i);
		ret.push_back(s1);
	}

	return ret;
}



bool XMLNode::isSplitCharacter(char c) const {
	switch (c) {
	case '=':
	case '\t':
	case ' ':
	case 10:
	case 13:
		return true;
	}

	return false;

}

bool XMLNode::attrparse(XMLParser *parser,std::vector<std::string> &ss) {


	std::vector<std::string>::iterator it = ss.begin();
	/*	for (;it!=ss.end();it++) {
			std::cout << *it << std::endl;
		}

		it = ss.begin();
	*/
	if (ss.size()) {
		nodename = *it;
		it++;
	} else {
		parser->AddError(XML_NO_ELEMENT_NAME);
		return false;
	}

	while (it!=ss.end()) {
		std::string attval;

		std::string attname = *it;
		it++;

		if (it!=ss.end()) {
			if (std::string("=") != *it) {
				parser->AddError(XML_EQUAL_EXPECTED,attname);
				return false;
			}
		} else {
			parser->AddError(XML_ATTRIBUTE_WITHOUT_ARGUMENTS,attname);
			return false;
		}

		it++;
		if (it!=ss.end())
			attval = *it;
		else {
			parser->AddError(XML_ATTRIBUTE_WITHOUT_ARGUMENTS,attname);
			return false;
		}


		char c1 = attval.at(0);
		char c2 = attval.at(attval.length()-1);

		if (c1==c2 && (c1=='"' || c1=='\'')) {
			int aa= (int)attval.length()-2;
			attribute[attname]= attval.substr(1,aa);
		} else {
			parser->AddError(XML_ATTRIBUTE_NOT_QUOTED,attval);
			return false;
		}

		it++;
	}

	return true;
}


bool XMLNode::doctypeparse(XMLParser *parser, const std::string &a) {

	std::vector<std::string> ss = ToStringArray(a);

	// DOCTYPE
	if (ss[0]=="!DOCTYPE") {

		SetName(std::string(&ss[0][1]));

		if (ss.size()>1) {
			SetAttribute(std::string("id"),ss[1]); // Wurzel-Element
		} else {
			parser->AddError(XML_NO_ROOT);
			return false;
		}

		if (ss.size()>2) {
			int b=2;

			if (ss[2]=="SYSTEM") {
				SetAttribute(std::string("SYSTEM"),ss[3]);
				b=4;
			}
			if (ss[2]=="PUBLIC") {
				SetAttribute(std::string("PUBLIC"),ss[3]);
				SetAttribute(std::string("PUBLIC_URL"),ss[4]);
				b=5;
			}

			// Rest...
			std::string a = Join(std::vector<std::string>(ss.begin()+b,ss.end()));
			std::stringstream r;

			if (a.size()) {
				if (a[0] == '[') {
					SetAttribute(std::string("LOCAL"),a);



				}

			}

		}

	}


	// !ENTITY
	// !ELEMENT
	// !ATTLIST


	return true;
}


std::string replacebrackets(std::string &c) {
	// node.text evtl. für HTML ohne <>><<><

	std::string::iterator it = c.begin();
	int i=0;
	for (;it!=c.end();) {
		if (*it=='<') {
			c.replace(i,1,"&lt;");
			it=c.begin();
			i=0;
		} else if (*it=='>') {
			c.replace(i,1,"&gt;");
			it=c.begin();
			i=0;
		} else {
			it++;
			i++;
		}
	}

	return c;

}


std::string replacequote(std::string &c) {
	// node.text evtl. für HTML ohne <>><<><

	//replace(c.begin(),c.end(),"\"","\\\"");
	// replace(c.begin(),c.end(),'a','b');

	std::string::iterator it = c.begin();
	//int i=0;
	for (;it!=c.end();) {
		if (*it=='"') {
			c.replace(it,it+1,"\\\"");
			it++;
			it++;
			//it=c.begin();
			//i=0;
		} else {
			it++;
			//i++;
		}
	}

	return c;

}

void myReplace(std::string& s, std::string toReplace, std::string replaceWith) {
	size_t first = s.find(toReplace);
	if (first < s.length()) {
		s.replace(first, toReplace.length(), replaceWith);
	}
}

std::string XMLNode::GetAttribute(const std::string &att) const {
	/*std::cout << "Getting " << att << std::endl;
	for (std::map<std::string,std::string>::const_iterator it = attribute.begin(); it != attribute.end(); ++it) {
		std::cout << (*it).first << " , " << (*it).second << std::endl;
	}*/
	std::map<std::string,std::string>::const_iterator it=attribute.find(att);
	if (it!=attribute.end()) {
		return it->second;
	}
	return "";
}

size_t XMLNode::CountAttributes() const {
		return attribute.size();
}

void XMLNode::SetText(const std::string& raw) {
	text = raw;
	size_t i = text.find("<");
	while (i != std::string::npos) {
		text = text.replace(i, 1, "&lt;");
		i = text.find("<");
	}
	i = text.find(">");
	while (i != std::string::npos) {
		text = text.replace(i, 1, "&gt;");
		i = text.find(">");
	}
}

void XMLNode::GetChildValueIfExists(const std::string &name, const std::string Map[], int &value) {

	XMLNode *n = GetFirstChild(name);
	if (n) {
		std::string s = n->GetText();
		for (int i=0;;i++) {
			if (Map[i][0]==0)
				break;
			if (s==Map[i]) {
				value = i;
				break;
			}
		}
	}
}
void XMLNode::GetChildValueIfExists(const std::string &name, const std::vector<char *> &Map, int &value) {

	XMLNode *n = GetFirstChild(name);
	if (n) {
		std::string s = n->GetText();
		std::vector<char*>::const_iterator it = Map.begin();
		int i=0;
		for (;it!=Map.end();it++,i++) {
			if ((*it)[0]==0)
				break;
			if (s==*it) {
				value = i;
				break;
			}
		}
	}
}
void XMLNode::GetChildValueIfExists(const std::string &name, const char Map[3][10], int &value) {


	XMLNode *n = GetFirstChild(name);
	if (n) {
		std::string s = n->GetText();
		for (int i=0;;i++) {
			if (Map[i][0]==0)
				break;
			if (s==Map[i]) {
				value = i;
				break;
			}
		}
	}
}




XMLNode* XMLNode::AddChildIfNotExists(const std::string &name, const std::string &text) {

	XMLNode *n = GetFirstChild(name);
	if (n)
		n->SetText(text);
	else
		n = AddChild(name,text);

	return n;
}

XMLNode* XMLNode::AddChildIfNotExists(const std::string &name) {

	XMLNode *n = GetFirstChild(name);
	if (!n)
		n = AddChild(name);

	return n;
}


XMLNode* XMLNode::Clone() {

	XMLNode *ret = new XMLNode(nodename,text);
	ret->attribute = this->attribute;
	ret->type = this->type;
	std::vector<XMLNode *>::iterator it = child.begin();
	for (;it!=child.end();it++) {
		ret->child.push_back((*it)->Clone());
	}

	return ret;
}


XMLNode *XMLNode::CreateRoot(const std::string &name) {

	XMLNode *root = new XMLNode(-1);

	XMLNode *xml = new XMLNode(1);
	xml->SetName("xml");
	xml->SetAttribute("encoding","UTF-8");
	xml->SetAttribute("standalone","yes");
	xml->SetAttribute("version","1.0");

	XMLNode *doctype = new XMLNode(3);
	doctype->SetName("DOCTYPE");
	doctype->SetAttribute("id", name);
	doctype->SetAttribute("SYSTEM","\"polyearth.dtd\"");

	XMLNode *main = new XMLNode(name);

	root->AddChild(xml);
	root->AddChild(doctype);
	root->AddChild(main);

	return root;
}





std::vector<XMLNode*>::const_iterator XMLNode::begin()  {

	it_tmp = child.begin();
	if (it_tmp!=child.end())
		return it_tmp;

	return child.end();
}

std::vector<XMLNode*>::const_iterator XMLNode::operator++() {
	++it_tmp;
	if (it_tmp!=child.end())
		return it_tmp;

	return child.end();
}

std::vector<XMLNode*>::const_iterator XMLNode::end()  {

	return child.end();

}

bool XMLNode::operator==(const char *name) {
	return GetName()==name;
}

void XMLNode::InsertHere(XMLNode *n) {
	if (it_tmp != child.end()) {
		it_tmp = child.insert(it_tmp, n);
	}
	else {
		child.push_back(n);
		it_tmp = child.end();
	}
}
