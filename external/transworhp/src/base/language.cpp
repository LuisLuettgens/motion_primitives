#include "language.h"

#include "conversion.h"
#include "../core/twstatus.h"

using namespace std;

string Language::empty="";

void Language::Init(XMLNode *xml) {

	id = xml->GetAttribute("id");
	text = xml->GetAttribute("text");


	XMLNode *n = xml->GetFirstChild("TEXT");
	while (n) {
		int a = ToInt(n->GetAttribute("id"));
		string s = n->GetText();

		data[a] = s;

		n = xml->GetNextChild("TEXT");
	}

	MyStatus("Language",
	         string("Loading Language ") + id + ": " + to_string(data.size()) + " entries",
	         Status::WARN);
}

const std::string &Language::GetText(int a) const {

	map<int,string>::const_iterator it = data.find(a);
	if (it!=data.end()) {return it->second;}
	else return empty;

}


void LanguageManager::Init(XMLNode *xml) {

	XMLNode *n = xml->GetFirstChild("LANGUAGE");
	maxlan=0;
	while (n) {
		l[maxlan].Init(n);
		maxlan++;
		n = xml->GetNextChild("LANGUAGE");

	}
}

const std::string &LanguageManager::GetText(int a) const {
	return l[curlan].GetText(a);
}

const std::string &LanguageManager::GetName(int a) const {
	return l[a].GetName();

}
