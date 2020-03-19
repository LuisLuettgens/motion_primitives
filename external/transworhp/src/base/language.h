#ifndef language_h
#define language_h
#include "defines.h"
#include <map>
#include <string>
#include "xmlio.h"


class Language {
public:
	Language() {}

	void Init(XMLNode *xml);
	const std::string &GetText(int a) const;
	const std::string &GetName() const {return text;}
private:
	std::map<int,std::string> data;
	std::string id;
	std::string text;

	static std::string empty;
};


class LanguageManager {
public:
	LanguageManager() : curlan(0), maxlan(0) {}

	void Init(XMLNode *xml);
	const std::string &GetText(int a) const;

	void SetLanguage(int i) {curlan=i;}

const std::string &GetName(int a) const;

int curlan;
int maxlan;
private:
	
	Language l[10];

};

#endif


