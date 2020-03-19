#include "font.h"

#include "conversion.h"

#ifdef WIN32
#include "windows.h"
#endif

const unsigned char FFont::psi = 130;
const unsigned char FFont::phi = 131;
const unsigned char FFont::lambda = 132;
const unsigned char FFont::alpha = 133;
const unsigned char FFont::beta = 134;
const unsigned char FFont::gamma = 135;
const unsigned char FFont::delta = 136;
const unsigned char FFont::epsilon = 137;
const unsigned char FFont::Omega = 138;
const unsigned char FFont::omega = 139;

const unsigned char FFont::UP = 16;
const unsigned char FFont::DOWN = 17;
const unsigned char FFont::UPUP = 18;
const unsigned char FFont::DOWNDOWN = 19;
const unsigned char FFont::STEP = 20;
const unsigned char FFont::BACK = 21;

const unsigned char FFont::CHECK = 15;
const unsigned char FFont::UNCHECK = 14;
const unsigned char FFont::GORIGHT = 13;


const unsigned char FFont::DOT = 144; // \220


const FFont::letter* FFont::Key(unsigned char s) const {

	std::vector<letter>::const_iterator it = letters.begin();
	for (;it!=letters.end();it++) {
		if (it->character==s) {
			const letter &l = *it;
			return &l;
		}
	}

	return nullptr;
}


int FFont::StringLength(const char *s) const {

	if (s==0) return 0;

	int ret = 0;
	for (int i = 0; s[i]; i++) {
		const letter* l = Key(s[i]);
		if (l) ret += l->width;
	}
	return ret;
}


/* ---------------- LETTER ------------- */

FFont::letter::letter(const std::string &text) {

	height = 0;

	std::vector<std::string> split = ToStringArray(text,std::string(","));

	if (split.size()>4) {
		std::string a = eatwhitespace(split[0]);
		character = a[1];
		width = ToInt(split[1]);

		for (int i=0;i<10;i++) {
			raster[i] = ToInt(split[2+i]);
		}
	}
}


FFont::letter::letter(XMLNode *node) {

	height = 0;

	//character = node->GetAttribute(std::string("char"))[0];
	character = ToInt(node->GetAttribute(std::string("char")));
	width = ToInt(node->GetAttribute(std::string("width")));
	height = ToInt(node->GetAttribute(std::string("height")));

	std::vector<int> split = ToIntArray(node->GetText(), std::string(","));

	//std::string a = eatwhitespace(split[0]);

	for (int i = 0; i < 10; i++) {
		raster[9-i] = split[i];
	}

}


std::ostream &operator<<(std::ostream &os, const FFont::letter &c) {

	for (int i=9;i>=0;i--) {
		for (int j=128;j;j/=2) {
			if (j & c.raster[i])
				os << "#";
			else
				os << "_";
		}
		os << std::endl;
	}
	return os;

}
