#ifdef WIN32
#include "windows.h"
#endif

#include "conversion.h"
#include "toolmenukey.h"
#include "tool.h"

using namespace std;


//void Message(int id);
//void Message(int id, const void *param);

namespace tw {

ToolMenuKey::ToolMenuKey(const std::string &s, int i) : id(i), ctrl(0), alt(0), shift(0) {

	vector<string> a = ToStringArray(s,"+");

	//m = KMOD_NONE;
	k = SDLK_CLEAR;

	for (unsigned int i=0;i<a.size();i++) {
		if (a[i] == "Ctrl") {
			ctrl = 1; //m = (SDLMod) KMOD_CTRL;
			continue;
		}
		if (a[i] == "Alt") {
			alt = 1; //m = (SDLMod) KMOD_CTRL;
			continue;
		}
		if (a[i] == "Shift") {
			shift = 1; //m = (SDLMod) KMOD_CTRL;
			continue;
		}
		if (a[i][0] == 'F' && a[i][1]>='0' && a[i][1]<='9') {
			k = (SDL_Keycode) (SDLK_F1 + ToInt(string(a[i],1,string::npos)) - 1);
			continue;
		}
		k = Tool::toKey(a[i][0]);
	}

	//cout << m << " " << k  << " " << id << endl;
}

bool ToolMenuKey::Check(SDL_Keycode key, Uint16 mod) {
	if (key == k) {
		int c = (ctrl==0?1:0) ^ ((mod & KMOD_CTRL)?1:0);
		int s = (shift==0?1:0) ^ ((mod & KMOD_SHIFT)?1:0);
		int a = (alt==0?1:0) ^ ((mod & KMOD_ALT)?1:0);

		if (c && s && a) {
			//if (((mod & m) || m==KMOD_NONE ) && (key == k)) {

			if (id) {
				Message(id)
				;
				return true;
			}
		}
	}
	return false;
}

}
