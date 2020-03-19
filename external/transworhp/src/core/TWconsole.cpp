#include "TWconsole.h"

#include "xmlio.h"

#include <string>

#ifdef WIN32
#include <windows.h>
#endif

TWconsole::TWconsole() : active(0), fw(6), fh(8), bw(132), bh(2000) {}

#ifndef TRANSWORHP_GRAPHICS
void TWconsole::Init() {}
void TWconsole::Resize(int /*bh_*/) {}
#else

#ifdef WIN32

void TWconsole::Init() {
	
	HANDLE outcon = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_FONT_INFOEX font;
	font.cbSize=sizeof(CONSOLE_FONT_INFOEX);
	GetCurrentConsoleFontEx(outcon, false, &font);
	font.dwFontSize.X = fw;
	font.dwFontSize.Y = fh;
	SetCurrentConsoleFontEx(outcon, false, &font);

	COORD dwSize;
	dwSize.X = bw;
	dwSize.Y = bh;
	SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), dwSize);
}

void TWconsole::Resize(int bh_) {
	HWND hWnd;
	RECT WinRect;

	std::string s("TransWORHP Log");
	SetConsoleTitleA(s.c_str());
	hWnd = FindWindowA(NULL,
	                   s.c_str());

	GetWindowRect(hWnd,
	              &WinRect);

	MoveWindow(hWnd,
	           WinRect.left,
	           WinRect.top,
	           bw*fw+40,
	           bh_*fh,
	           TRUE);

}
#else
void TWconsole::Init() {}
void TWconsole::Resize(int /*bh_*/) {}
#endif

#endif

void TWconsole::ParseXML(XMLNode *n, int *countParams, int *setParams) {

	if (n) {
		active = 1;

		(*countParams) += 2;
		XMLNode *xml = n->GetFirstChild("FONT");
		if (xml) {
			std::string s( xml->GetAttribute("width") );
			if (!s.empty()) {
				fw = std::stoi(s);
				(*setParams)++;
			}
			s = xml->GetAttribute("height");
			if (!s.empty()) {
				fh = stoi(s);
				(*setParams)++;
			}
		}

		(*countParams) += 2;
		std::string s( n->GetAttribute("width") );
		if (!s.empty()) {
			bw = std::stoi(s);
			(*setParams)++;
		}
		s = n->GetAttribute("height");
		if (!s.empty()) {
			bh = std::stoi(s);
			(*setParams)++;
		}

		(*countParams) += 1;
		s = n->GetAttribute("color");
		if (s == "CONSOLE") {
			standard_textoutputtype = CONSOLE;
			(*setParams)++;
		} else if (s == "ASCII") {
			standard_textoutputtype = ASCII;
			(*setParams)++;
		} else if (s == "HTML") {
			standard_textoutputtype = HTML;
			(*setParams)++;
		} else if (s == "RTF") {
			standard_textoutputtype = RTF;
			(*setParams)++;
		} else if (s == "SILENT") {
			standard_textoutputtype = SILENT;
			(*setParams)++;
		}

		Init();
	}
}

int TWconsole::Active() const {
	return active;
}
