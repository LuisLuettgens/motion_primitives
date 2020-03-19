#ifdef WIN32
#include <windows.h>
#endif

#include "textout.h"

TextOutputType_e standard_textoutputtype = CONSOLE;


std::ostream& beginfile ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	deb->boldfont_m=false;
	deb->colorfont_m=false;

	if (deb->GetTextOutputType()==CPP) {
		os << "std::string xmlfile = \"";
	}
	else if (deb->GetTextOutputType()==RTF)  {
		os << "{\\rtf1\\ansi\\ansicpg1252\\deff0\\deflang1031";
		os << "{\\fonttbl{\\f0\\fswiss\\fprq2\\fcharset0 Microsoft Sans Serif;}{\\f1\\fmodern\\fprq1\\fcharset0 Courier New;}}";
		os << "{\\colortbl ;\\red255\\green0\\blue0;\\red16\\green78\\blue139;\\red110\\green139\\blue61;\\red192\\green192\\blue192;\\red191\\green11\\blue11;\\red0\\green0\\blue0;}";
		os << "\\viewkind4\\uc1\\pard\\tx1500\\f0\\fs16 ";
	} else if (deb->GetTextOutputType()==HTML) {
		os << //"<html>\n\t<head>\n\t\t<title>debug</title>\n\t</head>\n\t<body>\n\t\t<font size=2>
		"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n"
		<< "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"de\" lang=\"de\">\n"
		<< "\t<head>\n"
		<< "\t\t<title>Sitemap</title>\n"
		<< "\t\t<style type=\"text/css\">\n"
		<< "\t\t\tbody {\n"
		<< "\t\t\tcolor: #000000;\n"
		<< "\t\t\tfont-family: Courier, monospace;\n"
		<< "\t\t\tfont-size: 1.em;\n"
		<< "\t\t\t}\n"
		<< "\t\t\tdiv.sub{\n"
		<< "\t\t\tmargin-left:3em;\n"
		//		<< "\t\t\tpadding:2px;\n"
		//		<< "\t\t\tborder-bottom:1px solid;\n"
		//		<< "\t\t\tborder-left:1px solid;\n"
		//		<< "\t\t\tborder-color:#ff0000;\n"
		//		<< "\t\t\tmin-height:20px;\n"
		<< "\t\t\t}\n"
		<< "\t\t\t.col1{\n"
		<< "\t\t\tcolor:#ff0000;\n"
		<< "\t\t\t}\n"
		<< "\t\t\t.col2{\n"
		<< "\t\t\tcolor:#104e8b;\n"
		<< "\t\t\t}\n"
		<< "\t\t\t.col3{\n"
		<< "\t\t\tcolor:#6e8b3d;\n"
		<< "\t\t\t}\n"
		<< "\t\t\t.col4{\n"
		<< "\t\t\tcolor:#cccccc;\n"
		<< "\t\t\t}\n"
		<< "\t\t\t.col5{\n"
		<< "\t\t\tcolor:#bf0b0b;\n"
		<< "\t\t\t}\n"
		<< "\t\t</style>\n"
		<< "\t</head>\n"
		<< "\t<body>";

	}
	return os;
}

std::ostream& endfile ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==CPP) {
		os << "\";\n";
	}
	else if (deb->GetTextOutputType()==RTF)  {
		os << "\\cf0 \\par }";
	} else if (deb->GetTextOutputType()==HTML) {
		os << "\n\t</body>\n</html>";

	}
	return os;
}

std::ostream& endline ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==CPP) {
		os << "\\n\" \\\n\t\"";
	}
	else if (deb->GetTextOutputType()==RTF)  {
		os << "\\par ";
	} else if (deb->GetTextOutputType()==HTML)  {
		os << "<br>\n";
	} else {
		os << "\n";
	}
	return os;
}


std::ostream& monospace ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==CPP) {
		os << "";
	}
	else if (deb->GetTextOutputType()==RTF)  {
		os << "\\f1 ";
	} else if (deb->GetTextOutputType()==HTML)  {
		os << "";
	} else {
		os << "";
	}
	return os;
}


std::ostream& proportional ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==CPP) {
		os << "";
	}
	else if (deb->GetTextOutputType()==RTF)  {
		os << "\\f0 ";
	} else if (deb->GetTextOutputType()==HTML)  {
		os << "";
	} else {
		os << "";
	}
	return os;
}

std::ostream& quote ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==CPP)  {
		os << "\\\"";
	} else {
		os << "\"";
	}
	return os;
}
std::ostream& lt ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==RTF)  {
		os << "<";
	} else if (deb->GetTextOutputType()==HTML)  {
		os << "&lt;";
	} else {
		os << "<";
	}
	return os;
}
std::ostream& gt ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==RTF)  {
		os << ">";
	} else if (deb->GetTextOutputType()==HTML)  {
		os << "&gt;";
	} else {
		os << ">";
	}
	return os;
}

std::ostream& tab ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);
	if (deb->GetTextOutputType() == RTF) {
		for (int i=0;i<deb->ident_m;i++)
			os << "\\tab ";
	} else if (deb->GetTextOutputType()==HTML)  {}
	else {
		for (int i=0;i<deb->ident_m;i++)
			os << "\t";
	}
	return os;
}

std::ostream& addtab ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);
	deb->ident_m++;
	if (deb->GetTextOutputType()==RTF)  {}
	else if (deb->GetTextOutputType()==HTML)  {
		if (deb->ident_m>0)
			os << "<div class=\"sub\">";
	} else {}
	return os;
}

std::ostream& removetab ( std::ostream& os ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);
	deb->ident_m--;
	if (deb->GetTextOutputType()==RTF)  {}
	else if (deb->GetTextOutputType()==HTML)  {
		if (deb->ident_m>=0)
			os << "</div>";
	} else {}
	return os;
}


std::ostream& textcolor ( std::ostream& os, int i ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==RTF)  {
		switch(i) {
		case 0: // normal
		case 1: // rot
		case 2: // blau
		case 3: // grün
		case 4: // grau
		case 5: // dark red
			os << "\\cf" << i << " ";
			break;
		}
	} else if (deb->GetTextOutputType()==HTML)  {
		switch(i) {
		case 0:
			//if (deb->boldfont_m)
			//	os << "</b>";
			if (deb->colorfont_m)
				os << "</span>";
			//	deb->boldfont_m=false;
			deb->colorfont_m=false;
			break;
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
			if (deb->colorfont_m)
				os << "</span>";

			os << "<span class=\"col" << i << "\">";
			deb->colorfont_m=true;
			break;

			/*case 4:
				os << "<font color=#ff0000>";
				deb->colorfont_m=true;
				break;
			case 5:
				os << "<font color=#6E8B3D>";
				deb->colorfont_m=true;
				break;*/
		}
	} else if (deb->GetTextOutputType()==CONSOLE) {

#ifdef WIN32

		HANDLE hConsole0;

		int cols[] = {
			FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED,
			FOREGROUND_RED | FOREGROUND_INTENSITY,
			FOREGROUND_GREEN | FOREGROUND_INTENSITY,
			FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED, //| BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED,
			FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_INTENSITY,
			FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_INTENSITY
		};

		hConsole0 = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hConsole0, cols[i]);

		//std::cout << i << std::endl;



		/*switch(i) {

		case 0:
			if (deb->colorfont_m) {
				SetConsoleTextAttribute(hConsole0, 7);
				deb->colorfont_m=true;
			}
			break;
		case 1:
			if (deb->colorfont_m) {
				SetConsoleTextAttribute(hConsole0, 12);
				deb->colorfont_m=true;
			}
			break;
		case 2:
			if (deb->colorfont_m) {
				SetConsoleTextAttribute(hConsole0, 10);
				deb->colorfont_m=true;
			}
			break;
		case 3:
			if (deb->colorfont_m) {
				SetConsoleTextAttribute(hConsole0, 9);
				deb->colorfont_m=true;
			}
			break;
		case 4:
			if (deb->colorfont_m) {
				SetConsoleTextAttribute(hConsole0, 15);
				deb->colorfont_m=true;
			}
			break;
		case 5:
			if (deb->colorfont_m) {
				SetConsoleTextAttribute(hConsole0, 14);
				deb->colorfont_m=true;
			}
			break;

		//SetColor(i);
		}*/

		//std::cout << "::" << i << "::" << std::endl;



#else
		switch(i) {
		case 0:
			os << "\033[0m";
			deb->colorfont_m=true;
			break;
		case 1:
			os << "\033[31m";
			deb->colorfont_m=true;
			break;
		case 2:
			os << "\033[32m";
			deb->colorfont_m=true;
			break;
		case 3:
			os << "\033[34m";
			deb->colorfont_m=true;
			break;
		case 4:
			os << "\033[37m";
			deb->colorfont_m=true;
			break;
		case 5:
			os << "\033[33m";
			deb->colorfont_m=true;
			break;
		}

#endif

	}
	return os;
}


std::ostream& textsize ( std::ostream& os, int i ) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->GetTextOutputType()==RTF)  {
		switch(i) {
		case 0:
			os << "\\fs16 ";
			break;	// normal
		case 1:
			os << "\\fs24 ";
			break;	// large
		case 2:
			os << "\\fs32 ";
			break;	// large

		}
	} else if (deb->GetTextOutputType()==HTML)  {
		/*switch(i) {
		case 0:
			if (deb->boldfont_m) os << "</b>";
			if (deb->colorfont_m) os << "</font>";
			deb->boldfont_m=false;
			deb->colorfont_m=false;
			break;
		case 1:
			os << "<font color=#ff0000><b>"; deb->colorfont_m=true; deb->boldfont_m=true; break;
		case 2:
			os << "<font color=#104E8B>"; deb->colorfont_m=true; break;
		case 3:
			os << "<b>"; deb->boldfont_m=true; break;
		case 4:
			os << "<font color=#ff0000>"; deb->colorfont_m=true; break;
		case 5:
			os << "<font color=#6E8B3D>"; deb->colorfont_m=true; break;
		}*/
	}
	return os;
}



std::ostream& bold ( std::ostream& os) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	deb->boldfont_m = true;
	if (deb->GetTextOutputType()==RTF)  {
		os << "\\b ";
	} else if (deb->GetTextOutputType()==HTML)  {
		os << "<b>";
	} else if (deb->GetTextOutputType()==CONSOLE)  {
		os << "\033[4m";
	}
	return os;
}

std::ostream& oblique ( std::ostream& os) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	deb->obliquefont_m = true;
	if (deb->GetTextOutputType()==RTF)  {
		os << "\\i ";
	} else if (deb->GetTextOutputType()==HTML)  {
		os << "<i>";
	} else if (deb->GetTextOutputType()==CONSOLE)  {
		os << "\033[7m";
	}

	return os;
}

std::ostream& normal ( std::ostream& os) {

	DebugStream* deb = dynamic_cast<DebugStream*>(&os);

	if (deb->boldfont_m) {
		if (deb->GetTextOutputType()==RTF)  {
			os << "\\b0 ";
		} else if (deb->GetTextOutputType()==HTML)  {
			os << "</b>";
		} else if (deb->GetTextOutputType()==CONSOLE)  {
			os << "\033[24m";
		}

		deb->boldfont_m = false;
	}
	if (deb->obliquefont_m) {
		if (deb->GetTextOutputType()==RTF)  {
			os << "\\i0 ";
		} else if (deb->GetTextOutputType()==HTML)  {
			os << "</i>";
		} else if (deb->GetTextOutputType()==CONSOLE)  {
			os << "\033[27m";
		}

		deb->obliquefont_m = false;
	}
	return os;
}

textcolor_int textcolor( int l ) {
	return textcolor_int ( textcolor, l );
}

std::ostream& operator<<(std::ostream& s, const textcolor_int & sm) {
	(*sm._fp)(s,sm._tp);
	return s;
}


textsize_int textsize( int l ) {
	return textsize_int ( textsize, l );
}

std::ostream& operator<<(std::ostream& s, const textsize_int & sm) {
	(*sm._fp)(s,sm._tp);
	return s;
}



