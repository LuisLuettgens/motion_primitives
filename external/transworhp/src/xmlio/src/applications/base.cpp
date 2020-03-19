#include "base.h"

#include <cstring>

TextOutputType_e GetTextOutputType(char *s) {

	TextOutputType_e tot = ASCII;

	if (!strcmp(s, "CONSOLE"))
		tot = CONSOLE;
	else if (!strcmp(s, "HTML"))
		tot = HTML;
	else if (!strcmp(s, "RTF"))
		tot = RTF;
	else if (!strcmp(s, "CPP"))
		tot = CPP;

	return tot;
}

