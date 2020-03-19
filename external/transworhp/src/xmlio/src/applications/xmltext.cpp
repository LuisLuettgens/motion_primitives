#include <iostream>
#include "xmlio.h"
#include "textout.h"
//#include <cstring>

#include "base.h"


std::string a= "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone = 'yes' ?>\n" \
"<?xml-stylesheet type=\"text/xsl\" href=\"kochbuch.xsl\"?>\n" \
"<!DOCTYPE SIMULATION SYSTEM \"simulation.dtd\" >\n"\
"<!-- Start -->\n"\
"<SIMULATION>"\
"	<RAUM id=\"Halle1\">"\
"		<DIMENSION x=\"-10 40\" y=\"-.01 15\" z=\"-20 30\"/>"\
"	</RAUM> "\
"<!-- sdf ><>> <<<<<<<<>>>>><<<<kommentar -->"\
"	<TRANSFASTERTYPE id=\"Modell 1\" lam=\"hubtisch\">"\
"	<!-- . -->"\
"		<WIDTH type=\"katze\">3.2</WIDTH>"\
"		<MOVEMENT   type=\"leer\">8 4 3</MOVEMENT>"\
"		<MOVEMENT    type=\"voll\"   >1 2 3</MOVEMENT>"\
"	</TRANSFASTERTYPE>"\
"	<TRANSFASTER id=\"Homer\" schieneid=\"oben\" connectionid=\"fm1\" "\
"	typeid=\"Modell 1\""\
"	"\
"	color=\"10,1,125\">"\
"		<SEILDIFF>	1.5		</SEILDIFF>"\
"		<HOEHE>1.5</HOEHE>"\
"	</TRANSFASTER>"\
"</SIMULATION>";


int main(int argc, char *argv[]) {

	if (argc>1) {

		XMLParser p;

		TextOutputType_e tot = ASCII;
		if (argc > 1) {
			tot = GetTextOutputType(argv[1]);
		}

		std::cout << a;

		XMLNode *node = p.ParseString(a); 	// Parses the file and returns root node.

		// std::cout << beginfile;		// Wrong! Don't write output directly to std::cout.

		DebugStream d(tot); 			// Open DebugStream.
		d << beginfile;
		if (node) {
			d << *node; 			// Write whole file.
			//d << *p.GetDocument(); 	// Write only main part.

		} else {
			p.GetError(d);			// Show all errors.
			//d << *node;			// No! We don't have this node,
			//d << *p.GetDocument(); 	// and no document at all.
		}
		d << endfile;

		std::cout << d.GetString();		// Finally write DebugStream anywhere you want.

	} else {
		std::cout << "Usage: " << argv[0] << " xml-file [type]" << std::endl;
		std::cout << "with type = [ASCII | HTML | RTF | CONSOLE] " << std::endl;
	}
	return 0;
}

