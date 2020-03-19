#include <iostream>
#include "xmlio.h"
#include "textout.h"
//#include <cstring>

#include "base.h"


int main(int argc, char *argv[]) {

	if (argc>1) {

		XMLParser p;

		std::string filename(argv[1]);

		TextOutputType_e tot = ASCII;
		if (argc > 2) { 
			tot = GetTextOutputType(argv[2]); 
		}
		

		XMLNode *node = p.Parse(filename); 	// Parses the file and returns root node.

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

	} 
	else {
		std::cout << "Usage: " << argv[0] << " xml-file [type]" << std::endl;
		std::cout << "with type = [ASCII | HTML | RTF | CONSOLE] " << std::endl;
	}
}

