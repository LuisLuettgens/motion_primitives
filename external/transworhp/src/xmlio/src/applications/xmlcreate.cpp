#include <iostream>
#include "xmlio.h"
#include "textout.h"
#include <cstring>

#include "base.h"


int main(int argc, char *argv[]) {

	TextOutputType_e tot = ASCII;
	if (argc > 1) {
		tot = GetTextOutputType(argv[1]);
	}

	XMLNode *node = XMLNode::CreateRoot("ADRESSEN");
	XMLNode *n = node->GetTypedChild(0);


	n->AddChild("PERSON")->SetAttribute("name", "Matthias");
	n->AddComment("HALLO");
	n->AddChild("PERSON", "Matthias");

	DebugStream d(tot); 			// Open DebugStream.
	d << beginfile;
	if (node) {
		d << *node; 			// Write whole file.
		//d << *p.GetDocument(); 	// Write only main part.
	}

	d << endfile;

	std::cout << d.GetString();		// Finally write DebugStream anywhere you want.

}

