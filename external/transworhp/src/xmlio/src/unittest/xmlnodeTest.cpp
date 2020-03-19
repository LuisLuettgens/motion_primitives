#include "../xmlio/xmlio.h"

#include <boost/test/unit_test.hpp>

#include <math.h>
#include <algorithm>
#include <string>
#include <utility>


std::string a = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone = 'yes' ?>\n" \
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

BOOST_AUTO_TEST_CASE(xmlnodeTest) {

	XMLNode* node = nullptr;
	XMLNode* clone = nullptr;

	{
		XMLParser p;
		node = p.ParseString(a);

		int i = 0;
		XMLNode *xml = node->GetFirstChild();
		while (xml) {
			i++;
			xml = node->GetNextChild();
		}
		BOOST_CHECK_EQUAL(i, 5);

		clone = node->Clone();

		// node von p wird gelöscht
	}

	node = nullptr;


	int i = 0;
	XMLNode *xml = clone->GetFirstChild();
	while (xml) {
		i++;
		xml = clone->GetNextChild();
	}
	BOOST_CHECK_EQUAL(i, 5);

	delete clone;
}

BOOST_AUTO_TEST_CASE(addChildAfterTest) {

	XMLNode* node = nullptr;
	int oldSize = 0;

	XMLParser p;
	node = p.ParseString(a);



	XMLNode *root = node->GetFirstChild("SIMULATION");

	XMLNode *xml = root->GetFirstChild();
	while (xml) {
		oldSize++;
		xml = root->GetNextChild();
	}
	BOOST_CHECK_EQUAL(oldSize, 4);

	xml = root->GetFirstChild("RAUM");
	XMLNode *ins = new XMLNode("RAUM");
	ins->SetAttribute("id", "inserted");
	auto ret = root->AddChildAfter(ins, xml);
	// Iterator was resetted, thus GetFirstChild again
	xml = root->GetFirstChild("RAUM");
	// Check if it is the next AddChildAfter
	xml = root->GetNextChild("RAUM");
	BOOST_CHECK_EQUAL(xml->GetAttribute("id"), ret->GetAttribute("id"));
	// Check if new size is correct
	int newSize = 0;
	xml = root->GetFirstChild();
	while (xml) {
		newSize++;
		xml = root->GetNextChild();
	}
	BOOST_CHECK_EQUAL(newSize, 5);
}
