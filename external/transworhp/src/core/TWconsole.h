#pragma once

class XMLNode;

/** Manipulation der Win-Console */
class TWconsole {

public:
	TWconsole();

	/** XML lesen.
	 * @param xmlmain XML-Knote
	 * @param countParams zum Parameter zeahlen: Gesamtzahl aller Parameter
	 * @param setParams #gesetzter Parameter
	 */
	void ParseXML(XMLNode *xmlmain, int *countParams=nullptr, int *setParams=nullptr);

	void Init();
	void Resize(int bh_);

	int Active() const;

private:
	int active, fw, fh, bw, bh;
};
