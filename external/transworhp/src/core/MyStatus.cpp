#include "twstatus.h"

#include "conversion.h"
#include "textout.h"

#include <iomanip>
#include <vector>

std::ostream *transworhpOutput = &std::cout;
tw_output newoutputfunction = nullptr;


void SetMyStatusStream(std::ostream &os) {
	transworhpOutput = &os;
}

Status operator|(Status a, Status b) {
	return static_cast<Status>(static_cast<int>(a) | static_cast<int>(b));
}

bool operator&(Status a, Status b) {
	return static_cast<bool>(static_cast<int>(a) & static_cast<int>(b));
}

void MyStatusLine(const std::string& tag, const std::string& msg, Status flag) {

	if (standard_textoutputtype == SILENT && flag != Status::ERR) {
		return;
	}

	if (newoutputfunction) {
		(*newoutputfunction) (tag.c_str(), msg.c_str(), flag);
		return;
	}

	if (!(flag & Status::CONT_PREV)) {
		DebugStream d1(standard_textoutputtype);
		d1 << textcolor(4) << std::setw(12) << tag;
		*transworhpOutput << d1.GetString() << " ";

		DebugStream d2(standard_textoutputtype);
		d2 <<  textcolor(0)<<  "| ";
		*transworhpOutput << d2.GetString();
	}
	// textcolor CODE:
	// 0 = normal (schwarz bzw. weiß)
	// 1 = rot
	// 2 = gruen
	// 3 = blau
	// 4 = grau
	// 5 = gelb

	DebugStream d3(standard_textoutputtype);
	if (flag == Status::NORMAL) {
		d3 << textcolor(0);
	}
	else if (flag & Status::RED && flag & Status::GREEN) {
		d3 << textcolor(5);
	}
	else if (flag & Status::RED) {
		d3 << textcolor(1);
	}
	else if (flag & Status::GREEN) {
		d3 << textcolor(2);
	}
	else if (flag & Status::BLUE) {
		d3 << textcolor(3);
	}

	d3 << msg;
	*transworhpOutput << d3.GetString();

	if (!(flag & Status::CONT_NEXT)) {
		DebugStream d4(standard_textoutputtype);
		d4 <<  textcolor(0) << std::endl;
		*transworhpOutput << d4.GetString();
	} else {
		transworhpOutput->flush();
	}
}

void MyStatus(const std::string& tag, const std::string& message, Status flag) {

	const std::vector<std::string> s = ToStringArray(message, "\n");

	if (!s.empty()) {
		auto it = s.cbegin();

		MyStatusLine(tag, *it, flag);
		it++;

		for (; it != s.end(); it++) {
			MyStatusLine("|", *it, flag);
		}
	}
}

void DllExport SetMyOutputFunction(tw_output f) {
	newoutputfunction = f;
}
