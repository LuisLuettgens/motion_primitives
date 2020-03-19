#include "toolstatus.h"

#include "tool.h"
#include "conversion.h"

using namespace std;

namespace tw {

void ToolStatus::Draw(int height) const {

	int aa = (int)statusline.size();
	if (statusline.size()) {

		int yi = -height+(int)statusline.size()*12;

		int bb = statusline.size();
		int ww = 0;
		std::vector<Statusline>::const_iterator it = statusline.begin();
		for (;it!=statusline.end();it++) {
			int a = Tool::font->StringLength((it)->text.c_str());
			if (a>ww)
				ww=a;

			bb = statusline.size();
			if (bb<aa)
				break;
		}

		glDisable(GL_LINE_SMOOTH);
		Tool::drawRect(Point<int>(5,-yi-2-12),Point<int>(ww+25,height-5),-.2f,
		               Tool::FILL2 | Tool::OUTSET);

		it = statusline.begin();
		for (;it!=statusline.end();it++) {

			if ((it)->text.substr(0,3)=="---") {
				Point <int>p1(5+15,    -yi-5);
				Point<int> p2(ww+25-15, -yi-4);
				Tool::drawRect(p1,p2,-.1, Tool::FILL2 | Tool::INSET);

			} else if ((it)->text.substr(0,2)=="  "
			           && (it)->text.substr((it)->text.length()-2,2)=="  ") {

				glColor4fv(Tool::colors[5].Dim(.8).GetData());
				int a = Tool::font->StringLength((it)->text.c_str());
				Tool::font->printString((it)->text.c_str(),(ww-a)/2+15,yi,0);
			} else if ((it)->text.find("|")) {

				glColor4fv(Tool::colors[5].Dim(.8).GetData());

				int i=15;
				vector<string> a = ToStringArray(it->text,"|");
				vector<string>::iterator it2 = a.begin();
				for (;it2!=a.end();it2++) {
					Tool::font->printString(it2->c_str(),i,yi,0);
					i+=5 * it2->size();
				}
			} else {
				glColor4fv(Tool::colors[5].Dim(.8).GetData());
				Tool::font->printString((it)->text.c_str(),15,yi,0);
			}
			yi-=12;
		}

	}

}


void ToolStatus::Timer(int time) {

	if (statusline.size()) {
		if (time-statusline[0].time>5000)
			statusline.erase(statusline.begin());
	}

}

void ToolStatus::Clear() {
	statusline.clear();
}

void ToolStatus::Add(const std::string &s, int time) {

	vector<string> a = ToStringArray(s,"\n");
	vector<string>::iterator it = a.begin();
	for (;it!=a.end();it++) {
		statusline.push_back(Statusline(*it,time));
	}
}

}
