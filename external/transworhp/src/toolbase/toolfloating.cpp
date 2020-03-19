#ifdef WIN32
#include "windows.h"
#endif

#include "toolfloating.h"
#include "toolitem.h"

ToolFloating::ToolFloating() : item(0), moving(0),rect0(0,0,0,0) {}


ToolFloating::~ToolFloating() {

	delete item;
}


bool ToolFloating::MouseMotion(SDL_MouseMotionEvent &m) {

	if (moving==1) {

		if (m.state) {
			item->rect.x = rect0.x + m.x-mouse.x;
			item->rect.y = rect0.y + m.y-mouse.y;
			return true;
		}
		else moving = 0;

	}


	if (moving==2) {

		if (m.state) {
			item->rect.y = rect0.y + m.y-mouse.y;
			item->rect.height = rect0.height - (m.y-mouse.y);
			if (item->rect.height<30) {
				item->rect.y = rect0.y + rect0.height - 30;
				item->rect.height=30;
			}
			return true;
		}
		else moving = 0;

	}

	if (moving==3) {

		if (m.state) {
			//item->rect.y = rect0.y + m.y-mouse.y;
			item->rect.height = rect0.height + (m.y-mouse.y);
			if (item->rect.height<30) {
				item->rect.height=30;
			}
			return true;
		}
		else moving = 0;

	}

	return false;
}

bool ToolFloating::MouseButton(SDL_MouseButtonEvent &m) {

 int dist = 5;
// oben

	if (m.x>item->rect.x-dist
			&& m.x<item->rect.x+item->rect.width+dist
			&& m.y>item->rect.y-dist
			&& m.y<item->rect.y+item->rect.height+dist) {


		if (m.type == SDL_MOUSEBUTTONDOWN) {
// in den Ecken
			if ((m.y<item->rect.y+dist || m.y>item->rect.y+item->rect.height-dist)
					&& (m.x<item->rect.x+dist || m.x>item->rect.x+item->rect.width-dist)) {

				mouse = Point<int>(m.x,m.y);
				rect0 = item->rect;
				moving=1;
				return true;
			}

			else if (m.y<item->rect.y+dist) {

				mouse = Point<int>(m.x,m.y);
				rect0 = item->rect;
				moving=2;
				return true;
			}


			else if (m.y>item->rect.y+item->rect.height-dist) {

				mouse = Point<int>(m.x,m.y);
				rect0 = item->rect;
				moving=3;
				return true;
			}
		}
		else {
			moving=0;
		}



		item->MouseClick(m.button,m.type,Point<int>(m.x,m.y));

		return true;
	}

	return false;
}


void ToolFloating::Draw(int HT, int WD, const Point<int> &mouse) {

	item->Draw(1,mouse);

}
