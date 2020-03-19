#include "sdlcursor.h"

namespace tw {

const
#include "arrow1.xpm"
const
#include "crosshair.xpm"
//const
//#include "info.xpm"
const
#include "arrow_right.xpm"
const
#include "arrow_left.xpm"
const
#include "arrow_up.xpm"
const
#include "arrow_down.xpm"
const
#include "arrow_left_up.xpm"
const
#include "arrow_left_down.xpm"
const
#include "arrow_right_up.xpm"
const
#include "arrow_right_down.xpm"


SDLCursor cursor[24];

void SDLCursor::Init(const char *image[], const Point<int> &hot) {
	int i, row, col;
	Uint8 data[4*32];
	Uint8 mask[4*32];

	i = -1;
	for ( row=0; row<32; ++row ) {
		for ( col=0; col<32; ++col ) {
			if ( col % 8 ) {
				data[i] <<= 1;
				mask[i] <<= 1;
			} else {
				++i;
				data[i] = mask[i] = 0;
			}
			switch (image[4+row][col]) {
			case '.':
				data[i] |= 0x01;
				mask[i] |= 0x01;
				break;
			case '+':
				mask[i] |= 0x01;
				break;
// 	case '@':
//           data[i] |= 0x01;
//           break;
			case ' ':
				break;
			}
		}
	}
// sscanf(image[4+row], "%d,%d", &hot_x, &hot_y);
	sdl_cursor = SDL_CreateCursor(data, mask, 32, 32, hot.x, hot.y);
}

void SDLCursor::Clean() {
//std::	cout << "freecursor" << std::endl;

//	SDL_FreeCursor(sdl_cursor);
}

void SDLCursor::Set() {
	SDL_SetCursor(sdl_cursor);
}



//extern SDLCursor cursor[24];

/*
int cursormodes[] = {GLUT_CURSOR_CROSSHAIR,
                     GLUT_CURSOR_INFO,
                     GLUT_CURSOR_LEFT_SIDE,
                     0,
                     GLUT_CURSOR_RIGHT_SIDE,
                     0,
                     0,
                     0,
                     GLUT_CURSOR_TOP_SIDE,
                     0,
                     GLUT_CURSOR_TOP_LEFT_CORNER,
                     0,
                     GLUT_CURSOR_TOP_RIGHT_CORNER,
                     0,
                     0,
                     0,
                     GLUT_CURSOR_BOTTOM_SIDE,
                     0,
                     GLUT_CURSOR_BOTTOM_LEFT_CORNER,
                     0,
                     GLUT_CURSOR_BOTTOM_RIGHT_CORNER
                    };
*/

void SDLCursor::CreateCursors() {

	cursor[0].Init(arrow1_xpm,   Point<int>(0,0));
	cursor[1].Init(crosshair_xpm,        Point<int>(15,15));
	cursor[2].Init(arrow_left_xpm,  Point<int>(0,7));
	//cursor[3].Init(arrow_right_xpm, Point<int>(31,16));
	cursor[4].Init(arrow_right_xpm, Point<int>(15,7));

	//  cursor[5].Init(arrow_right_xpm,Point<int>(31,16));
	//  cursor[6].Init(arrow_right_xpm,Point<int>(31,16));
	//  cursor[7].Init(arrow_right_xpm,Point<int>(31,16));
	cursor[8].Init(arrow_up_xpm,Point<int>(7,0));
	//cursor[9].Init(arrow_down_xpm,Point<int>(31,16));
	cursor[10].Init(arrow_left_up_xpm,Point<int>(0,0));
	//cursor[11].Init(arrow_right_xpm,Point<int>(31,16));
	cursor[12].Init(arrow_right_up_xpm,Point<int>(31,0));
	//cursor[13].Init(arrow_right_xpm,Point<int>(31,16));
	//  cursor[14].Init(arrow_right_xpm,Point<int>(31,16));
	// cursor[15].Init(arrow_right_xpm,Point<int>(31,16));
	cursor[16].Init(arrow_down_xpm,Point<int>(7,15));
	// cursor[17].Init(arrow_right_xpm,Point<int>(31,16));
	cursor[18].Init(arrow_left_down_xpm,Point<int>(0,31));

	cursor[20].Init(arrow_right_down_xpm,Point<int>(31,31));
}

void SDLCursor::CleanCursors() {
	
	cursor[0].Set();
	
	//std::cout << "CleanCursors" << std::endl;
	//cursor[0].Clean();
	cursor[1].Clean();
	cursor[2].Clean();
	cursor[4].Clean();
	cursor[8].Clean();
	cursor[10].Clean();
	cursor[12].Clean();
	cursor[16].Clean();
	cursor[18].Clean();
	cursor[20].Clean();
	
}

}
