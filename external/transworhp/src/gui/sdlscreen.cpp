#include "sdlscreen.h"

#include "sdlframe.h"
#include "sdlthread.h"
#include "icon-transworhp.h"

#include "conversion.h"

#include "../base/point.h"

#include "../core/twstatus.h"

#include <string>
#include <cmath>

using namespace std;

namespace tw {

TWwindow::TWwindow() : size(800,600), resolution(800,600),
	reference(0,0),
	fullscreen(false),
	maximized(false),
	multisamplebuffers(0), multisamplesamples(0),
	windowmode(0) {

}


void TWwindow::ParseXML(XMLNode *n, int *countParams, int *setParams) {

	if (n) {
		(*countParams) += 6;
		string a( n->GetAttribute("width") );
		if (!a.empty()) {
			size.x = stoi(a);
			(*setParams)++;
		}
		a = n->GetAttribute("height");
		if (!a.empty()) {
			size.y = stoi(a);
			(*setParams)++;
		}

		a = n->GetAttribute("mode");
		if (a == "quiet") {
			windowmode = 1;
			(*setParams)++;
		} else if (a == "fullscreen") {
			fullscreen = true;
			(*setParams)++;
		} else if (a == "maximized") {
			maximized = true;
			(*setParams)++;
		} else if (a == "normal") {
			(*setParams)++;
		}

		a = n->GetAttribute("resolution");
		if (a == "default") {
			resolution = size;
			(*setParams)++;
		} else if (!a.empty()) {
			vector<int> res = ToIntArray(a,"x");
			resolution.x = res[0];
			resolution.y = res[1];
			(*setParams)++;
		}

		a = n->GetAttribute("x");
		if (!a.empty()) {
			reference.x = stoi(a);
			(*setParams)++;
		}

		a = n->GetAttribute("y");
		if (!a.empty()) {
			reference.y = stoi(a);
			(*setParams)++;
		}

		(*countParams) += 2;
		XMLNode *n2 = n->GetFirstChild("OPENGL");
		if (n2) {
			string a( n2->GetAttribute("multisamplebuffers") );
			if (!a.empty()) {
				multisamplebuffers = stoi(a);
				(*setParams)++;
			}
			a = n2->GetAttribute("multisamplesamples");
			if (!a.empty()) {
				multisamplesamples = stoi(a);
				(*setParams)++;
			}
		}
	}
}

SDL_Window* createWindow(const TWwindow& twwindow) {
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, twwindow.multisamplebuffers);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, twwindow.multisamplesamples);

	Uint32 flags = SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE;
	if (twwindow.fullscreen) {
		flags |= SDL_WINDOW_FULLSCREEN_DESKTOP;
	} else if (twwindow.maximized) {
		flags |= SDL_WINDOW_MAXIMIZED;
	}

	return SDL_CreateWindow("This will surely be overwritten", SDL_WINDOWPOS_CENTERED,
	                        SDL_WINDOWPOS_CENTERED, twwindow.resolution.x, twwindow.resolution.y,
	                        flags);
}

SDLScreen::SDLScreen(TWwindow *twwindow_)
		: image(nullptr), twwindow(*twwindow_) {


	SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_NOPARACHUTE);

	/*
	SDL_GL_STEREO             Enable or disable stereo (left and right) buffers (0 or 1).
	SDL_GL_MULTISAMPLEBUFFERS Number of multisample buffers (0 or 1). *
	SDL_GL_MULTISAMPLESAMPLES Number of samples per pixel when multisampling is enabled. *
	SDL_GL_ACCELERATED_VISUAL Guarantee hardware acceleration (0 or 1) **
	*/

	SDL_GL_SetAttribute(SDL_GL_RED_SIZE,            8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,          8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE,           8);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,          8);

	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,         16);
	SDL_GL_SetAttribute(SDL_GL_BUFFER_SIZE,        32);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,        1);

	SDL_GL_SetAttribute(SDL_GL_ACCUM_RED_SIZE,      8);
	SDL_GL_SetAttribute(SDL_GL_ACCUM_GREEN_SIZE,    8);
	SDL_GL_SetAttribute(SDL_GL_ACCUM_BLUE_SIZE,     8);
	SDL_GL_SetAttribute(SDL_GL_ACCUM_ALPHA_SIZE,    8);

	window = createWindow(twwindow);

	if (!window && (twwindow.multisamplebuffers || twwindow.multisamplesamples)) {
		twwindow.multisamplebuffers = twwindow.multisamplesamples = 0;
		window = createWindow(twwindow);
		if (window) {
			MyStatus("TransWORHP", "Multisample anti-aliasing not supported", Status::WARN);
		}
	}
	if (!window) {
		clog << "Couldn't set " << twwindow.resolution.x << "x" << twwindow.resolution.y
		     << " GL video mode: " << SDL_GetError() << endl
		     << "  - Invoke program with flag -p to optimize without graphic output" << endl;
		SDL_Quit();
		exit(2);
	}


	// Create an OpenGL context associated with the window.
	glcontext = SDL_GL_CreateContext(window);

	if (glcontext == nullptr){
		clog << "OpenGL context could not be created! SDL Error: " << SDL_GetError();
		SDL_Quit();
		exit(3);
	} else {

		// Init GLEW
		glewExperimental = GL_TRUE;
		GLenum glewError = glewInit();
		if(glewError != GLEW_OK) {
			clog << "Error initializing GLEW! " << glewGetErrorString(glewError);
			SDL_Quit();
			exit(3);
		} else {
			//use shaders
			//initGL();
		}
	}

	SDL_SetWindowTitle(window,"TransWORHP");

	SetIcon();

	SDL_version v;
	SDL_GetVersion(&v);
	const string version("SDL "+to_string(v.major)+"."+to_string(v.minor)+"."+to_string(v.patch));
	MyStatus("Version", version, Status::NORMAL);

	stringstream a;
	const GLubyte *c = glGetString(GL_VERSION);
	if (c) {
		a << "OpenGL " << c;
		MyStatus("Version", a.str(), Status::NORMAL);
	} else {
		MyStatus("Version", "No OpenGL Version found!", Status::NORMAL);
	}

	thethread = new SDLThread();
	thethread0 = thethread;
}


SDLScreen::~SDLScreen() {

	SDL_GL_DeleteContext(glcontext);

	delete thethread;
	thethread = nullptr;

	SDL_FreeSurface(image);

	//cout << "QUIT" << endl;
	SDL_Quit();
}

/*
 * soll in der Zukunft dazu genutzt werden um shader zu erstellen
 * vorhanden: 3 verschiedene "porgramme": normal, fuer textur, fuer spasePlot
 */
void SDLScreen::initGL() {

	//Generate program
	textureProgram = glCreateProgram();

	//Create vertex shader
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	//Get vertex source
	const GLchar* vertexShaderSource[] = {
//	  "#version 140\n"
//	  "in vec2 LVertexPos2D;"
//	  "uniform mat4 ProjectionModelviewMatrix;"
	  "void main() {"
//	  "gl_Position = vec4(LVertexPos2D.x, LVertexPos2D.y, 0, 1);"
	  "gl_FrontColor = gl_Color;"
	  "gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;"
	  "gl_TexCoord[0] = gl_MultiTexCoord0;"
	  "}"
	};

	//Set vertex source
	glShaderSource(vertexShader, 1, vertexShaderSource, NULL);
	//Compile vertex source
	glCompileShader(vertexShader);
	//Check vertex shader for errors
	GLint vShaderCompiled = GL_FALSE;
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &vShaderCompiled );
	if( vShaderCompiled != GL_TRUE) {
		clog << "Unable to compile vertex shader " << vertexShader;
		SDL_Quit();
		exit(3);
	}else {
		//Attach fragment shader to program
		glAttachShader(textureProgram, vertexShader);
	}

	//Create fragment shader
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	//Get fragment source
	const GLchar* fragmentShaderSource[] = {
	  "#version 120\n"
	  "uniform sampler2D tex;"
//	  "out vec4 LFragment;"
	  "void main() {"
//	  "LFragment = vec4( 1.0, 1.0, 1.0, 1.0 );"
	  "gl_FragColor = texture2D(tex, gl_TexCoord[0].st) * gl_Color;"
	  "}"
	};

	//Set fragment source
	glShaderSource(fragmentShader, 1, fragmentShaderSource, NULL );
	//Compile fragment source
	glCompileShader(fragmentShader);
	//Check fragment shader for errors
	GLint fShaderCompiled = GL_FALSE;
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &fShaderCompiled );
	if( fShaderCompiled != GL_TRUE ) {
		clog << "Unable to compile fragment shader " << fragmentShader;
		SDL_Quit();
		exit(3);
	} else {
		//Attach fragment shader to program
		glAttachShader(textureProgram, fragmentShader);
	}

	//Link program
	glLinkProgram(textureProgram);
	//Check for errors
	GLint programSuccess = GL_TRUE;
	glGetProgramiv(textureProgram, GL_LINK_STATUS, &programSuccess);
	if(programSuccess != GL_TRUE) {
		clog << "Error linking program " << textureProgram;
		SDL_Quit();
		exit(3);
	}



	///////////////////////////////////////////////////////////////////////////////////

	//Generate program
	simpleProgram = glCreateProgram();

	//Create vertex shader
	GLuint vertexShader2 = glCreateShader(GL_VERTEX_SHADER);
	//Get vertex source
	const GLchar* vertexShaderSource2[] = {
//	  "#version 140\n"
//	  "in vec2 LVertexPos2D;"
//	  "uniform mat4 ProjectionModelviewMatrix;"
	  "void main() {"
//	  "gl_Position = vec4(LVertexPos2D.x, LVertexPos2D.y, 0, 1);"
	  "gl_FrontColor = gl_Color;"
	  "gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;"
	  "}"
	};

	//Set vertex source
	glShaderSource(vertexShader2, 1, vertexShaderSource2, NULL);
	//Compile vertex source
	glCompileShader(vertexShader2);
	//Check vertex shader for errors
	GLint vShaderCompiled2 = GL_FALSE;
	glGetShaderiv(vertexShader2, GL_COMPILE_STATUS, &vShaderCompiled2 );
	if( vShaderCompiled2 != GL_TRUE) {
		clog << "Unable to compile vertex shader " << vertexShader2;
		SDL_Quit();
		exit(3);
	}else {
		//Attach fragment shader to program
		glAttachShader(simpleProgram, vertexShader2);
	}

	//Create fragment shader
	GLuint fragmentShader2 = glCreateShader(GL_FRAGMENT_SHADER);
	//Get fragment source
	const GLchar* fragmentShaderSource2[] = {
//	  "#version 140\n"
//	  "out vec4 LFragment;"
	  "void main() {"
//	  "LFragment = vec4( 1.0, 1.0, 1.0, 1.0 );"
	  "gl_FragColor = gl_Color;"
	  "}"
	};

	//Set fragment source
	glShaderSource(fragmentShader2, 1, fragmentShaderSource2, NULL );
	//Compile fragment source
	glCompileShader(fragmentShader2);
	//Check fragment shader for errors
	GLint fShaderCompiled2 = GL_FALSE;
	glGetShaderiv(fragmentShader2, GL_COMPILE_STATUS, &fShaderCompiled2 );
	if( fShaderCompiled2 != GL_TRUE ) {
		clog << "Unable to compile fragment shader " << fragmentShader2;
		SDL_Quit();
		exit(3);
	} else {
		//Attach fragment shader to program
		glAttachShader(simpleProgram, fragmentShader2);
	}

	//Link program
	glLinkProgram(simpleProgram);
	//Check for errors
	GLint programSuccess2 = GL_TRUE;
	glGetProgramiv(simpleProgram, GL_LINK_STATUS, &programSuccess2);
	if(programSuccess2 != GL_TRUE) {
		clog << "Error linking program " << simpleProgram;
		SDL_Quit();
		exit(3);
	}

	////////////////////////////////////////////////////////////////


	//Generate program
	sparsePlotProgram = glCreateProgram();


	//Create vertex shader
	GLuint vertexShader3 = glCreateShader(GL_VERTEX_SHADER);
	//Get vertex source
	const GLchar* vertexShaderSource3[] = {
//	  "#version 150\n"
	  "in vec2 in_Position;"
	  "out vec3 gl_Position;"
//	  "uniform mat4 ProjectionModelviewMatrix;"
	  "void main() {"
//	  "gl_Position = vec4(LVertexPos2D.x, LVertexPos2D.y, 0, 1);"
//	  "gl_FrontColor = gl_Color;"
	  "gl_Position = gl_ModelViewProjectionMatrix * vec3(in_Position,0.0);"
	  "}"
	};

	//Set vertex source
	glShaderSource(vertexShader3, 1, vertexShaderSource3, NULL);
	//Compile vertex source
	glCompileShader(vertexShader3);
	//Check vertex shader for errors
	GLint vShaderCompiled3 = GL_FALSE;
	glGetShaderiv(vertexShader3, GL_COMPILE_STATUS, &vShaderCompiled3 );
	if( vShaderCompiled3 != GL_TRUE) {
		clog << "Unable to compile vertex shader " << vertexShader3;

		GLint maxLength = 0;
		glGetProgramiv(sparsePlotProgram, GL_INFO_LOG_LENGTH, &maxLength);

		//The maxLength includes the NULL character
		std::vector<GLchar> infoLog(maxLength);
		glGetProgramInfoLog(sparsePlotProgram, maxLength, &maxLength, &infoLog[0]);

		//The program is useless now. So delete it.
		glDeleteProgram(sparsePlotProgram);

		for (auto ele : infoLog) {
			cout << ele;
		}
		cout << endl;

		SDL_Quit();
		exit(3);
	}else {
		//Attach fragment shader to program
		glAttachShader(sparsePlotProgram, vertexShader3);
	}

	//Create fragment shader
	GLuint fragmentShader3 = glCreateShader(GL_FRAGMENT_SHADER);
	//Get fragment source
	const GLchar* fragmentShaderSource3[] = {
	  "#version 150\n"
	  "flat in int in_Color;"
	  "out vec4 gl_FragColor;"
	  "void main() {"
	  "if(in_Color==3) {"
	  "gl_FragColor = vec4(.5,0.5,.5,.5);"
	  "} else {gl_FragColor = vec4(.7,0.7,.5,.5);}"
	  "}"
	};

	//Set fragment source
	glShaderSource(fragmentShader3, 1, fragmentShaderSource3, NULL );
	//Compile fragment source
	glCompileShader(fragmentShader3);
	//Check fragment shader for errors
	GLint fShaderCompiled3 = GL_FALSE;
	glGetShaderiv(fragmentShader3, GL_COMPILE_STATUS, &fShaderCompiled3 );
	if( fShaderCompiled3 != GL_TRUE ) {
		clog << "Unable to compile fragment shader " << fragmentShader3;
		SDL_Quit();
		exit(3);
	} else {
		//Attach fragment shader to program
		glAttachShader(sparsePlotProgram, fragmentShader3);
	}

	glBindAttribLocation(sparsePlotProgram, 0, "in_Position");
	glBindAttribLocation(sparsePlotProgram, 1, "in_Color");

	//Link program
	glLinkProgram(sparsePlotProgram);
	//Check for errors
	GLint programSuccess3 = GL_TRUE;
	glGetProgramiv(sparsePlotProgram, GL_LINK_STATUS, &programSuccess3);
	if(programSuccess3 != GL_TRUE) {
		clog << "Error linking program " << sparsePlotProgram;

		GLint maxLength = 0;
		glGetProgramiv(sparsePlotProgram, GL_INFO_LOG_LENGTH, &maxLength);

		//The maxLength includes the NULL character
		std::vector<GLchar> infoLog(maxLength);
		glGetProgramInfoLog(sparsePlotProgram, maxLength, &maxLength, &infoLog[0]);

		//The program is useless now. So delete it.
		glDeleteProgram(sparsePlotProgram);

		for (auto ele : infoLog) {
			cout << ele;
		}
		cout << endl;

		SDL_Quit();
		exit(3);
	}


	////////////////////////////////////////////////////////////////

	//Bind program
	glUseProgram(simpleProgram);
}




void SDLScreen::ToggleFullScreen(SDLFrame */*parent*/) {
// MATTHIAS

	/*
	if (twwindow.fullscreen) {
		screen = SDL_SetVideoMode(twwindow.resolution.x,twwindow.resolution.y,
								  16, SDL_OPENGL | SDL_RESIZABLE);
		if (screen) {
		    cout << "ToggleF" << screen->w << " "<< screen->h << endl;
			parent->Reshape(screen->w, screen->h);
		}
		SDL_ShowCursor(SDL_ENABLE);

		twwindow.fullscreen = false;
	}
	else {
		screen = SDL_SetVideoMode(twwindow.resolution.x,twwindow.resolution.y,
								  16, SDL_OPENGL | SDL_RESIZABLE);
		if (screen) {
		    cout << "ToggleF" << screen->w << " "<< screen->h << endl;
			parent->Reshape(screen->w, screen->h);
		}

		SDL_WM_ToggleFullScreen(screen);

		SDL_ShowCursor(SDL_DISABLE);
		twwindow.fullscreen = true;
	}*/
}



void SDLScreen::SetIcon() {

	image = SDL_CreateRGBSurfaceFrom(transworhplogo, 32, 32, 24, 32*3,
                                 0xff, 0x00ff00, 0xff0000, 0x0);

	SDL_SetWindowIcon(window, image);
}

int SDLScreen::Width() const {
	return twwindow.size.x;
}
int SDLScreen::Height() const {
	return twwindow.size.y;
}

}
