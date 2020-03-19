#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#ifndef DllExport
	#ifdef _WIN32
		#define DllExport __declspec( dllexport )
	#else
		#define DllExport
	#endif
#endif

/** @defgroup textout Text Output
 *   Formatted text output.
 * @{
 */

enum TextOutputType_e {ASCII, RTF, HTML, CONSOLE, CPP, SILENT};

extern TextOutputType_e standard_textoutputtype;

/** @name Initializers
 * @{
 */
/** Write file header. */
std::ostream DllExport & beginfile( std::ostream& os );
/** Write file footer. */
std::ostream DllExport & endfile( std::ostream& os );
// @}

/** @name Page Manipulators
 * @{
 */
/** New line. */
std::ostream DllExport & endline ( std::ostream& os );
/** Tabulator */
std::ostream DllExport &tab ( std::ostream& os );
/** Add indentation by one level. */
std::ostream DllExport &addtab ( std::ostream& os );
/** Remove indentation by one level. */
std::ostream DllExport &removetab ( std::ostream& os );
// @}

/** @name Text Style Manipulators
 * @{
 */
/** Change text color. */
std::ostream DllExport &textcolor( std::ostream& os, int l );
/** Change to bold text. */
std::ostream DllExport &bold( std::ostream& os);
/** Change to bold text. */
std::ostream DllExport &oblique( std::ostream& os);
/** Change to normal text. */
std::ostream DllExport &normal( std::ostream& os);
/** Change text size. */
std::ostream DllExport &textsize( std::ostream& os, int l);
/** Change text size. */
std::ostream DllExport &proportional( std::ostream& os);
/** Change text size. */
std::ostream DllExport &monospace( std::ostream& os);
// @}

/** @name Special Characters
 * @{
 */
/** Wrute ". */
std::ostream DllExport &quote( std::ostream& os );
/** Write <. */
std::ostream DllExport &lt( std::ostream& os);
/** Write >. */
std::ostream DllExport &gt( std::ostream& os);
// @}

/** @brief A formatted ostream class
 *
 * A special std::ostream to write out simple formatted text (color, bold, etc.)
 */
class DebugStream : public std::ostream {

public:
	/** Constructor for file output. */
	DebugStream(const std::string &filename, TextOutputType_e totype) :
		std::ostream(nullptr),
		boldfont_m(false),
		obliquefont_m(false),
		colorfont_m(false),
		ident_m(0),
		filestream(filename.c_str()),
		TextOutputType_m(totype)
	{
		rdbuf(filestream.rdbuf());
		ident_m=0;
	}

	/** Constructor for console output. */
	DebugStream(TextOutputType_e totype) :
		std::ostream(nullptr),
		boldfont_m(false),
		obliquefont_m(false),
		colorfont_m(false),
		ident_m(0),
		datastream(),
		TextOutputType_m(totype)
	{
		rdbuf(datastream.rdbuf());
		ident_m=0;
	}

	/** Destructor */
	virtual ~DebugStream() {}

	/** Get std::string associated with this stream. */
	std::string GetString() const {
		return std::string(datastream.str());
	}

	/** Get text output type. */
	TextOutputType_e GetTextOutputType() {
		return TextOutputType_m;
	}

	bool boldfont_m;
	bool obliquefont_m;
	bool colorfont_m;
	int ident_m;

protected:
	std::ostringstream datastream;
	std::ofstream filestream;
	TextOutputType_e TextOutputType_m;

};

// @}

/* Textcolor Implementation (1 Argument manipulator) */
class textcolor_int {
public:
	textcolor_int(std::ostream& (*f)(std::ostream&,int), int t) : _fp(f), _tp(t) {}

	std::ostream& (* _fp)(std::ostream&,int);
	int _tp;
};

textcolor_int textcolor( int l );
std::ostream& operator<<(std::ostream& s, const textcolor_int & sm);



/* Textsize Implementation (1 Argument manipulator) */
class textsize_int {
public:
	textsize_int(std::ostream& (*f)(std::ostream&,int), int t) : _fp(f), _tp(t) {}

	std::ostream& (* _fp)(std::ostream&,int);
	int _tp;
};

textsize_int textsize( int l );
std::ostream DllExport & operator<<(std::ostream& s, const textsize_int & sm);




// for convenience

//#if 0
#define TEXTOUT(A) { if (debstr) (*debstr) << A ;}
#define TEXTOUT2(A) { if (collstr) (*collstr) << A ;}
//else
//#define TEXTOUT(A)
//#endif
