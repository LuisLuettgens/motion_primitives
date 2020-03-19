#include "xopt_eps.h"

#include "xopt_data.h"

#include "conversion.h"

#include "../core/twstatus.h"

#include <iostream>
#include <fstream>
#include <sstream>

#ifdef WIN32
#include "windows.h"
#else
#include <pwd.h>
#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)<(B)?(A):(B))
#define fabs(A) ((A)>0?(A):-(A))
#endif

using namespace std;

namespace tw {

extern "C" double calltimefunc_(int &i,double *X, double *UNBE, double *T,
	                                int &len, int &ndgl, int &nunbe);


EpsWriter *globEpswriter=nullptr;

//extern int colormax;
int colormax = 10;
double stopmark;

EpsWriter::EpsWriter(const std::string &filename) : //_LBorder(40)
		Width2(312), Height2(217), _LBorder(40), _RBorder(20), _TBorder(20), _BBorder(20),
		absetzen(false),linewidth(1) {

	file.open(filename.c_str());

}

EpsWriter::~EpsWriter() {
	file.close();
}


void EpsWriter::Header(const std::string &Title) {

#ifndef WIN32

	time_t t;
	time(&t);

	//char buf[30]="";
	//sprintf(buf,"%02d.%02d.%04d", timeptr.tm_mday,timeptr.tm_mon, timeptr.tm_year);
	//"12.09.2003";
	//cout << buf << endl;

	//char hostname[30]="";
	//gethostname(hostname,sizeof(hostname));
	//passwd *p = getpwuid(geteuid());

	int sizes[] = {0,0,Width2,Height2};

	file << "%!PS-Adobe-3.0 EPSF-3.0\n";
	file << "%%Title: " << Title <<"\n";
	file << "%%Creator: XOPT Version 0.1 Patchlevel 2\n";
	file << "%%CreationDate: " << ctime(&t);//<<"\n";
	//file << "%%Pages: " << 0 << "\n";
	//file << "%%For: " << p->pw_name << "@" << hostname << " (" << p->pw_gecos<<")\n";
	file << "%%BoundingBox: "<< sizes[0] << " " << sizes[1] << " " << sizes[2]
	<< " " <<sizes[3] <<"\n";
	file << "%Magnification: 1.0000\n";
	//*os << "%%Comments: @ARGV\n";
	file << "%%EndComments\n";
#else

	//time_t t;
	//time(&t);

	//char buf[30]="";
	//sprintf(buf,"%02d.%02d.%04d", timeptr.tm_mday,timeptr.tm_mon, timeptr.tm_year);
	//"12.09.2003";
	//cout << buf << endl;

	//char hostname[30];
	//gethostname(hostname,sizeof(hostname));
	//passwd *p = getpwuid(geteuid());

	int sizes[] = {0,0,Width2,Height2};

	file << "%!PS-Adobe-3.0 EPSF-3.0\n";
	file << "%%Title: " << Title <<"\n";
	file << "%%Creator: XOPT Version 0.1 Patchlevel 2\n";
	//	file << "%%CreationDate: " << ctime(&t);//<<"\n";
	//file << "%%Pages: " << 0 << "\n";
	//	file << "%%For: " << p->pw_name << "@" << hostname << " (" << p->pw_gecos<<")\n";
	file << "%%BoundingBox: "<< sizes[0] << " " << sizes[1] << " " << sizes[2]
	<< " " <<sizes[3] <<"\n";
	file << "%Magnification: 1.0000\n";
	//*os << "%%Comments: @ARGV\n";
	file << "%%EndComments\n";
#endif



	//Preview
	/*file << "%%BeginPreview: 80 24 1 24\n";
	for (int i=0;i<8;i++) {
		file << "%FFFFFFFFFFFFFFFFFFFF\n";
	}
	for (int i=0;i<8;i++) {
		file << "%FF0000000000000000FF\n";
	}
	for (int i=0;i<8;i++) {
		file << "%FFFFFFFFFFFFFFFFFFFF\n";
	}
	file << "%%EndPreview\n";*/
	//file << "%%EndProlog\n";

	//file << "%%Page: 1 1\n";


	file << "/$F2psDict 200 dict def\n";
	file << "$F2psDict begin\n";
	file << "$F2psDict /mtrx matrix put\n";
	file << "/col-1 {0 setgray} bind def\n";

	//enum cols {Black, White, Red, Grey, Green, Cyan};
	file << "/col0 {0.000 0.000 0.000 srgb} bind def\n";
	file << "/col1 {0.800 0.800 0.800 srgb} bind def\n";
	//file << "/col2 {0.750 0.000 0.000 srgb} bind def\n";
	file << "/col2 {0.0235294 0.4 0.709804 srgb} bind def\n";
	file << "/col3 {0.700 0.700 0.700 srgb} bind def\n";
	file << "/col4 {0.000 0.750 0.000 srgb} bind def\n";
	//file << "/col5 {1.000 0.000 1.000 srgb} bind def\n";
	//file << "/col6 {1.000 1.000 0.000 srgb} bind def\n";
	file << "/col5 {.486 1 .062 srgb} bind def\n";
	file << "/col6 {1.000 .25 0.08 srgb} bind def\n";
	file << "/col7 {1.000 1.000 1.000 srgb} bind def\n";
	file << "/col8 {0.000 0.000 0.560 srgb} bind def\n";
	//file << "/col9 {0.000 0.000 0.890 srgb} bind def\n";
	file << "/col9 {.68 .49 .31 srgb} bind def\n";
	file << "/col21 {.3 .3 1. srgb} bind def\n";
	file << "/col22 {.7 .7 .7 srgb} bind def\n";
	file << "/col23 {.3 .8 .3 srgb} bind def\n";
	file << "/col24 {.8 .3 .3 srgb} bind def\n";

	// 173 125 80 = .68 .49 .31
	// 200 172 146
	//.93,.8,.68}  = 237 204 173
	//237 176 119 = .93 .69 .47
	//237 164 95 = .93 .64 .37
	//237 140 49 = .93 .55 .19
	// 69 130 181

	// 35 116 181

	for (int i=0;i<=colormax;i++) {
		double a = i/(1.*colormax);
		double c1[] ={6/255.,102/255.,181/255.};
		double c2[] ={(35+255)/510.,(116+255)/510.,(181+255)/510.};
		file << "/col" << (i+10) << " {"
		<< (1-a)*c1[0] + a*c2[0] << " "
		<< (1-a)*c1[1] + a*c2[1] << " "
		<< (1-a)*c1[2] + a*c2[2] << " "
		<< " srgb} bind def\n";
	}

	for (int i=0;i<=colormax;i++) {
		double a = i/(1.*colormax);
		double c1[] ={255/255.,32/255.,32/255.};
		double c2[] ={(255)/255.,(177)/255.,(43)/255.};
		file << "/col" << (i+110) << " {"
		<< (1-a)*c1[0] + a*c2[0] << " "
		<< (1-a)*c1[1] + a*c2[1] << " "
		<< (1-a)*c1[2] + a*c2[2] << " "
		<< " srgb} bind def\n";
	}

	/*file << "/col11 {0.000 0.000 .350 srgb} bind def\n";
	file << "/col12 {0.000 0.000 .400 srgb} bind def\n";
	file << "/col13 {0.000 0.000 .45000 srgb} bind def\n";
	file << "/col14 {0.000 0.000 .5000 srgb} bind def\n";
	file << "/col15 {0.000 0.000 .55000 srgb} bind def\n";
	file << "/col16 {0.000 0.000 .6000 srgb} bind def\n";
	file << "/col17 {0.000 0.000 .65000 srgb} bind def\n";
	file << "/col18 {0.000 0.000 .7000 srgb} bind def\n";
	file << "/col19 {0.000 0.000 .75000 srgb} bind def\n";
	file << "/col20 {0.000 0.000 .8000 srgb} bind def\n";
	file << "/col21 {0.000 0.000 .85000 srgb} bind def\n";
	file << "/col22 {0.000 0.000 .9000 srgb} bind def\n";
	file << "/col23 {0.000 0.000 .95000 srgb} bind def\n";
	file << "/col24 {0.000 0.000 1.000 srgb} bind def\n";*/
	/*	file << "/col26 {0.750 0.380 0.000 srgb} bind def\n";
		file << "/col27 {1.000 0.500 0.500 srgb} bind def\n";
		file << "/col28 {1.000 0.630 0.630 srgb} bind def\n";
		file << "/col29 {1.000 0.750 0.750 srgb} bind def\n";
		file << "/col30 {1.000 0.880 0.880 srgb} bind def\n";
		file << "/col31 {1.000 0.840 0.000 srgb} bind def\n";*/
	file << "end\n";
	file << "save\n";

	Rectangle((float)sizes[0],(float)sizes[3],(float)sizes[2],(float)sizes[1],true);

	file << "newpath\n";
	file << "0 " << Height2 <<" translate\n";

	file << 1. <<" "<< -1. <<" scale\n";

	file << "\n";
	/*	file << "/cp {closepath} bind def\n";
		file << "/ef {eofill} bind def\n";
		file << "/gr {grestore} bind def\n";
		file << "/gs {gsave} bind def\n";
		file << "/sa {save} bind def\n";
		file << "/rs {restore} bind def\n";
		file << "/l {lineto} bind def\n";
		file << "/m {moveto} bind def\n";
		file << "/rm {rmoveto} bind def\n";
		file << "/n {newpath} bind def\n";
		file << "/s {stroke} bind def\n";
		file << "/sh {show} bind def\n";*/
	/*file << "/rot {rotate} bind def\n";
	file << "/sd {setdash} bind def\n";

	file << "/sw {stringwidth} bind def\n";
	file << "/tr {translate} bind def\n";*/

	file << "%-------------- Procedures ---------------\n";
	file << "/srgb {setrgbcolor} bind def\n";
	file << "/textout {/yalign exch def /xalign exch def /text exch def \n";
	file << "  /Helvetica findfont 12.00 scalefont setfont\n";
	file << "  moveto\n"; // x y moveto
	file << "  gsave newpath 0 0 moveto text true charpath pathbbox grestore\n";
	file << "  /ury exch def /urx exch def\n";
	file << "  /lly exch def /llx exch def\n";

	file << "  xalign 0 eq { llx urx sub 2 div llx sub } if\n";
	file << "  xalign 1 eq { 0 } if\n";
	file << "  xalign 2 eq { -4 urx sub } if\n";

	file << "  yalign 0 eq { 4 } if\n"; //"lly ury sub 60 div lly sub\n";
	file << "  yalign 1 eq { 20 ury sub } if\n";
	file << "  yalign 2 eq { 0 } if\n";

	file << "  rmoveto gsave 1 -1 scale text col0 show grestore} bind def\n";


	file << "/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def\n";
	file << "/$F2psEnd {$F2psEnteredState restore end} def\n";
	file << "\n";
	file << "$F2psBegin\n";



}



void EpsWriter::LineWidth(float width) {
	linewidth = width;
	file << width << " setlinewidth\n";
}

//enum Colors {black=0, white,red,green,grey,darkgrey};

void EpsWriter::SetLineColor(TWcolor color)  {

	std::stringstream s;
	s << "col" << static_cast<int>(color);
	lineColor = s.str();
}

void EpsWriter::SetLineColorI(int  color)  {

	/*
		if (mode==0) {
			SetColor(Black);
		} else if (mode==1) {
			SetColor(White);
		} else if (mode==2) {
			SetColor(Green);
		} else if (mode==4) {
			SetColor(Cyan);
			glLineWidth(2.f);
		} else if (mode==3) {
			SetColor(White);
			glLineWidth(2.f);
		} else if (mode>=110) {
			SetColor(Hue2);
			glLineWidth(2.f);
		} else if (mode>=10) {
			SetColor(Hue);
			glLineWidth(2.f);
		}
	*/
	std::stringstream s;
	if (color==0)
		s << "col0";
	else if (color==1)
		s << "col9";
	else if (color==2)
		s << "col4";
	else if (color==3)
		s << "col1";
	else if (color==4)
		s << "col4";
	else if (color==5)
		s << "col5";
	else
		s << "col" << color;

	lineColor = s.str();
}


void EpsWriter::SetTextColor(TWcolor color)  {

	std::stringstream s;
	s << "col" << static_cast<int>(color);
	textColor = s.str();
}


void EpsWriter::Line(float x1, float y1, float x2, float y2) {

	file << "% Polyline\n";

	file << "newpath " << x1 << " " << y1 << " moveto "
	<< x2 << " " << y2 << " lineto gsave " << lineColor <<" stroke grestore\n";
}


void EpsWriter::Line() {

	file << "% Polyline <" << xx.size() <<">\n";

	auto x = xx.begin();
	auto y = yy.begin();

	file << "newpath " << *x << " " << *y << " moveto\n";

	x++;
	y++;

	for (;x!=xx.end();x++,y++)
		file << "  " << *x << " " << *y << " lineto\n";
	file << "  gsave " << lineColor << " stroke grestore\n";
}

void EpsWriter::XLine(const std::vector<float> &x, float y1, float y2) {

	file << "% XLine <" << x.size() <<">\n";

	auto it = x.begin();
	file << "[";
	for (;it!=x.end();it++) {
		if (it!=x.begin())
			file << " ";
		file << *it;
	}

	file << "]\n  {";

	file << "newpath dup " << y1 << " moveto " << y2
	<< " lineto gsave " << lineColor << " stroke grestore} forall\n";
}

void EpsWriter::YLine(const std::vector<float> &y, float x1, float x2) {

	file << "% YLine <" << y.size() <<">\n";

	auto it = y.begin();
	file << "[";
	for (;it!=y.end();it++) {
		if (it!=y.begin())
			file << " ";
		file << *it;
	}

	file << "]\n  {";

	file << "newpath dup " << x1 << " exch moveto " << x2
	<< " exch lineto gsave " << lineColor << " stroke grestore} forall\n";
}

void EpsWriter::Dot(float x1, float y1, float sz) {

	file << "% Dot\n";

	float lastw = linewidth;
	LineWidth(sz);

	file << "newpath " << x1 << " " << y1 << " moveto "
	<< x1 << " " << y1 << " lineto gsave " << lineColor <<" stroke grestore\n";

	LineWidth(lastw);
}


void EpsWriter::Dot(float sz) {

	file << "% Dot <" << xx.size() <<">\n";

	float lastw = linewidth;
	LineWidth(sz);

	auto x = xx.begin();
	auto y = yy.begin();

	for (;x!=xx.end();x++,y++)
		file << "newpath " << *x << " " << *y << " moveto " << *x << " " << *y
		<< " lineto gsave " << lineColor << " stroke grestore\n";

	LineWidth(lastw);
}


void EpsWriter::Rectangle(float x1, float y1, float x2, float y2, bool clip) {

	if (!clip)
		file << "% Rectangle\n";

	file << "newpath "
	<< x1 << " " << y1 << " moveto "
	<< x1 << " " << y2 << " lineto "
	<< x2 << " " << y2 << " lineto "
	<< x2 << " " << y1 << " lineto closepath ";

	if (clip)
		file << "clip\n";
	else
		file << "gsave " << lineColor << " stroke grestore\n";
}





void EpsWriter::XTextAlign(const std::map<float, std::string> &x, float y,
                           int xalign, int yalign) {

	file << "% XText <" << x.size() <<">\n";

	auto it = x.begin();
	file << "[";
	for (;it!=x.end();it++) {
		if (it!=x.begin())
			file << " ";
		file << "[" <<it->first << " (" << it->second << ")]";
	}

	file << "]\n  {";

	file << "aload pop " << int(y) << " exch "
	<< xalign << " " << yalign << " textout} forall\n";
}


void EpsWriter::YTextAlign(const std::map<float, std::string> &y, float x,
                           int xalign, int yalign) {

	file << "% YText <" << y.size() <<">\n";

	auto it = y.begin();
	file << "[";
	for (;it!=y.end();it++) {
		if (it!=y.begin())
			file << " ";
		file << "[" <<it->first << " (" << it->second << ")]";
	}

	file << "]\n  {";

	file << "aload pop " << int(x) << " 3 1 roll "
	<< xalign << " " << yalign << " textout} forall\n";
}





void EpsWriter::TextAlign(float x, float y, const std::string &text, int xalign, int yalign)  {

	file << "% Text\n";

	string text2 = text;


	if (text2[0]==char(132)) {
		text2.replace(0,1,"lambda");
	}

	file << int(x) << " " << int(y) << " (" << text2 << ") "
	<< xalign << " " << yalign << " textout\n";

	/*	file << "/Helvetica ff 12.00 scf sf " << int(x) << " " << int(y) << " m\n";

		file << "  gs newpath 0 0 m ("<< text << ") true charpath pathbbox gr\n";
		file << "  /ury exch def /urx exch def\n";
		file << "  /lly exch def /llx exch def\n";
		if (align & 1) {
			file << "  0\n";
		} else if (align & 2) {
			file << "  -4 urx sub\n";
		} else {
			file << "  llx urx sub 2 div llx sub\n";
		}

		if (align & 4) {
			file << "  20 ury sub\n";
		} else if (align & 8) {
			file << "  0\n";
		} else {
			file << "  4\n";//"lly ury sub 60 div lly sub\n";
		}
		file << "  rm gs 1 -1 sc (" << text << ") " << textColor <<" sh gr\n";
		*/
}



void EpsWriter::Footer() {
	file << "$F2psEnd\n";
	file << "restore\n";
	file << "showpage\n";
	file << "%%Trailer\n";
	file << "%EOF\n";
}

void BasePlot::Matlab() {

	static int i=0;
	i++;

	string filename = "pix_" + ToString(i) + "_" + epstitle + ".m";
	MyStatus("MATLAB", "Writing " + filename, Status::NORMAL);

	ofstream of(filename.c_str());
	matlab(of);
}


void BasePlot::Print(DataStorage &ds, DataStorage &dstop) {

	static int i=0;
	i++;

	string filename = "pix_" + ToString(i) + "_" + epstitle + ".eps";
	MyStatus("EPS", "Writing " + filename, Status::NORMAL);
	EpsWriter epswriter(filename);

	epswriter.High = High;
	epswriter.Low = Low;
	epswriter.High2 = High2;
	epswriter.Low2 = Low2;
	epswriter.Mirror = Mirror;


	epswriter.Header(epstitle);
	epswriter.Frame(epstitle, this);

	epsCompareCurve(&epswriter);
	epsData(ds, dstop, &epswriter);

	epswriter.Footer();

	//delete epswriter;
}




void EpsWriter::Frame(const std::string &title, BasePlot *pl) {

	file << "% Begin Frame\n";
	file << "0 setlinejoin\n";
	file << "2 setlinecap\n";
	file << "10 setmiterlimit\n";
	LineWidth(1);

	SetLineColor(TWcolor::Grey);// = "col1";
	SetTextColor(TWcolor::Black);// = "col0";

	float xzero=epsMapToX(1,3);
	TextAlign(xzero, 0,title, 0,1);

	std::vector<float> x_lines;
	std::map<float, std::string> x_text;

	std::vector<float> y_lines;
	std::map<float, std::string> y_text;


	double step = 30*(High2-Low2)/Width2;
	step = pl->roundStep(step); // MATTHIAS


	for (double i=0;i<=High2+.001;i+=step) {
		xzero=epsMapFloatToX(i);
		if (xzero>_LBorder+2 && xzero<Width2-_RBorder-2) {
			if (i>0) {
				x_lines.push_back(xzero);
			}
		}

		if (xzero>=_LBorder-2 && xzero<=Width2-_RBorder+2) {
			x_text[xzero] = ticktext(step,i);
		}
	}

	for (double i=-step;i>=Low2-.001;i-=step) {
		xzero=epsMapFloatToX(i);

		if (xzero>_LBorder+2 && xzero<Width2-_RBorder-2) {
			x_lines.push_back(xzero);
		}

		if (xzero>=_LBorder-2 && xzero<=Width2-_RBorder+2) {
			x_text[xzero] = ticktext(step,i);
		}
	}

	XLine(x_lines, (float)_TBorder, (float)(Height2-_BBorder));
	XTextAlign(x_text,(float)(Height2-_BBorder), 0, 1);


	step = 20*(High-Low)/Height2;
	step = pl->roundStep(step);

	double l = max(0.,Low);
	l = (int)(l/step)*step;

	for (double i=l;i<=High+.001;i+=step) {
		float yzero=epsMapToY(i);

		if (yzero>_TBorder+2 && yzero<Height2-_BBorder-2) {
			if (i>0) {
				y_lines.push_back(yzero);
			}
		}

		if (yzero>=_TBorder-2 && yzero<=Height2-_BBorder+2) {
			y_text[yzero] = ticktext2(step,i);
		}
	}

	l = min(0.,High);
	l = (int)(l/step)*step;
	if (fabs(l)<.001)
		l= -step;

	for (double i=l;i>=Low-.001;i-=step) {
		float yzero=epsMapToY(i);

		if (yzero>_TBorder+2 && yzero<Height2-_BBorder-2) {
			y_lines.push_back(yzero);
		}

		if (yzero>=_TBorder-2 && yzero<=Height2-_BBorder+2) {
			y_text[yzero] = ticktext2(step,i);
		}
	}

	YLine(y_lines, (float)_LBorder,(float)(Width2-_RBorder));
	YTextAlign(y_text,(float)_LBorder, 2,0);


	LineWidth(1.2f);


	SetLineColor(TWcolor::Black);


	xzero=epsMapFloatToX(0);
	if (xzero>_LBorder+2 && xzero<Width2-_RBorder-2) {
		Line((float)xzero,(float)_TBorder,(float)xzero,(float)(Height2-_BBorder));
	}

	XLine(x_lines,(float)_TBorder,_TBorder+4.f);
	XLine(x_lines,(float)(Height2-_BBorder),Height2-_BBorder-4.f);


	float yzero=epsMapToY(0);
	if (yzero>_TBorder+2 && yzero<Height2-_BBorder-2) {
		Line(_LBorder,yzero,Width2-_RBorder,yzero);
	}

	YLine(y_lines,(float)_LBorder,(float)(_LBorder+4));
	YLine(y_lines,(float)(Width2-_RBorder),(float)(Width2-_RBorder-4));


	SetLineColor(TWcolor::Black);
	Rectangle((float)_LBorder,(float)_TBorder,(float)(Width2-_RBorder),(float)(Height2-_BBorder));

	Rectangle((float)_LBorder+1,(float)_TBorder+1,(float)(Width2-_RBorder-1),(float)(Height2-_BBorder-1),true);

	file << "% End Frame\n";

	file << "1 setlinejoin\n";
	file << "1 setlinecap\n";
	LineWidth(1.2);

}







void EpsWriter::AddPoint(float x, float y) {
	xx.push_back(x);
	yy.push_back(y);
}

void EpsWriter::ClearPoints() {
	xx.clear();
	yy.clear();
}


void BasePlot::epsCompareCurve(EpsWriter */*epsw*/) const {
/*
	if (comp) {

		epsw->SetLineColor(Grey); // "col4";

		if (compstep==0) {
			float cury = epsw->epsMapToY(comp[0]);
			epsw->Line(0,cury,epsw->_LBorder,cury);
			epsw->Line(epsw->Width2-epsw->_RBorder,cury,epsw->Width2,cury);
		} else if (comptime==0) {
			for (int i=0;i<gl->ds.GetLength();i++) {
				epsw->AddPoint(epsw->epsMapTimeToX(i), epsw->epsMapToY(comp[compstep*i]));
			}
			epsw->Line();
			epsw->ClearPoints();
		} else {

			for (int i=0;i<gl->ds.GetLength();i++) {
				epsw->AddPoint(epsw->epsMapFloatToX(comptime[i]), epsw->epsMapToY(comp[compstep*i]));
			}
			epsw->Line();
			epsw->ClearPoints();


			/ *	SetColor(Grey);
				//std::cout << "COMPAREC" << comptime[40] <<  std::endl;

				glBegin(GL_LINE_STRIP);
				for (int i=0;i<gl->ds.GetLength();i++) {

					//	std::cout << "COMPAREC" << comp[compstep*i]
					//	<< " " <<  comptime[i] <<  std::endl;
					cury = MapToY(comp[compstep*i]);
					curx = MapFloatToX(comptime[i]);

					glVertex3f(curx, cury,0.5);

				}
				glEnd();* /

		}
	}*/
}


float EpsWriter::epsMapToY(double d) const {
	double yptscale = 1./(High-Low)*((double)(Height2-_TBorder-_BBorder));
	if (Mirror)
		return (float)(((d-Low)*yptscale)+_TBorder);
	else
		return (float)(((High-d)*yptscale)+_TBorder);
}

float EpsWriter::epsMapToX(int x, int ndis) const {
	return (float)(((double)x/(double)(ndis-1)*
	        ((double)(Width2-_LBorder-_RBorder))-.5)+_LBorder);
}

float EpsWriter::epsMapFloatToX(double d) const {
	return (float)(((d-Low2)/(High2-Low2)*
	        ((double)(Width2-_LBorder-_RBorder))-.5)+_LBorder);
}



float EpsWriter::epsMapTimeToX(DataStorage &ds, int i) {

	float ww = Width2-_LBorder-_RBorder;
	float w2 = _LBorder;

	if (ds.timemode==DataStorage::time_as_func) {
		double s = i; // calltimefunc_(i,ds.X.D,ds.UNKNOWN.D,ds.T.D,
		                //         ds.X.ndis,ds.X.size,ds.UNKNOWN.size) / ds.GetTF();
		return (float)(s*ww+.5)+w2;
	}

	if (ds.timemode==DataStorage::time_is_time) {
		double s = ds.getData(Selector('t',0),i)/ds.getData(Selector('t',0),ds.getLength()-1);
		return (float)(s*ww+.5)+w2;
	}

	if (ds.timemode==DataStorage::index_as_time) {

		int a = ds.getTF();
		//cout << "i " << i << " " << a << " " << i-1
		return (float)((i)/(a-1.)*ww+.5)+w2;

		/*return (int)((gl->ds.GetData(0,i)-1)/(gl->ds.GetLength()-1)
		             *(ww)-.5)+w2;*/
	}

	if (ds.timemode==DataStorage::dgl2_as_time) {

		absetzen=false;
		double t = ds.getData(Selector('t',0),i)*ds.getData(Selector('x',DataStorage::timedgl),i);
		if (i>=DataStorage::nsplitdis) {

			int block = ((i-DataStorage::nsplitdis)/DataStorage::nsplitdis2);
			int b = DataStorage::nsplitindex[block]-1;
			double dt = ds.getData(Selector('t',0),b)*ds.getData(Selector('x',DataStorage::timedgl),b);

			if ((i-DataStorage::nsplitdis)%DataStorage::nsplitdis2==0)
				absetzen=true;
			else
				absetzen=false;
			//	std::cout << "XXX " << i << " "<< block << " " << b << " " << t << " " <<dt<< std::endl;
			t +=dt;

		}
		return (short)(t/(ds.getTF())*((double)(Width2-_LBorder-_RBorder))+.5)+_LBorder;
	}

	return (ds.getData(Selector('t',0),i)*((double)(Width2-_LBorder-_RBorder))-.5)+_LBorder;

}

float EpsWriter::epsMapStopTimeToX(DataStorage &ds,DataStorage &dstop, int i) const {

	float a = dstop.getData(Selector('t',0),i)*dstop.getTF()+stopmark;
	return (a/ds.getTF()*((double)(Width2-_LBorder-_RBorder))-.5)+_LBorder;

}


 /*
void GeneralPlot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	epsw->SetLineColor(Red);// "col2";

	for (int i=0;i<ds.GetLength();i++) {
		epsw->AddPoint(epsw->epsMapTimeToX(ds,i), epsw->epsMapToY(ds.GetData(Dgl,i)));
	}
	if (!epsw->absetzen)
		epsw->Line();
	else {
		epsw->SetLineColor(Green);
		epsw->Dot(ds.GetLength());
		epsw->LineWidth(.6);
	}
	epsw->ClearPoints();


	if (dstop.GetLength()) {
		epsw->SetLineColor(Green);// "col3";
		epsw->LineWidth(1.2);

		for (int i=0;i<dstop.GetLength();i++) {
			epsw->AddPoint(epsw->epsMapStopTimeToX(ds,dstop, i),epsw->epsMapToY(dstop.GetData(Dgl,i)));
		}
		epsw->Line();
		epsw->ClearPoints();

	}

}


void PhasePlot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	epsw->SetLineColor(Red);// "col2";

	for (int i=0;i<ds.GetLength();i++) {
		epsw->AddPoint(epsw->epsMapFloatToX(ds.GetData(Dgl2,i)),
		               epsw->epsMapToY(ds.GetData(Dgl,i)));
	}

	epsw->Line();
	epsw->ClearPoints();


	if (dstop.GetLength()) {

		epsw->SetLineColor(Green);// "col2";

		for (int i=0;i<dstop.GetLength();i++) {
			epsw->AddPoint(epsw->epsMapFloatToX(dstop.GetData(Dgl2,i)),
			               epsw->epsMapToY(dstop.GetData(Dgl,i)));
		}
		epsw->Line();
		epsw->ClearPoints();
	}
}





void UserPlot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	globEpswriter = epsw;

	epsw->SetLineColor(Red);//="col2";
	int f = 0; //=*dpy->GetFrame();

	/ *if (dpy->GetTime()>0)
		f = (int)(dpy->GetTime()*1000.);
* /

	if (Func)
		ds.Call(Func,index,f,userBorder);
	else
		ds.Call(FuncI,index,f,userBorder,iindex);

	if (dstop.GetLength()) {
		epsw->SetLineColor(Grey);//epscol="col3";
		if (Func)
			dstop.Call(Func,index,f,userBorder);
		else
			dstop.Call(FuncI,index,f,userBorder,iindex);
	}

	globEpswriter = 0;

}


void DataPlot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	std::vector<int>::iterator it = indices.begin();
	for (;it!=indices.end();it++) {

		epsw->SetLineColor(Red);// "col2";
		for (int i=0;i<ds.GetLength();i++) {
			epsw->AddPoint(epsw->epsMapTimeToX(ds, i),
			               epsw->epsMapToY(ds.Call(Func,i,*it)));
		}
		epsw->Line();
		epsw->ClearPoints();


		if (dstop.GetLength()) {
			epsw->SetLineColor(Green);// "col3";

			for (int i=0;i<dstop.GetLength();i++) {
				epsw->AddPoint(epsw->epsMapStopTimeToX(ds, dstop, i),
				               epsw->epsMapToY(dstop.Call(Func,i,*it)));
			}
			epsw->Line();
			epsw->ClearPoints();
		}
	}
}

*/

}
