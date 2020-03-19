#pragma once

#include <string>
#include <iostream>

#ifndef DllExport
#ifdef _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif
#endif

#ifdef _MSC_VER
#pragma warning(disable : 4251)
#include <windows.h>
#endif

/** Zaehlen und Timen von Aufrufen */
class DllExport Timing {
public:
	Timing();
	void Start();
	void Stop();
	void Reset();
	static void GetTime(std::string &str);
	double beforeUser;
	int count;

#ifdef _MSC_VER
	timeval beforeReal;
#else
	double beforeReal;
#endif
	double sumUser;
	double sumReal;

};

#define TIMING_CT 9
enum WORHPTiming_e {
	TIME_WORHP = 0,
	TIME_OUTPUT,
	TIME_F,
	TIME_G,
	TIME_DF,
	TIME_DG,
	TIME_HM,
	TIME_FIDIF,
	TIME_ALL
};

class DllExport WORHP_Timing {
public:
	WORHP_Timing();

	Timing timing[TIMING_CT];

	void Start(WORHPTiming_e t) {timing[t].Start();}
	void Stop(WORHPTiming_e t) {timing[t].Stop();}

	void PrintTimes(std::ostream &os);
	void Bar(std::ostream &os, double z);
};
