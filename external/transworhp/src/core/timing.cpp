#include "timing.h"

#include <ctime>
#include <sstream>
#include <iomanip>

#ifdef WIN32
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif


#ifdef _MSC_VER
struct timezone
{
	int tz_minuteswest; /* minutes W of Greenwich */
	int tz_dsttime;     /* type of dst correction */
};
#endif


int gettimeofday1(struct timeval *tv, struct timezone *tz)
{
	FILETIME ft;
	unsigned __int64 tmpres = 0;
	static int tzflag;

	if (NULL != tv)
	{
		GetSystemTimeAsFileTime(&ft);

		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;

		/*converting file time to unix epoch*/
		tmpres -= DELTA_EPOCH_IN_MICROSECS;
		tmpres /= 10;  /*convert into microseconds*/
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	}

	if (NULL != tz)
	{
		if (!tzflag)
		{
			_tzset();
			tzflag++;
		}
		tz->tz_minuteswest = _timezone / 60;
		tz->tz_dsttime = _daylight;
	}

	return 0;
}

#else
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

Timing::Timing() : count(0), sumUser(0.0), sumReal(0.0) {}

void Timing::Start() {

	count++;

	//beforeCPU =  (double) clock() / (double) CLOCKS_PER_SEC;

#ifdef WIN32

	FILETIME createTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;
	if ( GetProcessTimes( GetCurrentProcess( ), &createTime, &exitTime, &kernelTime, &userTime ) != -1 ) {

		SYSTEMTIME userSystemTime;
		if ( FileTimeToSystemTime( &userTime, &userSystemTime ) != -1 )
			beforeUser=(double)userSystemTime.wHour * 3600.0 +
			(double)userSystemTime.wMinute * 60.0 +
			(double)userSystemTime.wSecond +
			(double)userSystemTime.wMilliseconds / 1000.0;

		/*	SYSTEMTIME kernelSystemTime;
		if ( FileTimeToSystemTime( &kernelTime, &kernelSystemTime ) != -1 )
		beforeKernel=(double)kernelSystemTime.wHour * 3600.0 +
		(double)kernelSystemTime.wMinute * 60.0 +
		(double)kernelSystemTime.wSecond +
		(double)kernelSystemTime.wMilliseconds / 1000.0;*/

	}


	//#ifdef _OPENMP
	//beforeReal = omp_get_wtime();
	gettimeofday1(&beforeReal, NULL);

	//#endif
#else
	beforeUser =  static_cast<double>(clock()) / static_cast<double>(CLOCKS_PER_SEC);

#ifdef _OPENMP
	beforeReal = omp_get_wtime();
#else
	beforeReal = 0.0;
#endif
#endif

}

void Timing::Stop() {
	//sumCPU += (double) clock() / (double) CLOCKS_PER_SEC - beforeCPU;


#ifdef WIN32
	FILETIME createTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;
	if ( GetProcessTimes( GetCurrentProcess( ),
		&createTime, &exitTime, &kernelTime, &userTime ) != -1 )
	{
		SYSTEMTIME userSystemTime;
		if ( FileTimeToSystemTime( &userTime, &userSystemTime ) != -1 ) {
			double t=(double)userSystemTime.wHour * 3600.0 +
				(double)userSystemTime.wMinute * 60.0 +
				(double)userSystemTime.wSecond +
				(double)userSystemTime.wMilliseconds / 1000.0;
			sumUser += t-beforeUser;
		}
		/*		SYSTEMTIME kernelSystemTime;
		if ( FileTimeToSystemTime( &kernelTime, &kernelSystemTime ) != -1 ) {
		double t=(double)kernelSystemTime.wHour * 3600.0 +
		(double)kernelSystemTime.wMinute * 60.0 +
		(double)kernelSystemTime.wSecond +
		(double)kernelSystemTime.wMilliseconds / 1000.0;

		sumKernel += t-beforeKernel;
		}*/
	}
	//#ifdef _OPENMP
	timeval endReal;
	gettimeofday1(&endReal,0);
	sumReal += (endReal.tv_sec - beforeReal.tv_sec)  + ( endReal.tv_usec - beforeReal.tv_usec)*1e-6;
	//sumReal += //omp_get_wtime() - beforeReal;
	//#endif

#else

	sumUser += static_cast<double>(clock()) / static_cast<double>(CLOCKS_PER_SEC) - beforeUser;

#ifdef _OPENMP
	sumReal += omp_get_wtime() - beforeReal;
#else
	sumReal += 0.0;
#endif

#endif

}

void Timing::Reset() {
	sumUser = 0.0;
	sumReal = 0.0;
	count = 0;
}


void Timing::GetTime(std::string &str) {

#ifdef WIN32
	SYSTEMTIME SystemTime;
	GetLocalTime(&SystemTime);

	int h = SystemTime.wHour;
	int m = SystemTime.wMinute;
	int s = SystemTime.wSecond;
	int ms = SystemTime.wMilliseconds/10;

	char buf[100];
	sprintf(buf,"[%02d:%02d:%02d.%02d]  ", h, m, s, ms );

	str = std::string(buf);
#else
	str = "";
#endif
}


const std::string WORHP_Timing_Titles[] = {"WORHP", "Output", "Objective F", "Constraints G", "Gradient DF", "Jacobian DG", "Hessian HM", "Finite Diff", "TransWORHP"};

WORHP_Timing::WORHP_Timing() {

	for (auto & elem : timing) {
		elem.Reset();
	}
}

void WORHP_Timing::PrintTimes(std::ostream &os) {

	os << "                     " << std::setw(15) << " "
		<< std::setw(10) << "count"
		<< std::setw(12) << "Real Time"
		<< std::setw(12) << "CPU Time" << std::endl;
	os << std::setprecision(3);
	for (int i = 0; i < TIMING_CT; i++) {
		os << "                     ";
		os.setf(std::ios::left, std::ios::adjustfield);
		os << std::setw(15) << WORHP_Timing_Titles[i];
		os.setf(std::ios::fixed | std::ios::right, std::ios::adjustfield);
		os.setf(std::ios::fixed);
		os << std::setw(10) << timing[i].count;
		os << std::setw(12) << timing[i].sumReal;
		os << std::setw(12) << timing[i].sumUser;
		Bar(os,timing[i].sumUser / timing[TIME_ALL].sumUser);
		os << std::endl;
	}

}

void WORHP_Timing::Bar(std::ostream &os, double z) {
	os << " [";
	const int NMAX = 20;
	for (int i = 0; i < NMAX; i++) {
		if (i*100./NMAX < z*100) {
			char b = '#';
			os << b ;
		}
		else { os << "-";}
	}
	os << "]";
}
