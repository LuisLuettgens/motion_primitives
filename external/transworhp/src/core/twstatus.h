#pragma once

#include <string>

#ifndef DllExport
#ifdef _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif
#endif

#ifdef _MSC_VER
#pragma warning(disable : 4251)
#endif

enum class Status {
	BOLD      = 1,
	RED       = 2,
	GREEN     = 4,
	BLUE      = 8,
	CONT_NEXT = 16,
	CONT_PREV = 32,
	NORMAL    = Status::BLUE | Status::GREEN | Status::RED,
	INFO      = Status::BOLD | Status::BLUE,
	WARN      = Status::BOLD | Status::GREEN | Status::RED,
	ERR       = Status::BOLD | Status::RED,
};

Status operator|(Status a, Status b);

using tw_output = void (*)(const char* tag, const char* s, Status flag);

void DllExport MyStatus(const std::string &tag, const std::string &msg, Status flag);
//void DllExport MyOutputFunction(ext_io_func f);
void DllExport SetMyStatusStream(std::ostream &os);
void DllExport SetMyOutputFunction(tw_output f);
