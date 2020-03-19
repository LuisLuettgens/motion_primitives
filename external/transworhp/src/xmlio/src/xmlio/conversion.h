#pragma once

#include <string>
#include <vector>
#include <stdint.h>

#ifndef DllExport
	#ifdef _WIN32
		#define DllExport __declspec( dllexport )
	#else
		#define DllExport
	#endif
#endif

/** @defgroup convert Conversion routines
 *  Convert freely between numbers and std::strings.
 * @{
 */

/** @name Convert std::string to numbers.
* @{
*/
/** Convert string to integer. */
int DllExport ToInt(const std::string &s);
/** Convert string to double. */
double DllExport ToDouble(const std::string &s);
float DllExport ToFloat(const std::string &s);
/** Convert string to bool. */
bool DllExport ToBool(const std::string &s);

void DllExport InitLocale();
// @}

/** @name Convert numbers to std::string.
* @{
*/
/** Convert integer to string. */
std::string DllExport ToString(int a);
std::string DllExport ToString(size_t a);

#ifdef WIN32
std::string DllExport ToString(unsigned long a);
#endif

std::string DllExport ToHexString(unsigned long a);
/** Convert double to string. */
std::string DllExport ToString(double a);
std::string DllExport ToScientificString(double a);
std::string DllExport ToString(double a, int width, int prec);
std::string DllExport ToString(char a);

/** Add PointZero if no point is in number */
std::string DllExport ToStringFakePoint(double a);
// @}

/** @name Split and join std::strings.
* @{
*/
/** Split string by space or tab. */
std::vector<std::string> DllExport ToStringArray(const std::string &s);
/** Auch leere Strings */
std::vector<std::string> DllExport ToStringArrayEmpty(const std::string &s, const std::string &split);
/** Split string by any char from split. */
std::vector<std::string> DllExport ToStringArray(const std::string &s, const std::string &split);
/** Split string by any char from split returning the first part and the rest (max 2 entries) */
std::vector<std::string> DllExport ToStringArray2(const std::string &s, const std::string &split);
/** Split string by space, tab or comma and convert to double. */
std::vector<double> DllExport ToDoubleArray(const std::string &s);
/** Split string by any char from split. */
std::vector<double> DllExport ToDoubleArray(const std::string &s, const std::string &split);
std::vector<double> DllExport ToDoubleArrayEmpty(const std::string &s, const std::string &split);
/** Split string by space, tab or comma and convert to float. */
std::vector<float> DllExport ToFloatArray(const std::string &s);
/** Split string by any char from split. */
std::vector<float> DllExport ToFloatArray(const std::string &s, const std::string &split);
/** Split string by space, tab or comma and convert to int. */
std::vector<int> DllExport ToIntArray(const std::string &s);
/** Split string by any char from split. */
std::vector<int> DllExport ToIntArray(const std::string &s, const std::string &split);
/** Join a vector of strings, j is the glue. */
std::string DllExport Join(const std::vector<std::string> &s, const std::string &j=std::string(" "));
std::string DllExport Join(const std::vector<int> &s, const std::string &j = std::string(" "));
std::string DllExport Join(const std::vector<double> &s, const std::string &j = std::string(" "));
// @}
std::string DllExport Replace(std::string &s, const std::string &r1, const std::string &r2);

std::string DllExport eatwhitespace(std::string &c);
std::string DllExport replacebrackets(std::string &c);
std::string DllExport replacequote(std::string &c);

void DllExport myReplace(std::string& s, std::string toReplace, std::string replaceWith);

std::vector<double> ToMatlabArray(const std::string &achse);

// @}