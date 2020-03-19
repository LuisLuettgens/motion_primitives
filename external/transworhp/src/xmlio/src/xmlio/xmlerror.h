#ifndef xmlerror_h
#define xmlerror_h

#include <string>
#include <vector>

/** @defgroup error Error Codes
 *  Error Codes for XML parsing
 * @{
 */
/** ewr */
#define XML_ATTRIBUTE_NOT_QUOTED 10
/** ewr */
#define XML_NO_ELEMENT_NAME 11
/** ewr */
#define XML_ATTRIBUTE_WITHOUT_ARGUMENTS 12
/** ewr */
#define XML_EQUAL_EXPECTED 13
/** ewr */
#define XML_NO_ROOT 14
/** ewr */
#define XML_ELEMENT_NOT_CLOSED 15
/** ewr */
#define XML_WRONG_ROOT_NAME 16
/** ewr */
#define XML_UNEXPECTED_CHARACTER 17
/** ewr */
#define XML_FILE_NOT_FOUND 18
// @}


/** @ingroup xml
 *  @brief An exception class.
 *
 *  Error tracking in XML files.
 */
class XMLError {
public:
	/** Constructor
	 * @param c Code
	 * @param l Line number
	 * @param ca Character (unused)
	 * @param s Additional information for error tracking
	 */
	XMLError(int c, int l, int ca, const std::string &s) :
	code(c), line(l), character(ca), info(s) {}

	/** Destructor */
	~XMLError() {}

	/** Get message for internal value "code". */
	std::string GetErrorText() const;

	/** ostream this class: "XMLError: error text (line number)". */
	friend std::ostream &operator << (std::ostream &os, const XMLError &e);

protected:
	int code;
	int line;
	int character;
	std::string info;

};


#endif
