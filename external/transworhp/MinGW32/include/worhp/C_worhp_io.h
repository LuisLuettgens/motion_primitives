#ifndef C_WORHP_IO_H
#define C_WORHP_IO_H

#include "macros.h"

void init_worhp_io(void);
int worhp_kbhit(void);
void set_keyboard(void);
void reset_keyboard(void);
void WorhpPrintWithLevel(const char message[], const int prn);


typedef void (*worhpPrintFunc) (int mode, const char s[]); 

DLL_PUBLIC void WorhpPrint(const int mode, const char message[]);
DLL_PUBLIC void SetWorhpPrint(worhpPrintFunc f);

/**
 * C printf replacement for WORHP.
 * Redirects to WorhpPrint, therefore adds the Fortran-like initial blank
 * and appends a newline to every string (when using the builtin WorhpPrint).
 *
 * Restriction: Resulting string length is limited to STRING_LENGTH.
 * Exceeding it is safe, but will result in truncation.
 */
DLL_PRIVATE void worhprintf(const char * format, ...);


/**
 * Prints an informative message including its origin, or continues a message.
 * @note Fortran-independent C implementation, since string interoperability
 * is inconvenient.
 * @see std::WorhpMessage
 */
DLL_PUBLIC void WorhpMessage(const char *message, const char *source, int prn);


/**
 * Prints a warning message including its origin, or continues a warning
 * message.
 * @note Fortran-independent C implementation, since string interoperability
 * is inconvenient.
 * @see std::WorhpError
 */
DLL_PUBLIC void WorhpWarning(const char *message, const char *source, int prn);


/**
 * Prints an error message including its origin, or continues an error
 * message.
 * @note Fortran-independent C implementation, since string interoperability
 * is inconvenient.
 * @see std::WorhpError
 */
DLL_PUBLIC void WorhpError(const char *message, const char *source, int prn);

enum {
  WORHP_PRINT_MESSAGE = 1,
  WORHP_PRINT_WARNING = 2,
  WORHP_PRINT_ERROR = 4
};

#endif
