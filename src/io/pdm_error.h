#ifndef __PDM_ERROR_H__
#define __PDM_ERROR_H__

/*-----------------------------------------------------------------------------*/

/* Standard C library headers */

#include <stdarg.h>

/* BFT library headers */

#include "pdm_config.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

typedef void (PDM_error_handler_t) (const char    *const file_name,
                                    const int            line_num,
                                    const int            sys_error_code,
                                    const char    *const format,
                                    va_list              arg_ptr);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Calls the error handler (set by PDM_error_handler_set() or default).
 *
 * With the default error handler, PDM_print_flush() is called, an error
 * message is output to stderr, and the current process exits with an
 * EXIT_FAILURE code.
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   ... :           <-- variable arguments based on format string.
 */

void
PDM_error(const char  *const file_name,
          const int          line_num,
          const int          sys_error_code,
          const char  *const format,
          ...);

/*
 * Returns the error handler associated with the PDM_error() function.
 *
 * returns:
 *  pointer to the error handler function.
 */

PDM_error_handler_t *
PDM_error_handler_get(void);

/*
 * Associates an error handler with the PDM_error() function.
 *
 * parameters:
 *   handler: <-- pointer to the error handler function.
 */

void
PDM_error_handler_set(PDM_error_handler_t  *const handler);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_ERROR_H__ */
