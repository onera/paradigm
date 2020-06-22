/*============================================================================
 * CUDA error handling
 *============================================================================*/

/*
  This file is part of the CWIPI library.

  Copyright (C) 2017 ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------------*/

#include <assert.h>

#include "pdm_cuda_error.cuh"
#include "pdm_cuda.cuh"
#include "pdm_error.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default error handler.
 *
 * An error message is output to stderr (stdout is flushed first),
 * and the current process is terminated.
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as PDM_printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

__device__
static void
_PDM_error_handler_default(const char     *const file_name,
                           const int             line_num,
                           const int             sys_error_code,
                           const char     *const format,
                           va_list               arg_ptr);

/*-----------------------------------------------------------------------------
* Local static variable definitions
*-----------------------------------------------------------------------------*/

__device__ static PDM_error_handler_t  *_PDM_error_handler = (_PDM_error_handler_default);

/*-----------------------------------------------------------------------------
* Local function definitions
*-----------------------------------------------------------------------------*/

/*
* Default error handler.
*
* An error message is output to stderr (after pdm_print_flush() is called),
* and the current process exits with an EXIT_FAILURE code.
*
* parameters:
*   file_name:      <-- name of source file from which error handler called.
*   line_num:       <-- line of source file from which error handler called.
*   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
*   format:         <-- format string, as PDM_printf() and family.
*   arg_ptr:        <-> variable argument list based on format string.
*/

__device__
static void
_PDM_error_handler_default(const char     *const file_name,
                           const int             line_num,
                           const int             sys_error_code,
                           const char     *const format,
                           va_list               arg_ptr)
{
  //PDM_printf_flush_GPU();

  //PDM_printf_GPU("\n");
  printf("\n");

  if (sys_error_code != 0)
    printf("\nSystem error\n");

  printf("\n%s:%d: Fatal error.\n\n", file_name, line_num);

  //PDM_printf_GPU(format, arg_ptr);

  printf("\n\n");

  assert(0);

  //exit(EXIT_FAILURE);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Calls the error handler (set by PDM_error_handler_set() or default).
 *
 * With the default error handler, PDM_print_flush() is called, an error
 * message is output to stderr, and the current process exits with an
 * EXIT_FAILURE code.
 *
 * \param [in] file_name      name of source file from which error handler
 *                            called.
 * \param [in] line_num       line of source file from which error handler
 *                            called.
 * \param [in] sys_error_code error code if error in system or libc call,
 *                            0 otherwise.
 * \param [in] format         format string, as PDM_printf() and family.
 * \param [in] ...            variable arguments based on format string.
 */


__device__
void
PDM_error_GPU(const char  *const file_name,
              const int          line_num,
              const int          sys_error_code,
              const char  *const format,
              ...)
{
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  _PDM_error_handler(file_name, line_num, sys_error_code, format, arg_ptr);

  va_end(arg_ptr);
}

/*!
* \brief Returns the error handler associated with the PDM_error() function.
*
* \return pointer to the error handler function.
*/

inline
__device__
PDM_error_handler_t *
PDM_error_handler_get_GPU(void)
{
  return _PDM_error_handler;
}

/*!
* \brief Associates an error handler with the PDM_error() function.
*
* \param [in] handler pointer to the error handler function.
*/

inline
__device__
void
PDM_error_handler_set_GPU(PDM_error_handler_t  *const handler)
{
  _PDM_error_handler = handler;
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
