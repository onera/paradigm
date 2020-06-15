#ifndef __PDM_CUDA_ERROR_H__
#define __PDM_CUDA_ERROR_H__

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

/* Standard C library headers */

#include <stdarg.h>
#include <iostream>

/* BFT library headers */

#include "pdm_config.h"

/*-----------------------------------------------------------------------------*/

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
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

__device__
void
PDM_error_GPU(const char  *const file_name,
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

__device__
PDM_error_handler_t *
PDM_error_handler_get_GPU(void);

/*
 * Associates an error handler with the PDM_error() function.
 *
 * parameters:
 *   handler: <-- pointer to the error handler function.
 */

__device__
void
PDM_error_handler_set_GPU(PDM_error_handler_t  *const handler);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_CUDA_ERROR_H__ */
