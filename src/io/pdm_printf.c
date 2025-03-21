/*============================================================================
 * Base user-definable PDM_printf() wrapper or replacement.
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

/*
 * Standard C library headers
 */
#include <stdarg.h>
#include <stdio.h>

/*
 * Optional library and BFT headers
 */

#include "pdm_printf.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/* Associated typedef documentation (for PDM_printf.h) */

/*!
 * \typedef PDM_printf_proxy_t
 *
 * \brief Function pointer for PDM_printf() type functions.
 *
 * \param [in] format       format string, as PDM_printf() and family.
 * \param [in, out] arg_ptr pointer to variable argument list based on format
 *                          string.
 */

/*!
 * \typedef PDM_printf_flush_proxy_t
 *
 * \brief Function pointer for fflush(stdout) type functions.
 */

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default PDM_printf_flush() proxy.
 *
 * returns:
 *   return code of fflush(stdout).
 */

static int
_PDM_printf_flush_proxy_default(void);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static PDM_printf_proxy_t        *_PDM_printf_proxy = vprintf;
static PDM_printf_flush_proxy_t  *_PDM_printf_flush_proxy
                                    = _PDM_printf_flush_proxy_default;

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Default PDM_printf_flush() proxy.
 *
 * returns:
 *   return code of fflush(stdout).
 */

static int
_PDM_printf_flush_proxy_default(void)
{
  return fflush(stdout);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

int
PDM_printf(const char *const format,
           ...)
{
  int  retval;
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  retval = _PDM_printf_proxy(format, arg_ptr);

  va_end(arg_ptr);

  return retval;
}

int
PDM_printf_flush(void)
{
  return _PDM_printf_flush_proxy();
}


PDM_printf_proxy_t *
PDM_printf_proxy_get(void)
{
  return _PDM_printf_proxy;
}

void
PDM_printf_proxy_set(PDM_printf_proxy_t *const fct)
{
  _PDM_printf_proxy = fct;
}

PDM_printf_flush_proxy_t *
PDM_printf_flush_proxy_get(void)
{
  return _PDM_printf_flush_proxy;
}

void
PDM_printf_flush_proxy_set(PDM_printf_flush_proxy_t *const fct)
{
  _PDM_printf_flush_proxy = fct;
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
