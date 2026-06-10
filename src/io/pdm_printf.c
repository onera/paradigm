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
