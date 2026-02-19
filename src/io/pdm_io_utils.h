/*
 * \file
 */

#ifndef __PDM_IO_UTILS_H__
#define __PDM_IO_UTILS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stddef.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief  Get file extension from a path
 *
 * \note The returned extension is not copied
 *
 * \param [in] filename Path to file
 *
 * \return File extension (not including the dot)
 */
const char *
PDM_io_utils_file_extension
(
  const char *filename
);


/**
 * \brief Get file name from a path
 *
 * \note The returned name is not copied
 *
 * \param [in] path Path to file
 *
 * \return File name (including the extension)
 */
const char *
PDM_io_utils_file_name_from_path
(
  const char *path
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_IO_UTILS_H__ */
