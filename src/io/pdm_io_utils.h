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


/**
 * \brief Compare two files line by line
 *
 * \param [in] filename1 Path to first file
 * \param [in] filename2 Path to second file
 *
 * \return (1-based) number of first different line if the files differ,
 *          0 if the files are identical
 *         -1 if the first file could not be opened
 *         -2 if the first file could not be opened
 */
int
PDM_io_utils_diff_files
(
  const char *filename1,
  const char *filename2
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_IO_UTILS_H__ */
