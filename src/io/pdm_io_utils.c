/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_io_utils.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

const char *
PDM_io_utils_file_extension
(
  const char *filename
)
{
  // https://stackoverflow.com/questions/5309471/getting-file-extension-in-c
  const char *dot = strrchr(filename, '.');
  if (!dot || dot == filename) {
    return "";
  }
  else {
    return dot + 1;
  }
}


const char *
PDM_io_utils_file_name_from_path
(
  const char *path
)
{
  // https://stackoverflow.com/questions/3288006/are-there-any-c-apis-to-extract-the-base-file-name-from-its-full-path-in-linux
  char delimiter = '/'; // use `\\' on Windows
  char *s = strrchr(path, delimiter);
  if (s == NULL) {
    return path;
  }
  else {
    return s+1;
  }
}

#ifdef  __cplusplus
}
#endif
