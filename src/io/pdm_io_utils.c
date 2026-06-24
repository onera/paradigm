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


int
PDM_io_utils_diff_files
(
  const char *filename1,
  const char *filename2
)
{
  // https://www.wscubetech.com/resources/c-programming/programs/compare-two-files
  FILE *f1 = fopen(filename1, "r");
  if (f1 == NULL) {
    printf("PDM_io_utils_diff_files: Could not open %s\n", filename1);
    return -1;
  }

  FILE *f2 = fopen(filename2, "r");
  if (f2 == NULL) {
    printf("PDM_io_utils_diff_files: Could not open %s\n", filename2);
    return -2;
  }

  char line1[999], line2[999];

  int diff   = 0;
  int i_line = 1;

  // Compare each line from both files
  while (fgets(line1, sizeof(line1), f1) != NULL &&
         fgets(line2, sizeof(line2), f2) != NULL) {

    if (strcmp(line1, line2) != 0) {
      // We found different lines
      printf("Difference at line %d:\n", i_line);
      printf("  %s : %s\n", filename1, line1);
      printf("  %s : %s\n", filename2, line2);
      diff = i_line;
      break;
    }

    i_line++;
  }

  fclose(f1);
  fclose(f2);

  return diff;
}

