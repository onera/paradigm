/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_version.h"
#include "pdm_config.h"
#include "pdm_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*----------------------------------------------------------------------------
 *
 * Get Version
 *
 * parameters:
 *
 * return:
 *   Version
 *
 *----------------------------------------------------------------------------*/

char *
PDM_version_get
(
 void
)
{
  char *_version;
  PDM_malloc(_version,(strlen(PDM_VERSION) + 1),char);

  strcpy(_version, PDM_VERSION);
  return _version;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
