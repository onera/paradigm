#ifndef __PDM_PART_EXTENSION_PRIV_H__
#define __PDM_PART_EXTENSION_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

// #include "pdm_multipart.h"
#include "pdm_part_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_part_extension_t
 * \brief  Distributed cube
 *
 * _dcube_t define a distributed mesh of a cube
 *
 */

struct _pdm_part_extension_t {
  PDM_MPI_Comm     comm;            /*!< MPI communicator                          */
  PDM_ownership_t  owner;           /*!< Which have the responsabilities of results*/

  int             n_domain;
  const int      *n_part;

  /* Plus tard --> Un pdm_multipart */
  _part_t  **parts;

};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_PRIV_H__ */
