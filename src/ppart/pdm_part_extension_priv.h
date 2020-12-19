#ifndef __PDM_PART_EXTENSION_PRIV_H__
#define __PDM_PART_EXTENSION_PRIV_H__

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

  int             n_block;
  const int      *n_part;
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_PRIV_H__ */
