#ifndef __PDM_HASH_TAB_PRIV_H__
#define __PDM_HASH_TAB_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_error.h"
#include "pdm_hash_tab.h"
#include "pdm_printf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _hash_tab_t
 * \brief  Hash table
 *
 * \ref _hash_tab_t defines a hash table structure
 *
 */
struct _hash_tab_t {

  PDM_hash_tab_key_t    t_key;        /*!< Type of key                  */
  int                  *n_data_key;   /*!< Number of data for each key  */
  PDM_g_num_t           key_max;      /*!< Key max                      */
  void               ***data;         /*!< Data                         */
  int                  *m_data_key;   /*!< Max data for each key        */
  int                   n_key_info;   /*!< Number keys with information */
  int                   l_key_info;   /*!< Size of \ref key_info        */
  PDM_g_num_t          *key_info;     /*!< list of Keys with info       */

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HASH_TAB_PRIV_H__ */
