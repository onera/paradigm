
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_hash_tab.h"
#include "pdm.h"
#include "pdm_hash_tab_priv.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"
#include "pdm_printf.h"
#include "pdm_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_hash_tab_t *
PDM_hash_tab_create
(
const PDM_hash_tab_key_t  t_key,
void                     *key_max
)
{
  PDM_hash_tab_t *ht;
  PDM_malloc(ht, 1, PDM_hash_tab_t);
  const int n_data_default = 0;

  ht->t_key = t_key;

  if (t_key == PDM_HASH_TAB_KEY_INT) {
    ht->key_max = (PDM_g_num_t) *((int *) key_max);
  } else if (t_key == PDM_HASH_TAB_KEY_LONG) {
    ht->key_max = *((PDM_g_num_t *) key_max);
  } else {
    PDM_error("Error PDM_hash_tab_create : Unknown hey type");
  }

  ht->l_key_info = PDM_MAX (ht->key_max/10, 2);
  PDM_malloc(ht->key_info, ht->l_key_info, PDM_g_num_t);
  ht->n_key_info = 0;

  PDM_malloc(ht->data      , ht->key_max, void **);
  PDM_malloc(ht->n_data_key, ht->key_max, int    );
  PDM_malloc(ht->m_data_key, ht->key_max, int    );
  for (int i = 0; i < ht->key_max; i++) {
    ht->n_data_key[i] = 0;
    ht->m_data_key[i] = n_data_default;
    PDM_malloc(ht->data[i], n_data_default, void *);
    for (int j = 0; j < n_data_default; j++) {
      ht->data[i][j] = NULL;
    }
  }

  return (PDM_hash_tab_t *) ht;
}

void
PDM_hash_tab_data_add
(
PDM_hash_tab_t *ht,
void           *key,
void           *data
)
{
  PDM_g_num_t _key = -1;

  if (ht->t_key == PDM_HASH_TAB_KEY_INT) {
    _key = (PDM_g_num_t) *((int *) (key));
  } else if (ht->t_key == PDM_HASH_TAB_KEY_LONG) {
    _key = *((PDM_g_num_t *) (key));
  } else {
	  PDM_error("PDM_hash_tab_data_add error : unknown PDM_hash_tab_key_t\n");
	}

  assert ((PDM_g_num_t) _key < ht->key_max);

  if (ht->n_data_key[_key] >= ht->m_data_key[_key]) {
    ht->m_data_key[_key] += PDM_MAX (1, ht->m_data_key[_key]);
    PDM_realloc(ht->data[_key], ht->data[_key], ht->m_data_key[_key], void *);
  }

  if (ht->n_data_key[_key] == 0) {
    if (ht->n_key_info >= ht->l_key_info) {
      ht->l_key_info += PDM_MAX (1, ht->l_key_info/3);
      PDM_realloc(ht->key_info, ht->key_info, ht->l_key_info, PDM_g_num_t);
    }
    ht->key_info[ht->n_key_info] = _key;
    ht->n_key_info += 1;
  }

  ht->data[_key][ht->n_data_key[_key]++] = data;

}

void
PDM_hash_tab_data_free
(
PDM_hash_tab_t *ht,
void           *key
)
{
  PDM_g_num_t _key = -1;

  if (ht->t_key == PDM_HASH_TAB_KEY_INT) {
    _key = (PDM_g_num_t) *((int *) (key));
  } else if (ht->t_key == PDM_HASH_TAB_KEY_LONG) {
    _key = *((PDM_g_num_t *) (key));
  } else {
	  PDM_error("PDM_hash_tab_data_free error : unknown PDM_hash_tab_key_t\n");
	}

  assert ((PDM_g_num_t) _key < ht->key_max);

  for (int i = 0; i < ht->n_data_key[_key]; i++) {
    if (ht->data[_key][i] != NULL) {
      PDM_free(ht->data[_key][i]);
    }
    ht->data[_key][i] = NULL;
  }

  ht->n_data_key[_key] = 0;

}

int
PDM_hash_tab_n_data_get
(
PDM_hash_tab_t *ht,
void           *key
)
{

  PDM_g_num_t _key = -1;

  if (ht->t_key == PDM_HASH_TAB_KEY_INT) {
    _key = (PDM_g_num_t) *((int *) (key));
  } else if (ht->t_key == PDM_HASH_TAB_KEY_LONG) {
    _key = *((PDM_g_num_t *) (key));
  } else {
	  PDM_error("PDM_hash_tab_data_get error : unknown PDM_hash_tab_key_t\n");
	}

  assert ((PDM_g_num_t) _key < ht->key_max);

  return ht->n_data_key[_key];
}


void **
PDM_hash_tab_data_get
(
PDM_hash_tab_t *ht,
void           *key
)
{
  PDM_g_num_t _key = -1;

  if (ht->t_key == PDM_HASH_TAB_KEY_INT) {
    _key = (PDM_g_num_t) *((int *) (key));
  } else if (ht->t_key == PDM_HASH_TAB_KEY_LONG) {
    _key = *((PDM_g_num_t *) (key));
  } else {
	  PDM_error("PDM_hash_tab_data_get error : unknown PDM_hash_tab_key_t\n");
	}
  assert ((PDM_g_num_t) _key < ht->key_max);
  return ht->data[_key];
}


PDM_hash_tab_t *
PDM_hash_tab_free
(
PDM_hash_tab_t *ht
)
{
  for (int i = 0; i<ht->key_max; i++) {
    PDM_free(ht->data[i]);
  }
  PDM_free(ht->data);
  PDM_free(ht->n_data_key);
  PDM_free(ht->m_data_key);
  PDM_free(ht->key_info);

  PDM_free(ht);

  return NULL;
}


int
PDM_hash_tab_n_used_keys_get
(
PDM_hash_tab_t *ht
)
{
  return ht->n_key_info;
}


PDM_g_num_t *
PDM_hash_tab_used_keys_get
(
PDM_hash_tab_t *ht
)
{
  return ht->key_info;
}


void
PDM_hash_tab_purge
(
PDM_hash_tab_t *ht,
PDM_bool_t      remove_data
)
{

  if (remove_data) {
    for (int i = 0; i < ht->n_key_info; i++) {
      PDM_g_num_t _key = ht->key_info[i];
      PDM_hash_tab_data_free (ht, &_key);
    }
  }

  int *_n_data_key = ht->n_data_key;
  for (int i = 0; i < ht->n_key_info; i++) {
    PDM_g_num_t _key = ht->key_info[i];
    _n_data_key[_key] = 0;
  }

  ht->n_key_info = 0;
}


void *
PDM_hash_tab_key_max_get
(
PDM_hash_tab_t *ht
)
{
  return &(ht->key_max);
}


PDM_hash_tab_key_t
PDM_hash_tab_key_type_get
(
PDM_hash_tab_t *ht
)
{
  return ht->t_key;
}

void
PDM_hash_tab_dump
(
PDM_hash_tab_t *ht
)
{
 PDM_printf ("==== PDM_hash_tab_dump ==== \n");
  for (int i=0; i <ht->key_max; i++){
	  PDM_printf ("ht->n_data_key[%d] : %d\n", i, ht->n_data_key[i]);
  }
  PDM_printf ("ht->data = %d\n", ht->data);
  for (int key=0; key < ht->key_max; key++){
    PDM_printf ("ht->data[%d] = %d, ",key, ht->data[key]);
    int n_data_in_key = PDM_hash_tab_n_data_get(ht, &key );
    for(int i_data = 0; i_data < n_data_in_key; ++i_data){
      PDM_printf ("ht->data[%d][%i] = %d, ", key, i_data, ht->data[key][i_data]);
    }
    PDM_printf ("\n");
  }

  PDM_printf ("==== PDM_hash_tab_dump ==== terminated ====\n");
}

int
PDM_hash_tab_check_collision
(
 PDM_hash_tab_t *ht,
 const int       value,
 const int       key_max,
 int            *key
)
{
  *key = value % key_max;

  int n_data = PDM_hash_tab_n_data_get (ht, key);

  PDM_g_num_t **data = (PDM_g_num_t **) PDM_hash_tab_data_get (ht, key);
  for (int i = 0; i < n_data; i++) {
    if (*(data[i]) == value) {
      return 1;
    }
  }
  return 0;
}



#ifdef	__cplusplus
}
#endif
