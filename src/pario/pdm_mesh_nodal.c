
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * \brief Storage of mesh handles
 *
 */

static PDM_Handles_t *mesh_handles = NULL; 


/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 * 
 * \brief Free vtx structure
 *
 * \param[inout]  vtx    Vertices
 *
 * \return        NULL
 * 
 */

static 
PDM_Mesh_nodal_vtx_t *
_vtx_free
(
 PDM_Mesh_nodal_vtx_t *vtx
)
{
  if (vtx->parent != NULL) {
    vtx->parent =_vtx_free (vtx->parent);
  }

  if (vtx->coords != NULL) {
    free (vtx->coords);
    vtx->coords = NULL;
  }

  free (vtx);
  
  return NULL;
}


/**
 * 
 * \brief Initialize a mesh
 *
 * \param [inout]  mesh        Mesh
 * \param [in]     n_part      Number of partitions
 */

static void   
_mesh_init
(
PDM_Mesh_nodal_t *mesh,
const int        n_part
)
{

  mesh->n_som_abs                = 0;          
  mesh->n_elt_abs                = 0;          
  mesh->n_part                   = n_part;             

  mesh->vtx                      = malloc(n_part * sizeof(PDM_Mesh_nodal_vtx_t *));
  for (int i = 0; i < n_part; i++) {
    mesh->vtx[i] = malloc(sizeof(PDM_Mesh_nodal_vtx_t));
    mesh->vtx[i]->_coords = NULL;
    mesh->vtx[i]->_numabs = NULL;
    mesh->vtx[i]->n_vtx   = 0;
  }
  
  mesh->n_cell                   = malloc(n_part * sizeof(int));  

  mesh->blocks_std               = NULL;  
  mesh->blocks_poly2d            = NULL;            
  mesh->blocks_poly3d            = NULL;           

  mesh->pdm_mpi_comm             = PDM_MPI_COMM_NULL;
  mesh->prepa_blocks             = NULL;
  mesh->num_cell_parent_to_local = NULL;

}


/**
 * 
 * \brief Free partially a standard block
 *
 * \param [inout]  _bloc_std    Standard block
 *  
 */

static
void
_block_std_free_partial
(
PDM_Mesh_nodal_block_std_t *_block_std
)
{
  if (_block_std->_connec != NULL) {
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_connec[i] != NULL)
          free(_block_std->_connec[i]);
        _block_std->_connec[i] = NULL;
      }
    }
    free(_block_std->_connec);
    _block_std->_connec = NULL;
  }
  
  if (_block_std->_numabs != NULL) {
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_numabs[i] != NULL)
          free(_block_std->_numabs[i]);
        _block_std->_numabs[i] = NULL;
      }
    }
    free(_block_std->_numabs);
    _block_std->_numabs = NULL;
  }
  
  if (_block_std->_num_part != NULL) {
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_num_part[i] != NULL)
          free(_block_std->_num_part[i]);
        _block_std->_num_part[i] = NULL;
      }
    }
    free(_block_std->_num_part);
    _block_std->_num_part = NULL;
  }
}


/**
 * 
 * \brief Free a standard block
 *
 * \param [inout]  _bloc_std    Standard block
 *
 * \return         Null
 *   
 */

static
PDM_Mesh_nodal_block_std_t *
_block_std_free
(
PDM_Mesh_nodal_block_std_t *_block_std
)
{
  _block_std_free_partial(_block_std);

  if (_block_std->n_elt != NULL) {
    free(_block_std->n_elt);
    _block_std->n_elt = NULL;
  }

  if (_block_std->numabs_int != NULL) {
    for (int j = 0; j < _block_std->n_part; j++) {
      if (_block_std->numabs_int[j] != NULL) {
        free(_block_std->numabs_int[j]);
      }
    }
    free(_block_std->numabs_int);
    _block_std->numabs_int = NULL;
  }

  free(_block_std);
  return NULL;
}


/**
 * 
 * \brief Free partially a polygon block
 *
 * \param [inout]  _block_poly2d   polygon block
 *   
 */

static
void
_block_poly2d_free_partial
(
PDM_Mesh_nodal_block_poly2d_t *_block_poly2d
)
{

  if (_block_poly2d->_connec_idx != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_connec_idx[i] != NULL)
          free(_block_poly2d->_connec_idx[i]);
        _block_poly2d->_connec_idx[i] = NULL;
      }
    }
    free(_block_poly2d->_connec_idx);
    _block_poly2d->_connec_idx = NULL;
  }

  if (_block_poly2d->_connec != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_connec[i] != NULL)
          free(_block_poly2d->_connec[i]);
        _block_poly2d->_connec[i] = NULL;
      }
    }
    free(_block_poly2d->_connec);
    _block_poly2d->_connec = NULL;
  }
  
  if (_block_poly2d->_num_part != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_num_part[i] != NULL)
          free(_block_poly2d->_num_part[i]);
        _block_poly2d->_num_part[i] = NULL;
      }
    }
    free(_block_poly2d->_num_part);
    _block_poly2d->_num_part = NULL;
  }
  
  if (_block_poly2d->_numabs != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_numabs[i] != NULL)
          free(_block_poly2d->_numabs[i]);
        _block_poly2d->_numabs[i] = NULL;
      }
    }
    free(_block_poly2d->_numabs);
    _block_poly2d->_numabs = NULL;
  }

}


/**
 * 
 * \brief Free a polygon block
 *
 * \param [inout]  _bloc_poly2d    Polygon block
 *
 * \return         Null
 *   
 */

static
PDM_Mesh_nodal_block_poly2d_t *
_block_poly2d_free
(
PDM_Mesh_nodal_block_poly2d_t *_block_poly2d
)
{
  _block_poly2d_free_partial(_block_poly2d);
  
  if (_block_poly2d->n_elt != NULL) {
    free(_block_poly2d->n_elt);
    _block_poly2d->n_elt = NULL;
  }

  if (_block_poly2d->numabs_int != NULL) {
    for (int j = 0; j < _block_poly2d->n_part; j++) {
      if (_block_poly2d->numabs_int[j] != NULL) {
        free(_block_poly2d->numabs_int[j]);
      }
    }
    free(_block_poly2d->numabs_int);
    _block_poly2d->numabs_int = NULL;
  }

  free(_block_poly2d);

  return NULL;
}


/**
 * 
 * \brief Free partially a polyhedron block
 *
 * \param [inout]  _block_poly3d   polyhedron block
 *   
 */

static
void
_block_poly3d_free_partial
(
PDM_Mesh_nodal_block_poly3d_t *_block_poly3d
)
{
  
  if (_block_poly3d->_facvtx_idx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_facvtx_idx[i] != NULL)
          free(_block_poly3d->_facvtx_idx[i]);
        _block_poly3d->_facvtx_idx[i] = NULL;
      }
    }
    free(_block_poly3d->_facvtx_idx);
    _block_poly3d->_facvtx_idx = NULL;
  }
  
  if (_block_poly3d->_facvtx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_facvtx[i] != NULL)
          free(_block_poly3d->_facvtx[i]);
        _block_poly3d->_facvtx[i] = NULL;
      }
    }
    free(_block_poly3d->_facvtx);
    _block_poly3d->_facvtx = NULL;
  }
  
  if (_block_poly3d->_cellfac_idx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_cellfac_idx[i] != NULL)
          free(_block_poly3d->_cellfac_idx[i]);
        _block_poly3d->_cellfac_idx[i] = NULL;
      }
    }
    free(_block_poly3d->_cellfac_idx);
    _block_poly3d->_cellfac_idx = NULL;
  }
  
  if (_block_poly3d->_cellfac != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_cellfac[i] != NULL)
          free(_block_poly3d->_cellfac[i]);
        _block_poly3d->_cellfac[i] = NULL;
      }
    }
    free(_block_poly3d->_cellfac);
    _block_poly3d->_cellfac = NULL;
  }
  
  if (_block_poly3d->_numabs != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_numabs[i] != NULL)
          free(_block_poly3d->_numabs[i]);
        _block_poly3d->_numabs[i] = NULL;
      }
    }
    free(_block_poly3d->_numabs);
    _block_poly3d->_numabs = NULL;
  }
}


/**
 * 
 * \brief Free a polyhedron block
 *
 * \param [inout]  _block_poly3d    Polyhedron block
 *
 * \return         Null
 *   
 */

static
PDM_Mesh_nodal_block_poly3d_t *
_block_poly3d_free
(
PDM_Mesh_nodal_block_poly3d_t *_block_poly3d
)
{
  _block_poly3d_free_partial(_block_poly3d);

  if (_block_poly3d->n_elt != NULL) {
    free(_block_poly3d->n_elt);
    _block_poly3d->n_elt = NULL;
  }

  if (_block_poly3d->n_face!= NULL) {
    free(_block_poly3d->n_face);
    _block_poly3d->n_face= NULL;
  }

  if (_block_poly3d->numabs_int != NULL) {
    for (int j = 0; j < _block_poly3d->n_part; j++) {
      if (_block_poly3d->numabs_int[j] != NULL) {
        free(_block_poly3d->numabs_int[j]);
      }
    }
    free(_block_poly3d->numabs_int);
    _block_poly3d->numabs_int = NULL;
  }

  free(_block_poly3d);
  return NULL;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 *
 * \return       New mesh nodal handle
 *
 */

int 
PDM_Mesh_nodal_create
(
const int     n_part
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) malloc (sizeof(PDM_Mesh_nodal_t));
  
  _mesh_init (mesh, n_part);
  
  return PDM_Handles_store (mesh_handles, mesh);
}

/**
 * \brief Free partially a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_Mesh_nodal_partial_free
(
const int idx
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
	
  if (mesh->blocks_std != NULL) {
    const int n_blocks_std = PDM_Handles_n_get (mesh->blocks_std);
    const int *list_ind = PDM_Handles_idx_get (mesh->blocks_std);

    for (int i = 0; i < n_blocks_std; i++) {
      PDM_Mesh_nodal_block_std_t *_block_std = 
              (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_std, list_ind[i]);
      _block_std_free_partial(_block_std);
    }
  }
	
  if (mesh->blocks_poly2d != NULL) {
    const int n_blocks_poly2d = PDM_Handles_n_get (mesh->blocks_poly2d);
    const int *list_ind = PDM_Handles_idx_get (mesh->blocks_poly2d);

    for (int i = 0; i < n_blocks_poly2d; i++) {
      PDM_Mesh_nodal_block_std_t *_block_std = 
              (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_poly2d, list_ind[i]);
      _block_std_free_partial(_block_std);
    }
  }
	
  if (mesh->blocks_poly3d != NULL) {
    const int n_blocks_poly3d = PDM_Handles_n_get (mesh->blocks_poly3d);
    const int *list_ind = PDM_Handles_idx_get (mesh->blocks_poly3d);

    for (int i = 0; i < n_blocks_poly3d; i++) {
      PDM_Mesh_nodal_block_std_t *_block_std = 
              (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_poly3d, list_ind[i]);
      _block_std_free_partial(_block_std);
    }
  }
}


/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_Mesh_nodal_free
(
const int idx
)
{
  
  PDM_Mesh_nodal_partial_free (idx);
  
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh != NULL) {

    /* Free vertices */
  
    if (mesh->vtx != NULL) {
      for (int i = 0; i < mesh->n_part; i++) {
        mesh->vtx[i] = _vtx_free (mesh->vtx[i]);
      }

      free(mesh->vtx);
      mesh->vtx = NULL;
    }

    /* free standard blocks */

    if (mesh->blocks_std != NULL) {
      int n_blocks_std = PDM_Handles_n_get (mesh->blocks_std);
      const int *list_ind = PDM_Handles_idx_get (mesh->blocks_std);

      while (n_blocks_std > 0) {
        PDM_Mesh_nodal_block_std_t *_bloc_std = 
          (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_std, list_ind[0]);
        _block_std_free(_bloc_std);
        PDM_Handles_handle_free (mesh->blocks_std, list_ind[0], PDM_FALSE);
        n_blocks_std = PDM_Handles_n_get (mesh->blocks_std);
      }

      mesh->blocks_std = PDM_Handles_free (mesh->blocks_std); 
    }

    /* Free polygon blocks */ 

    if (mesh->blocks_poly2d != NULL) {
      int n_blocks_poly2d = PDM_Handles_n_get (mesh->blocks_poly2d);
      const int *list_ind = PDM_Handles_idx_get (mesh->blocks_poly2d);

      while (n_blocks_poly2d > 0) {
        PDM_Mesh_nodal_block_poly2d_t *_bloc_poly2d = 
          (PDM_Mesh_nodal_block_poly2d_t *) PDM_Handles_get (mesh->blocks_poly2d, list_ind[0]);
        _block_poly2d_free(_bloc_poly2d);
        PDM_Handles_handle_free (mesh->blocks_poly2d, list_ind[0], PDM_FALSE);
        n_blocks_poly2d = PDM_Handles_n_get (mesh->blocks_poly2d);
      }

      mesh->blocks_poly2d = PDM_Handles_free (mesh->blocks_poly2d); 
    }

    /* Free polyhedron blocks */ 

    if (mesh->blocks_poly3d != NULL) {
      int n_blocks_poly3d = PDM_Handles_n_get (mesh->blocks_poly3d);
      const int *list_ind = PDM_Handles_idx_get (mesh->blocks_poly3d);

      while (n_blocks_poly3d > 0) {
        PDM_Mesh_nodal_block_poly3d_t *_bloc_poly3d = 
          (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, list_ind[0]);
        _block_poly3d_free(_bloc_poly3d);
        PDM_Handles_handle_free (mesh->blocks_poly3d, list_ind[0], PDM_FALSE);
        n_blocks_poly3d = PDM_Handles_n_get (mesh->blocks_poly3d);
      }

      mesh->blocks_poly3d = PDM_Handles_free (mesh->blocks_poly3d); 
    }

    /* Free structure */ 

    if (mesh->num_cell_parent_to_local != NULL) {
      for (int ipart = 0; ipart < mesh->n_part; ipart++) {
        if (mesh->num_cell_parent_to_local[ipart] != NULL)
          free(mesh->num_cell_parent_to_local[ipart]);
      }
      free(mesh->num_cell_parent_to_local);
      mesh->num_cell_parent_to_local = NULL;
    }

    free(mesh->n_cell);
    mesh->n_cell = NULL;

    free(mesh);

    PDM_Handles_handle_free (mesh_handles, idx, PDM_FALSE);
  
    int n_mesh_array = PDM_Handles_n_get (mesh_handles); 

    if (n_mesh_array == 0) {
      mesh_handles = PDM_Handles_free (mesh_handles);
    }
  }
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  numabs    Global numbering
 *
 */

void
PDM_Mesh_nodal_coord_set
(
 const int          idx,
 const int          id_part, 
 const int          n_vtx,  
 const PDM_real_t  *coords,  
 const PDM_g_num_t *numabs
)
{
  PDM_Mesh_nodal_t * mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  if (id_part <= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");  
  } 
  
  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  if ((vtx->_coords != NULL) ||
      (vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "these partition vertices are already defined\n");
  }

  /* Mapping memoire */

  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;
  vtx->_numabs = numabs;

}


/**
 * \brief  Return number of vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Number of vertices
 *
 */

int
PDM_Mesh_nodal_n_vertices_get
(
 PDM_Mesh_nodal_t  *mesh
)
{
  return 0; 
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Coordinates of vertices
 *
 */

double *
PDM_Mesh_nodal_vertices_get
(
 PDM_Mesh_nodal_t  *mesh
)
{
  return NULL;
}


/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_vtx          Number of vertices
 * \param [in]  n_vtx_parent   Number of parent vertices
 * \param [in]  numabs         Global numbering (size = \ref n_vtx)
 * \param [in]  num_parent     Numbering in the parent numbering (size = \ref n_vtx)
 * \param [in]  coords_parent  Parent interlaced coordinates (size = 3 * \ref n_vtx_parent)
 * \param [in]  numabs_parent  Parent global numbering (size = \ref n_vtx_parent)
 *
 */

void
PDM_Mesh_nodal_coord_from_parent_set
(
 PDM_Mesh_nodal_t  *mesh,
 const int          id_part, 
 const int          n_vtx,  
 const int          n_vtx_parent,  
 const PDM_g_num_t *numabs,
 const int         *num_parent,
 const PDM_real_t  *coords_parent,  
 const PDM_g_num_t *numabs_parent
)
{
}


/**
 * \brief  Return number of blocks
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return  Number of blocks
 *
 */

int
PDM_Mesh_n_blocks_get
(
 PDM_Mesh_nodal_t  *mesh
)
{
 return 0;
}


/**
 * \brief  Return type of block
 *
 * \param [in]  mesh       Nodal mesh
 * \param [in]  id_block   Block identifier
 *
 * \return  Type of block
 *
 */

int
PDM_Mesh_block_type_get
(
 PDM_Mesh_nodal_t     *mesh,
 const int            id_bloc     
)
{
 return 0;
}


/**
 * \brief  Add a new block to the current mesh
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_block       Block identifier
 *
 * \return      
 *
 */

int 
PDM_Mesh_nodal_block_add 
(
PDM_Mesh_nodal_t            *mesh,
PDM_bool_t                   st_free_data,  
const PDM_Mesh_nodal_elt_t   t_elt
) 
{
 return 0;
}


/**
 * \brief Define a standard block
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x            
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 * \param [in]  numabs         Global numbering
 *
 */

void
PDM_Mesh_nodal_block_std_set 
(
PDM_Mesh_nodal_t    *mesh,
const int            id_block,     
const int            id_part, 
const int            n_elt,    
      PDM_l_num_t   *connec,   
      PDM_g_num_t   *numabs
) 
{
}


/**
 * \brief Return standard block description
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x            
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out]  n_elt          Number of elements
 * \param [out]  connect        Connectivity
 * \param [out]  numabs         Global numbering
 *
 */

void
PDM_Mesh_nodal_block_std_get 
(
PDM_Mesh_nodal_t    *mesh,
const int            id_block,     
const int            id_part, 
const int           *n_elt,    
      PDM_l_num_t  **connec,   
      PDM_g_num_t  **numabs
) 
{
}


/**
 * \brief Define a polygon block
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [in]  numabs         Global numbering
 *
 */
 
void
PDM_Mesh_nodal_block_poly2d_set 
(
PDM_Mesh_nodal_t    *mesh,
const int            id_block, 
const int            id_part, 
const PDM_l_num_t    n_elt,    
      PDM_l_num_t   *connec_idx,   
      PDM_l_num_t   *connec,
      PDM_g_num_t   *numabs
) 
{
}


/**
 * \brief Return a polygon block description
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] n_elt          Number of elements
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [out] numabs         Global numbering
 *
 */
 
void
PDM_Mesh_nodal_block_poly2d_get 
(
PDM_Mesh_nodal_t    *mesh,
const int            id_block, 
const int            id_part, 
const PDM_l_num_t   *n_elt,    
      PDM_l_num_t  **connec_idx,   
      PDM_l_num_t  **connec,
      PDM_g_num_t  **numabs
) 
{
}


/**
 * \brief Define a polyhedra block
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 * \param [in]  numabs         Global numbering
 *
 */

void
PDM_Mesh_nodal_block_poly3d_set 
(
PDM_Mesh_nodal_t    *mesh,
const int            id_block, 
const int            id_part, 
const PDM_l_num_t    n_elt,    
const PDM_l_num_t    n_face,   
      PDM_l_num_t   *facvtx_idx,   
      PDM_l_num_t   *facvtx,
      PDM_l_num_t   *cellfac_idx,   
      PDM_l_num_t   *cellfac,
      PDM_g_num_t   *numabs
) 
{
}


/**
 * \brief  Add some 3D cells from cell face conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection 
 * between old and new numbering.
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  face_vtx_idx   Index of face vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each face
 * \param [in]  face_vtx       Face vertex connectivity
 * \param [in]  cell_face_idx  Index of cell face connectivity
 * \param [in]  cell_face_nb   Number of faces for each cell
 * \param [in]  cell_face      Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [out] ind_num        old numbering to new numbering
 *
 */

void
PDM_Mesh_nodal_mesh_cell3d_cellface_add
(
PDM_Mesh_nodal_t *mesh,
const int         id_part, 
const int         n_elt,
const int         n_face,
PDM_l_num_t      *face_vtx_idx,
PDM_l_num_t      *face_vtx_nb,
PDM_l_num_t      *face_vtx,
PDM_l_num_t      *cell_face_idx,
PDM_l_num_t      *cell_face_nb,
PDM_l_num_t      *cell_face,
PDM_g_num_t      *numabs,
PDM_l_num_t      *ind_num
) 
{
}


/**
 * \brief  Add some 2D cells from cell edge conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection 
 * between old and new numbering.
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_edge         Number of edges used to describe polyhedra
 * \param [in]  edge_vtx_idx   Index of edge vertex connectivity
 * \param [in]  edge_vtx_nb    Number of vertices for each edge
 * \param [in]  edge_vtx       Edge vertex connectivity
 * \param [in]  cell_edge_idx  Index of cell edge connectivity
 * \param [in]  cell_edge_nb   Number of edges for each cell
 * \param [in]  cell_edge      Cell edge connectivity
 * \param [in]  numabs         Global numbering
 * \param [out] ind_num        old numbering to new numbering
 *
 */

void
PDM_Mesh_nodal_geom_cell2d_celledge_add
(
PDM_Mesh_nodal_t  *mesh,
const int          id_part, 
const int          n_elt,
const int          n_edge,
PDM_l_num_t       *edge_vtx_idx,
PDM_l_num_t       *edge_vtx_nb,
PDM_l_num_t       *edge_vtx,
PDM_l_num_t       *cell_edge_idx,
PDM_l_num_t       *cell_edge_nb,
PDM_l_num_t       *cell_edge,
PDM_g_num_t       *numabs,
PDM_l_num_t       *ind_num
) 
{
}


/**
 * \brief  Add some 2D cells from cell vertex connectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection 
 * between old and new numbering.
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_face         Number of polygon
 * \param [in]  face_vtx_idx   Index of edge vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each edge
 * \param [in]  face_vtx       Edge vertex connectivity
 * \param [in]  numabs         Global numbering
 * \param [out] ind_num        old numbering to new numbering
 *
 */

void
PDM_Mesh_nodal_geom_faces_facesom_add
(
PDM_Mesh_nodal_t *mesh,
const int         id_part, 
const int         n_face,
PDM_l_num_t      *face_vtx_idx,
PDM_l_num_t      *face_vtx_nb,
PDM_l_num_t      *face_vtx,
PDM_g_num_t      *numabs,
PDM_l_num_t      *ind_num
) 
{
}


#ifdef __cplusplus
}
#endif /* __cplusplus */

